/*
 * Copyright 2016-2019 The Hong Kong University of Science and Technology
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package proteomics;

import ProteomicsLibrary.MassTool;
import ProteomicsLibrary.PrepareSpectrum;
import ProteomicsLibrary.Types.SparseVector;
import org.apache.commons.math3.analysis.function.Exp;
import proteomics.Index.BuildIndex;
import proteomics.PTM.InferPTM;
import proteomics.Search.Search;
import proteomics.Segment.InferSegment;
import proteomics.Spectrum.PreSpectra;
import proteomics.Types.*;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;

import java.math.BigInteger;
import java.nio.ByteBuffer;
import java.nio.charset.Charset;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.locks.ReentrantLock;
import com.google.common.graph.*;

public class PreSearch implements Callable<PreSearch.Entry> {
    private static final int candisNum = 20;
    private int hashbits = 16;
    private final BuildIndex buildIndex;
    private final MassTool massTool;
    private final double ms1Tolerance;
    private final double leftInverseMs1Tolerance;
    private final double rightInverseMs1Tolerance;
    private final int ms1ToleranceUnit;
    private final double ms2Tolerance;
    private final double minPtmMass;
    private final double maxPtmMass;
    private final int localMaxMs2Charge;
    private final Map<String, Peptide0> peptide0Map;
    private final JMzReader spectraParser;
    private final double minClear;
    private final double maxClear;
    private final ReentrantLock lock;
    private final String scanId;
    private final int precursorCharge;
    private final double precursorMass;
    private final InferPTM inferPTM;
    private final PrepareSpectrum preSpectrum;
    private final String sqlPath;
    private final int scanNum;
    private final int precursorScanNo;
    private String truth;
    public int truthHash;

    public PreSearch(int scanNum, BuildIndex buildIndex, MassTool massTool, double ms1Tolerance, double leftInverseMs1Tolerance, double rightInverseMs1Tolerance
            , int ms1ToleranceUnit, double ms2Tolerance, double minPtmMass, double maxPtmMass, int localMaxMs2Charge
            , JMzReader spectraParser, double minClear, double maxClear, ReentrantLock lock, String scanId, int precursorCharge, double precursorMass
            , InferPTM inferPTM, PrepareSpectrum preSpectrum, String sqlPath, int precursorScanNo, String truth, int truthHash) {

        this.buildIndex = buildIndex;
        this.massTool = massTool;
        this.ms1Tolerance = ms1Tolerance;
        this.leftInverseMs1Tolerance = leftInverseMs1Tolerance;
        this.rightInverseMs1Tolerance = rightInverseMs1Tolerance;
        this.ms1ToleranceUnit = ms1ToleranceUnit;
        this.ms2Tolerance = ms2Tolerance;
        this.minPtmMass = minPtmMass;
        this.maxPtmMass = maxPtmMass;
        this.localMaxMs2Charge = localMaxMs2Charge;
        this.spectraParser = spectraParser;
        this.minClear = minClear;
        this.maxClear = maxClear;
        this.lock = lock;
        this.scanId = scanId;
        this.precursorCharge = precursorCharge;
        this.precursorMass = precursorMass;
        this.inferPTM = inferPTM;
        this.preSpectrum = preSpectrum;
        this.sqlPath = sqlPath;
        peptide0Map = buildIndex.getPeptide0Map();
        this.scanNum = scanNum;
        this.precursorScanNo = precursorScanNo;
        this.truth = truth;
        this.truthHash = truthHash;
    }

    @Override
    public Entry call() throws Exception {
        Map<Double, Double> rawPLMap;
        try {
            lock.lock();
            // Reading peak list.
            rawPLMap = spectraParser.getSpectrumById(scanId).getPeakList();
        } finally {
            lock.unlock();
        }
        // preprocess peak list
        TreeMap<Double, Double> plMap = preSpectrum.preSpectrumTopNStyle(rawPLMap, precursorMass, precursorCharge, minClear, maxClear, 15);

        if (plMap.isEmpty()) {
            return null;
        }

        // Coding
        InferSegment inferSegment = buildIndex.getInferSegment();
        TreeMap<Double, Double> finalPlMap = inferSegment.addVirtualPeaks(precursorMass, plMap);
        List<ThreeExpAA> expAaLists = inferSegment.inferSegmentLocationFromSpectrum(precursorMass, finalPlMap, scanNum);


        //check for NC tag
        List<ThreeExpAA> ncTags = new ArrayList<>();
        for (ThreeExpAA tag : expAaLists) {
            if ( tag.getHeadLocation() == MassTool.PROTON ) {
                tag.ncTag = ThreeExpAA.NC.N;
                ncTags.add(tag);
            }
            if ( tag.getHeadLocation() == MassTool.PROTON + massTool.H2O ) {
                tag.ncTag = ThreeExpAA.NC.C;
                ncTags.add(tag);

            }
        }

        Map<Double, Integer> mzIdMap = new HashMap<>();
        int i = 0;
        for (double mz : finalPlMap.keySet()) {
            mzIdMap.put(mz, i);
            i++;
        }
//        for (ThreeExpAA tag : expAaLists) {
//            System.out.println(tag.getPtmFreeAAString() + ","
//                    +mzIdMap.get(tag.getExpAAs()[0].getHeadLocation())+ ","+mzIdMap.get(tag.getExpAAs()[1].getHeadLocation())+ ","+mzIdMap.get(tag.getExpAAs()[2].getHeadLocation())+ ","+mzIdMap.get(tag.getExpAAs()[2].getTailLocation())
//                    + "," + tag.getExpAAs()[0].getHeadLocation()+ ","+tag.getExpAAs()[1].getHeadLocation()+ ","+tag.getExpAAs()[2].getHeadLocation()+ ","+tag.getExpAAs()[2].getTailLocation());
//
//        }
        TreeMap<Integer, Double> nodeMap = new TreeMap<>();
        Map<Integer, Map<Integer, Double>> inEdgeMap = new HashMap<>();
        Map<Integer, Map<Integer, Double>> outEdgeMap = new HashMap<>();
        for (ThreeExpAA tag : expAaLists) {
            for(ExpAA aa : tag.getExpAAs()) {
                if ((aa.getHeadLocation() == MassTool.PROTON + massTool.H2O) && (!"KR".contains(aa.getAA()))) {
                    continue;
                }
                int n1 = mzIdMap.get(aa.getHeadLocation());
                int n2 = mzIdMap.get(aa.getTailLocation());
                double inte1 = aa.getHeadIntensity();
                double inte2 = aa.getTailIntensity();
                nodeMap.put(n1, inte1);
                nodeMap.put(n2, inte2);
                if (outEdgeMap.containsKey(n1)) {
                    outEdgeMap.get(n1).put(n2, 0.5*(inte1+inte2));
                } else {
                    Map<Integer, Double> temp = new HashMap<>();
                    temp.put(n2, 0.5*(inte1+inte2));
                    outEdgeMap.put(n1, temp);
                }

                if (inEdgeMap.containsKey(n2)) {
                    inEdgeMap.get(n2).put(n1, 0.5*(inte1+inte2));
                } else {
                    Map<Integer, Double> temp = new HashMap<>();
                    temp.put(n1, 0.5*(inte1+inte2));
                    inEdgeMap.put(n2, temp);
                }

            }
        }
        for (int node : nodeMap.keySet()) {
            if (!inEdgeMap.containsKey(node) ){
                for (int n2 : outEdgeMap.get(node).keySet()) {
                    outEdgeMap.get(node).put(n2, outEdgeMap.get(node).get(n2)+nodeMap.get(node)/2);
                    inEdgeMap.get(n2).put(node, inEdgeMap.get(n2).get(node)+nodeMap.get(node)/2);
                }
            }
            if (!outEdgeMap.containsKey(node) ){
                for (int n1 : inEdgeMap.get(node).keySet()) {
                    outEdgeMap.get(n1).put(node, outEdgeMap.get(n1).get(node)+nodeMap.get(node)/2);
                    inEdgeMap.get(node).put(n1, inEdgeMap.get(node).get(n1)+nodeMap.get(node)/2);
                }
            }
        }
        // finish prepare graph

        Set<Edge> edgeToDel = new HashSet<>();
        Map<Integer, Integer> tempTodel = new HashMap<>();

        for (int node : nodeMap.descendingKeySet()){
            if (outEdgeMap.containsKey(node)) {
                Map<Integer, Double> tempOutMap = outEdgeMap.get(node);
                double maxOut = 0d;
                int maxOutPeak = -1;
                for (int n2 : tempOutMap.keySet()) {
                    if (tempOutMap.get(n2) > maxOut) {
                        maxOut = tempOutMap.get(n2);
                        maxOutPeak = n2;
                    }
                }
                for (int n2 : tempOutMap.keySet()) {
                    if (n2 != maxOutPeak) {
                        edgeToDel.add(new Edge(node,n2));
                        tempTodel.put(node,n2);

                    }
                }
                if (inEdgeMap.containsKey(node)) {

                    for (int n1 : inEdgeMap.get(node).keySet()) {
                        outEdgeMap.get(n1).put(node, outEdgeMap.get(n1).get(node) + maxOut);
                        inEdgeMap.get(node).put(n1, inEdgeMap.get(node).get(n1) + maxOut);

                    }
                }
            }
        }

        for (int n1 : tempTodel.keySet()) {
            int n2 = tempTodel.get(n1);
            inEdgeMap.get(n2).remove(n1);
            outEdgeMap.get(n1).remove(n2);
        }
        for (int node : nodeMap.keySet()){
            if (inEdgeMap.containsKey(node)) {
                Map<Integer, Double> tempInMap = inEdgeMap.get(node);
                double maxIn = 0d;
                int maxInPeak = -1;
                for (int n1 : tempInMap.keySet()) {
                    if (tempInMap.get(n1) > maxIn) {
                        maxIn = tempInMap.get(n1);
                        maxInPeak = n1;
                    }
                }
                for (int n1 : tempInMap.keySet()) {
                    if (n1 != maxInPeak) {
                        edgeToDel.add(new Edge(n1,node));

                    }
                }
                if (outEdgeMap.containsKey(node)) {

                    for (int n2 : outEdgeMap.get(node).keySet()) {
                        outEdgeMap.get(node).put(n2, outEdgeMap.get(node).get(n2) + maxIn);
                        inEdgeMap.get(n2).put(node, inEdgeMap.get(n2).get(node) + maxIn);

                    }
                }
            }
        }
//        truth = "#RFPAEDEFPD#SAHNNHMAK";
        int numCorrectTagsBefore = 0;
        int numTotalTagsAfter = 0;
        int numCorrectTagsAfter = 0;
        List<ThreeExpAA> denoisedTags = new ArrayList<>();
        for(ThreeExpAA tag: expAaLists) {
            String tagSeq = tag.getPtmFreeAAString();
            String revTagSeq = new StringBuilder(tagSeq).reverse().toString();

            if (truth.contains(tagSeq) || truth.contains(revTagSeq)) {
                numCorrectTagsBefore++;
            }

            boolean containsBadEdge = false;
            for (ExpAA aa : tag.getExpAAs()) {
                Edge tempEdge = new Edge(mzIdMap.get(aa.getHeadLocation()), mzIdMap.get(aa.getTailLocation()));
                if (edgeToDel.contains(tempEdge)) {
                    containsBadEdge = true;
                }
            }
            if (!containsBadEdge) {
                denoisedTags.add(tag);
                numTotalTagsAfter++;
                if (truth.contains(tagSeq) || truth.contains(revTagSeq)) {
                    numCorrectTagsAfter++;
                }
            }
        }
//        System.out.println("scan, "+scanNum+","+expAaLists.size() + "," + numCorrectTagsBefore +","+numTotalTagsAfter+","+numCorrectTagsAfter);


        if (!denoisedTags.isEmpty()) {
            Map<String, Double> tagWeightMap = new HashMap<>();
            SparseVector scanCode = inferSegment.generateSegmentIntensityVector(denoisedTags, tagWeightMap);
            if (scanNum == 1886) {
                int a =1 ;
            }
//            tagWeightMap.clear();
//            tagWeightMap.put("GES",1.0);//
//            tagWeightMap.put("HAV",1.0);//
//            tagWeightMap.put("EGT",1.0);//
//            tagWeightMap.put("ESV",1.0);//
//            tagWeightMap.put("ETK",1.14);
//            tagWeightMap.put("AVS",1.0);//
////            tagWeightMap.put("TKV",1.93);
//            tagWeightMap.put("GTK",1.0);//


            int scanHash = simhash32(tagWeightMap);

            Entry entry = new Entry();
            Search search = new Search(truthHash, truth, entry, scanHash, scanNum, buildIndex, precursorMass, scanCode, massTool, ms1Tolerance, leftInverseMs1Tolerance, rightInverseMs1Tolerance, ms1ToleranceUnit, minPtmMass, maxPtmMass, localMaxMs2Charge, ncTags);

//            Search search = new Search(entry, scanNum, buildIndex, precursorMass, scanCode, massTool, ms1Tolerance, leftInverseMs1Tolerance, rightInverseMs1Tolerance, ms1ToleranceUnit, minPtmMass, maxPtmMass, localMaxMs2Charge, ncTags);
//            entry.candidateList = search.candidatesList;
            //move search to here as presearch.

            return entry;
        } else {
            return null;
        }
    }

    public static long hash64(String doc) {
        byte[] buffer = doc.getBytes(Charset.forName("utf-8"));
        ByteBuffer data = ByteBuffer.wrap(buffer);
        return hash64(data, 0, buffer.length, 0);
    }

    public static long hash64(ByteBuffer key, int offset, int length, long seed) {
        long m64 = 0xc6a4a7935bd1e995L;
        int r64 = 47;

        long h64 = (seed & 0xffffffffL) ^ (m64 * length);

        int lenLongs = length >> 3;

        for (int i = 0; i < lenLongs; ++i) {
            int i_8 = i << 3;

            long k64 = ((long) key.get(offset + i_8 + 0) & 0xff)
                    + (((long) key.get(offset + i_8 + 1) & 0xff) << 8)
                    + (((long) key.get(offset + i_8 + 2) & 0xff) << 16)
                    + (((long) key.get(offset + i_8 + 3) & 0xff) << 24)
                    + (((long) key.get(offset + i_8 + 4) & 0xff) << 32)
                    + (((long) key.get(offset + i_8 + 5) & 0xff) << 40)
                    + (((long) key.get(offset + i_8 + 6) & 0xff) << 48)
                    + (((long) key.get(offset + i_8 + 7) & 0xff) << 56);

            k64 *= m64;
            k64 ^= k64 >>> r64;
            k64 *= m64;

            h64 ^= k64;
            h64 *= m64;
        }

        int rem = length & 0x7;

        switch (rem) {
            case 0:
                break;
            case 7:
                h64 ^= (long) key.get(offset + length - rem + 6) << 48;
            case 6:
                h64 ^= (long) key.get(offset + length - rem + 5) << 40;
            case 5:
                h64 ^= (long) key.get(offset + length - rem + 4) << 32;
            case 4:
                h64 ^= (long) key.get(offset + length - rem + 3) << 24;
            case 3:
                h64 ^= (long) key.get(offset + length - rem + 2) << 16;
            case 2:
                h64 ^= (long) key.get(offset + length - rem + 1) << 8;
            case 1:
                h64 ^= (long) key.get(offset + length - rem);
                h64 *= m64;
        }

        h64 ^= h64 >>> r64;
        h64 *= m64;
        h64 ^= h64 >>> r64;

        return h64;
    }
    public static long simhash64(Map<String, Double> tfMap) {
        int bitLen = 64;
        double[] bits = new double[bitLen];
        for (String tag : tfMap.keySet()) {
            double weight = tfMap.get(tag);
            long v = hash64(tag);
            for (int i = bitLen; i >= 1; --i) {
                if (((v >> (bitLen - i)) & 1) == 1)
                    bits[i - 1] += weight;
                else
                    bits[i - 1] -= weight;
            }
        }
        long hash = 0x0000000000000000;
        long one = 0x0000000000000001;
        for (int i = bitLen; i >= 1; --i) {
            if (bits[i - 1] > 0) {
                hash |= one;
            }
            one = one << 1;
        }
        return hash;
    }

    public int simhash32(Map<String, Double> tfMap) {
        int bitLen = 32;
        double[] bits = new double[bitLen];
        for (String tag : tfMap.keySet()) {
            double weight = tfMap.get(tag);
            int v = hash32(tag);
            for (int i = bitLen; i >= 1; --i) {
                if (((v >> (bitLen - i)) & 1) == 1)
                    bits[i - 1] += weight;
                else
                    bits[i - 1] -= weight;
            }
        }
        int hash = 0x00000000;
        int one = 0x00000001;
        for (int i = bitLen; i >= 1; --i) {
            if (bits[i - 1] > 1) {
                hash |= one;
            }
            one = one << 1;
        }
        return hash;
    }

    public static int hash32(String doc) {
        byte[] buffer = doc.getBytes(Charset.forName("utf-8"));
        ByteBuffer data = ByteBuffer.wrap(buffer);
        return hash32(data, 0, buffer.length, 0);
    }

    public static int hash32(ByteBuffer data, int offset, int length, int seed) {
        int m = 0x5bd1e995;
        int r = 24;

        int h = seed ^ length;

        int len_4 = length >> 2;

        for (int i = 0; i < len_4; i++) {
            int i_4 = i << 2;
            int k = data.get(offset + i_4 + 3);
            k = k << 8;
            k = k | (data.get(offset + i_4 + 2) & 0xff);
            k = k << 8;
            k = k | (data.get(offset + i_4 + 1) & 0xff);
            k = k << 8;
            k = k | (data.get(offset + i_4 + 0) & 0xff);
            k *= m;
            k ^= k >>> r;
            k *= m;
            h *= m;
            h ^= k;
        }

        // avoid calculating modulo
        int len_m = len_4 << 2;
        int left = length - len_m;

        if (left != 0) {
            if (left >= 3) {
                h ^= (int) data.get(offset + length - 3) << 16;
            }
            if (left >= 2) {
                h ^= (int) data.get(offset + length - 2) << 8;
            }
            if (left >= 1) {
                h ^= (int) data.get(offset + length - 1);
            }

            h *= m;
        }

        h ^= h >>> 13;
        h *= m;
        h ^= h >>> 15;

        return h;
    }
    private BigInteger simHash(Map<String, Double> tagWeightMap) {
        double[] v = new double[hashbits];
        for (String tag : tagWeightMap.keySet()) {
            double weight = tagWeightMap.get(tag);
            BigInteger t = this.hash(tag);
            for (int i = 0; i < hashbits; i++) {
                BigInteger bitmask = new BigInteger("1").shiftLeft(i);
                if (t.and(bitmask).signum() != 0) {
                    v[i] += weight;
                } else {
                    v[i] -= weight;
                }
            }
        }

        //binString to BigInteger
        BigInteger fingerprint = new BigInteger("0");
        for (int i = 0; i < hashbits; i++) {
            if (v[i] >= 0) {
                fingerprint = fingerprint.add(new BigInteger("1").shiftLeft(i));
            }
        }
        return fingerprint;
    }

    private BigInteger hash(String source) {
        if (source == null || source.length() == 0) {
            return new BigInteger("0");
        } else {
            /**
             * 当sourece 的长度过短，会导致hash算法失效，因此需要对过短的词补偿
             */
            while (source.length() < 3) {
                source = source + source.charAt(0);
            }
            char[] sourceArray = source.toCharArray();
            BigInteger x = BigInteger.valueOf(((long) sourceArray[0]) << 7);
            BigInteger m = new BigInteger("1000003");
            BigInteger mask = new BigInteger("2").pow(hashbits).subtract(new BigInteger("1"));
            for (char item : sourceArray) {
                BigInteger temp = BigInteger.valueOf((long) item);
                x = x.multiply(m).xor(temp).and(mask);
            }
            x = x.xor(new BigInteger(String.valueOf(source.length())));
            if (x.equals(new BigInteger("-1"))) {
                x = new BigInteger("-2");
            }
            return x;
        }
    }
    public class Edge {
        public int n1;
        public int n2;
        public double weight = 0d;
        Edge(int n1, int n2){
            this.n1 = n1;
            this.n2 = n2;
        };

        @Override
        public int hashCode(){
            return n1*n2;
        }
        @Override
        public boolean equals(Object other) {
            if (this.n1 == 31 && this.n2 == 36) {
                int a = 1;
            }
            Edge temp = (Edge) other;
            return (this.n1 == temp.n1) && (this.n2 == temp.n2);
        }
    }
    public class Entry {

        public int scanNum = PreSearch.this.scanNum;
        public double precursorMass = PreSearch.this.precursorMass;
        public List<Peptide> ptmOnlyList = new ArrayList<>();
        public List<Peptide> ptmFreeList = new ArrayList<>();
        public String candidateList = null;
        Entry() {
//            this.scanNum = scanNum;
        }
//        Entry(List<Peptide> ptmOnlyList, List<Peptide> ptmFreeList) {
//
//            this.ptmOnlyList = ptmOnlyList;
//            this.ptmFreeList = ptmFreeList;
//        }
    }
}
