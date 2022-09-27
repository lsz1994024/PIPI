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
import ProteomicsLibrary.SpecProcessor;
import ProteomicsLibrary.Types.SparseVector;
import org.apache.commons.math3.util.Pair;
import proteomics.Index.BuildIndex;
import proteomics.Search.Search;
import proteomics.Segment.InferSegment;
import proteomics.Spectrum.DatasetReader;
import proteomics.Types.*;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;

import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.locks.ReentrantLock;

public class PreSearch implements Callable<PreSearch.Entry> {
    private static final int candisNum = 20;
    private final BuildIndex buildIndex;
    private final MassTool massTool;
    private final double ms1Tolerance;
    private final double leftInverseMs1Tolerance;
    private final double rightInverseMs1Tolerance;
    private final int ms1ToleranceUnit;
    private final double minPtmMass;
    private final double maxPtmMass;
    private final int localMaxMs2Charge;
    private final JMzReader spectraParser;
    private final double minClear;
    private final double maxClear;
    private final ReentrantLock lock;
    private final String scanName;
    private final int precursorCharge;
    private final double precursorMass;
    private final SpecProcessor specProcessor;
    private final int scanNum;
    private String truth;


    public PreSearch(int scanNum, BuildIndex buildIndex, MassTool massTool, double ms1Tolerance, double leftInverseMs1Tolerance, double rightInverseMs1Tolerance
            , int ms1ToleranceUnit, double minPtmMass, double maxPtmMass, int localMaxMs2Charge
            , JMzReader spectraParser, double minClear, double maxClear, ReentrantLock lock, String scanName, int precursorCharge, double precursorMass
            , SpecProcessor specProcessor, String truth) {

        this.buildIndex = buildIndex;
        this.massTool = massTool;
        this.ms1Tolerance = ms1Tolerance;
        this.leftInverseMs1Tolerance = leftInverseMs1Tolerance;
        this.rightInverseMs1Tolerance = rightInverseMs1Tolerance;
        this.ms1ToleranceUnit = ms1ToleranceUnit;
        this.minPtmMass = minPtmMass;
        this.maxPtmMass = maxPtmMass;
        this.localMaxMs2Charge = localMaxMs2Charge;
        this.spectraParser = spectraParser;
        this.minClear = minClear;
        this.maxClear = maxClear;
        this.lock = lock;
        this.scanName = scanName;
        this.precursorCharge = precursorCharge;
        this.precursorMass = precursorMass;
        this.specProcessor = specProcessor;
        this.scanNum = scanNum;
        this.truth = truth;
    }

    @Override
    public Entry call() throws Exception {
        Map<Double, Double> rawPLMap;
        try {
            lock.lock();
            rawPLMap = spectraParser.getSpectrumById(scanName.split("\\.")[1]).getPeakList();
        } finally {
            lock.unlock();
        }
        if (this.scanNum == 1976 || this.scanNum == 1977) {
            int a = 1;
        }
        TreeMap<Double, Double> plMap = specProcessor.preSpectrumTopNStyleWithChargeLimit(rawPLMap, precursorMass, precursorCharge, minClear, maxClear, DatasetReader.topN);

        if (plMap.isEmpty()) {
            return null;
        }


        // Coding
        InferSegment inferSegment = buildIndex.getInferSegment();
        TreeMap<Double, Double> finalPlMap = inferSegment.addVirtualPeaks(precursorMass, plMap);

        List<ThreeExpAA> tag3List = inferSegment.inferSegmentLocationFromSpectrum(precursorMass, finalPlMap, scanNum);

        if (!tag3List.isEmpty()) {
            Entry entry = new Entry();
            //check for NC tag
            List<ThreeExpAA> ncTags = new ArrayList<>();

            Map<Double, Integer> mzIdMap = new HashMap<>();
            int i = 0;
            for (double mz : finalPlMap.keySet()) {
                mzIdMap.put(mz, i);
                i++;
            }
            TreeMap<Integer, Double> nodeMap = new TreeMap<>();
            Map<Integer, Map<Integer, Double>> inEdgeMap = new HashMap<>();
            Map<Integer, Map<Integer, Double>> outEdgeMap = new HashMap<>();
            for (ThreeExpAA tag : tag3List) {
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
            for(ThreeExpAA tag: tag3List) {
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

            SparseVector scanCode = inferSegment.generateSegmentIntensityVector(denoisedTags);

//            Search search = new Search(entry, scanNum, buildIndex, precursorMass, scanCode, massTool, ms1Tolerance, leftInverseMs1Tolerance, rightInverseMs1Tolerance
//                    , ms1ToleranceUnit, minPtmMass, maxPtmMass, localMaxMs2Charge, candiSet, "n"+truth+"c");
            Search search = new Search(entry, scanName, buildIndex, precursorMass, scanCode, massTool, ms1Tolerance, leftInverseMs1Tolerance, rightInverseMs1Tolerance
                    , ms1ToleranceUnit, minPtmMass, maxPtmMass, localMaxMs2Charge, ncTags);



            entry.scanName = this.scanName;

            return entry;
        } else {
//            System.out.println("nullLsz, "+ scanNum);entry = {PreSearch$Entry@2739}
            return null;
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
        public int precursorCharge = PreSearch.this.precursorCharge;
        public List<Peptide> ptmOnlyList = new ArrayList<>();
        public List<Peptide> ptmFreeList = new ArrayList<>();

        public String scanName;
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
