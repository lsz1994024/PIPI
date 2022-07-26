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
import it.unimi.dsi.fastutil.doubles.DoubleLinkedOpenCustomHashSet;
import org.apache.commons.math3.util.Pair;
import proteomics.FM.FMIndex;
import proteomics.FM.SearchInterval;
import proteomics.Index.BuildIndex;
import proteomics.PTM.InferPTM;
import proteomics.Segment.InferSegment;
import proteomics.Spectrum.PreSpectra;
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


    public PreSearch(int scanNum, BuildIndex buildIndex, MassTool massTool, double ms1Tolerance, double leftInverseMs1Tolerance, double rightInverseMs1Tolerance
            , int ms1ToleranceUnit, double ms2Tolerance, double minPtmMass, double maxPtmMass, int localMaxMs2Charge
            , JMzReader spectraParser, double minClear, double maxClear, ReentrantLock lock, String scanId, int precursorCharge, double precursorMass
            , InferPTM inferPTM, PrepareSpectrum preSpectrum, String sqlPath, int precursorScanNo, String truth) {

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
        TreeMap<Double, Double> plMap = preSpectrum.preSpectrumTopNStyle(rawPLMap, precursorMass, precursorCharge, minClear, maxClear, PreSpectra.topN);

        if (plMap.isEmpty()) {
            return null;
        }

        // Coding
        InferSegment inferSegment = buildIndex.getInferSegment();
        TreeMap<Double, Double> finalPlMap = inferSegment.addVirtualPeaks(precursorMass, plMap);

        List<ThreeExpAA> tag3List = inferSegment.inferSegmentLocationFromSpectrum(precursorMass, finalPlMap, scanNum);
        List<ThreeExpAA> tag4List = inferSegment.getAllTag4(precursorMass, finalPlMap, scanNum);
        List<ThreeExpAA> allLongTagList = inferSegment.getLongTag(finalPlMap, precursorMass - massTool.H2O + MassTool.PROTON, scanNum);
        tag4List.sort(Comparator.comparingDouble(ThreeExpAA::getTotalIntensity).reversed());

        if (allLongTagList.isEmpty()) {
            return null;
        }

        //which top long tags to uses
        double scoreCutOff = allLongTagList.get(0).getTotalIntensity()-1;
        int rankCutOff = 10;
        int validTagNum = 0;
        for (int ii = 0; ii < Math.min(rankCutOff, allLongTagList.size()); ii++){
            if (allLongTagList.get(ii).getTotalIntensity() > scoreCutOff) {
                validTagNum++;
            } else {
                break;
            }
        }


//        tag4List.sort(Comparator.comparingDouble(ThreeExpAA::getTotalIntensity).reversed());
        Set<ThreeExpAA> tag8Set = new HashSet<>();
        Set<ThreeExpAA> tag7Set = new HashSet<>();
        Set<ThreeExpAA> tag6Set = new HashSet<>();
//        for (ThreeExpAA longTag : allLongTagList) {
//            if (longTag.size() >= 7) {
//                for (int i = 0; i <= longTag.size() - 7; i++) {
////                    List<ExpAA> expAAList = longTag.expAAList.subList(i, i+5)
//                    tag7Set.add(new ThreeExpAA(longTag.expAAList.subList(i, i+7)));
//                }
//            }
//            if (longTag.size() >= 6) {
//                for (int i = 0; i <= longTag.size() - 6; i++) {
//                    tag6Set.add(new ThreeExpAA(longTag.expAAList.subList(i, i+6)));
//                }
//            }
//            if (longTag.size() >= 8) {
//                for (int i = 0; i <= longTag.size() - 8; i++) {
//                    tag8Set.add(new ThreeExpAA(longTag.expAAList.subList(i, i+8)));
//                }
//            }
//        }
        List<ThreeExpAA> tag7List = new LinkedList<>(tag7Set);
        List<ThreeExpAA> tag8List = new LinkedList<>(tag8Set);
        List<ThreeExpAA> tag6List = new LinkedList<>(tag6Set);
        //check for NC tag
        List<ThreeExpAA> ncTags = new ArrayList<>();
        for (ThreeExpAA tag : tag4List) {
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
        for (ThreeExpAA tag : allLongTagList) {
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
        Map<Integer, Integer> tempToDel = new HashMap<>();

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
                        tempToDel.put(node,n2);

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

        for (int n1 : tempToDel.keySet()) {
            int n2 = tempToDel.get(n1);
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
        Set<Pair<Integer,Integer>> pairToDel = new HashSet<>();
        for (Edge e : edgeToDel){
            pairToDel.add(new Pair(e.n1, e.n2));
        }
        List<ThreeExpAA> allLongDenoisedTagList = inferSegment.getLongDenoisedTag(finalPlMap, precursorMass - massTool.H2O + MassTool.PROTON, scanNum, pairToDel);
//        truth = "#RFPAEDEFPD#SAHNNHMAK";
        Set<String> correctTagsBefore = new HashSet<>();
        Set<String> totalTagsAfter = new HashSet<>();
        Set<String> correctTagsAfter = new HashSet<>();
        List<ThreeExpAA> denoisedTags = new ArrayList<>();
        for(ThreeExpAA tag: allLongTagList) {
            String tagSeq = tag.getPtmFreeAAString();
            String revTagSeq = new StringBuilder(tagSeq).reverse().toString();

            if (truth.contains(tagSeq) || truth.contains(revTagSeq)) {
                correctTagsBefore.add(tagSeq);
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
                totalTagsAfter.add(tag.getPtmFreeAAString());
                if (truth.contains(tagSeq) || truth.contains(revTagSeq)) {
                    correctTagsAfter.add(tag.getPtmFreeAAString());
                }
            }
        }
        if (scanNum == 1886) {
            int a = 1;
        }
        if (!allLongTagList.isEmpty()) {

            int[] debugScans = new int[]{3452,45474,45515,45532,45741,45780,46000,46034,50377,50650,52158,51975,52021,52031,52090,52180,52338,52335,52411,52438,53487};
            List<Integer> debugList = new LinkedList<>();
            for (int deScanNum : debugScans){
                debugList.add(deScanNum);
            }
            if (debugList.contains(scanNum)){
                int a = 1;
            }
            List<Pair<String, Double>> tagSeqList = new ArrayList<>();
            FMIndex fmIndex = buildIndex.fmIndex;

            ArrayList<Pair<Integer, Double>> posScoreList = new ArrayList<>();
            Set<String> usedLongTags = new HashSet<>();
            Map<Integer , String> posTagMap = new HashMap<>();
            Queue<Pair<String, Double>> tagQueue = new LinkedList<>();

            for (ThreeExpAA longTag : allLongTagList.subList(0, Math.min(validTagNum, allLongTagList.size()))) {
                if (longTag.size() < 6) continue;

                double[] intensityArray = longTag.intensityArray;
                String oriTagStr = longTag.getPtmFreeAAString().replace('#', 'L');
//                Queue<Pair<String, Double>> tagQueue = new LinkedList<>();

                tagQueue.offer(new Pair<>(oriTagStr, longTag.getTotalIntensity()));
//                tagList.add(new Pair<>(oriTagStr, longTag.getTotalIntensity()));
            }
//tagQueue.
            LinkedList<Pair<String, Double>> tagList = new LinkedList<>();
            Set<String> tagsThisRound = new HashSet<>();
            while (!tagQueue.isEmpty()) {
                Pair<String, Double> tagPair = tagQueue.poll();
                int ptnForwardCount = 0;
                int ptnBackwardCount = 0;

                String tagStr = tagPair.getFirst();
                if (!usedLongTags.contains(tagStr)) {
                    usedLongTags.add(tagStr);
                    char[] tagChar = tagStr.toCharArray();
                    SearchInterval searchForward = fmIndex.fmSearch(tagChar);
                    if (searchForward != null) {
                        ptnForwardCount = searchForward.ep - searchForward.sp + 1;
                        if (ptnForwardCount > 0) {
                            for (int ii = searchForward.sp; ii <= searchForward.ep; ii++) {
                                posScoreList.add(new Pair<>(fmIndex.SA[ii], tagPair.getSecond()/1));
                                posTagMap.put(fmIndex.SA[ii], tagStr);
                            }
                        }
                    }
                }
                String revTagStr = new StringBuilder(tagStr).reverse().toString();
                if (!usedLongTags.contains(revTagStr)) {
                    usedLongTags.add(revTagStr);
                    char[] revTagChar = revTagStr.toCharArray();
                    SearchInterval searchBackward = fmIndex.fmSearch(revTagChar);
                    if (searchBackward != null) {
                        ptnBackwardCount = searchBackward.ep - searchBackward.sp + 1;
                        if (ptnBackwardCount > 0) {
                            for (int ii = searchBackward.sp; ii <= searchBackward.ep; ii++) {
                                posScoreList.add(new Pair<>(fmIndex.SA[ii], tagPair.getSecond()/1));
                                posTagMap.put(fmIndex.SA[ii], revTagStr);
                            }
                        }
                    }
                }


                if (ptnForwardCount + ptnBackwardCount >= 1) {
                    break; // if it is processing original tags, dont break
                }
                if (tagStr.length() > 6) tagsThisRound.add(tagStr);

                if (tagQueue.isEmpty()) {
                    for (String tag : tagsThisRound) {
                        String subTagL = tag.substring(0, tag.length() - 1);
                        String subTagR = tag.substring(1);

                        if (!usedLongTags.contains(subTagL)){
                            tagQueue.offer(new Pair<>(subTagL, subTagL.length()-0.5));
                        }
                        if (!usedLongTags.contains(subTagR)){
                            tagQueue.offer(new Pair<>(subTagR, subTagR.length()-0.5));
                        }
                    }
                }


            }

            double highestScore = 0;
            Map<String, Double> protScoreMap = new HashMap<>();
            Map<String, Set<String>> protTagMap = new HashMap<>();
            Map<String, Integer> protTimesMap = new HashMap<>();
            for (Pair<Integer, Double> posScore : posScoreList) {
                int res = Arrays.binarySearch(buildIndex.dotPosArr, posScore.getFirst());
                if (res > 0) {
                    System.out.println("located on dot sign," + scanNum);
                }
                String prot = buildIndex.posProtMap.get(-res-2);
                if (protScoreMap.containsKey(prot) ) {
                    protScoreMap.put(prot, protScoreMap.get(prot)+posScore.getSecond());
                    protTimesMap.put(prot, protTimesMap.get(prot)+1);
                    protTagMap.get(prot).add(posTagMap.get(posScore.getFirst()));
                    if (protScoreMap.get(prot)+posScore.getSecond() > highestScore) highestScore = protScoreMap.get(prot)+posScore.getSecond();
//                    if (prot.contentEquals("sp|Q8WXG9|AGRV1_HUMAN")){
//                        System.out.println(scanNum + "," + posScore.getSecond() );
//                    }
                } else {
                    protScoreMap.put(prot, posScore.getSecond());
                    protTimesMap.put(prot, 1);
                    Set<String> temp = new HashSet<>();
                    temp.add(posTagMap.get(posScore.getFirst()));
                    protTagMap.put(prot,temp);
                    if (posScore.getSecond() > highestScore) highestScore = posScore.getSecond();
//                    if (prot.contentEquals("sp|Q8WXG9|AGRV1_HUMAN")){
//                        System.out.println(scanNum + "," + posScore.getSecond() );
//                    }
                }
            }


            Entry entry = new Entry();
            entry.candidateList = tagSeqList;
            Map<String, Double> top1Pep = new HashMap<>();
            for (String prot : protScoreMap.keySet()) {
                if (protScoreMap.get(prot) > highestScore*0.9) {
                    top1Pep.put(prot, protScoreMap.get(prot));
                }
            }
            if (top1Pep.containsKey("sp|Q8WXG9|AGRV1_HUMAN")){
//                System.out.println(scanNum + "," + top1Pep.get("sp|Q8WXG9|AGRV1_HUMAN") );
            }
            entry.protScoreMap = top1Pep;
//            Search search = new Search(entry, scanNum, buildIndex, precursorMass, scanCode, massTool, ms1Tolerance, leftInverseMs1Tolerance, rightInverseMs1Tolerance
//                    , ms1ToleranceUnit, minPtmMass, maxPtmMass, localMaxMs2Charge, candiSet, "n"+truth+"c");


            return entry;
        } else {
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
        public List<Peptide> ptmOnlyList = new ArrayList<>();
        public List<Peptide> ptmFreeList = new ArrayList<>();
        public List<Pair<String, Double>> candidateList = new ArrayList<>();
        public Map<String, Double> protScoreMap = new HashMap<>();
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
