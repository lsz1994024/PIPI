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
import proteomics.FM.FMIndex;
import proteomics.FM.SearchInterval;
import proteomics.Index.BuildIndex;
import proteomics.PTM.InferPTM;
import proteomics.Segment.InferSegment;
import proteomics.Spectrum.DatasetReader;
import proteomics.Types.Peptide0;
import proteomics.Types.ThreeExpAA;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;

import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.locks.ReentrantLock;

public class SpecCoder implements Callable<SpecCoder.Entry> {
    private static final int candisNum = 20;
    private final BuildIndex buildIndex;
    private final MassTool massTool;
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


    public SpecCoder(int scanNum, BuildIndex buildIndex, MassTool massTool, JMzReader spectraParser, double minClear, double maxClear, ReentrantLock lock
            , String scanName, int precursorCharge, double precursorMass, SpecProcessor specProcessor, String truth) {

        this.buildIndex = buildIndex;
        this.massTool = massTool;
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
            rawPLMap = spectraParser.getSpectrumById(scanName.split("\\.")[1]).getPeakList();//fileId.scanId.scanNum.coId
        } finally {
            lock.unlock();
        }
        TreeMap<Double, Double> plMap = specProcessor.preSpectrumTopNStyle(rawPLMap, precursorMass, precursorCharge, minClear, maxClear, DatasetReader.topN);

        if (plMap.isEmpty()) {
            return null;
        }

        // Coding
        InferSegment inferSegment = buildIndex.getInferSegment();
        TreeMap<Double, Double> finalPlMap = inferSegment.addVirtualPeaks(precursorMass, plMap);
        Map<String, String> protIdSeqMap = buildIndex.proteinPeptideMap;

        List<ThreeExpAA> tag3List = inferSegment.inferSegmentLocationFromSpectrum(precursorMass, finalPlMap, scanNum);
        List<ThreeExpAA> allLongTagList = inferSegment.getLongTag(finalPlMap, precursorMass - massTool.H2O + MassTool.PROTON, scanNum);
//        List<ThreeExpAA> tag4List = inferSegment.getAllTag4(precursorMass, finalPlMap, scanNum);
//        tag4List.sort(Comparator.comparingDouble(ThreeExpAA::getTotalIntensity).reversed());
        if (!tag3List.isEmpty()) {
//            System.out.println(scanNum+",top");
            Entry entry = new Entry();
            List<ThreeExpAA> ncTags = new ArrayList<>();

            if (!allLongTagList.isEmpty()) {
                //which top long tags to uses
                double scoreCutOff = allLongTagList.get(0).getTotalIntensity() - 3;
                int rankCutOff = 10;
                int validTagNum = 0;
                for (int ii = 0; ii < Math.min(rankCutOff, allLongTagList.size()); ii++) {
                    if (allLongTagList.get(ii).getTotalIntensity() > scoreCutOff) {
                        validTagNum++;
                    } else {
                        break;
                    }
                }

//                Map<Double, Integer> mzIdMap = new HashMap<>();
//                int i = 0;
//                for (double mz : finalPlMap.keySet()) {
//                    mzIdMap.put(mz, i);
//                    i++;
//                }
//                TreeMap<Integer, Double> nodeMap = new TreeMap<>();
//                Map<Integer, Map<Integer, Double>> inEdgeMap = new HashMap<>();
//                Map<Integer, Map<Integer, Double>> outEdgeMap = new HashMap<>();
//                for (ThreeExpAA tag : allLongTagList) {
//                    for (ExpAA aa : tag.getExpAAs()) {
//                        if ((aa.getHeadLocation() == MassTool.PROTON + massTool.H2O) && (!"KR".contains(aa.getAA()))) {
//                            continue;
//                        }
//                        int n1 = mzIdMap.get(aa.getHeadLocation());
//                        int n2 = mzIdMap.get(aa.getTailLocation());
//                        double inte1 = aa.getHeadIntensity();
//                        double inte2 = aa.getTailIntensity();
//                        nodeMap.put(n1, inte1);
//                        nodeMap.put(n2, inte2);
//                        if (outEdgeMap.containsKey(n1)) {
//                            outEdgeMap.get(n1).put(n2, 0.5 * (inte1 + inte2));
//                        } else {
//                            Map<Integer, Double> temp = new HashMap<>();
//                            temp.put(n2, 0.5 * (inte1 + inte2));
//                            outEdgeMap.put(n1, temp);
//                        }
//
//                        if (inEdgeMap.containsKey(n2)) {
//                            inEdgeMap.get(n2).put(n1, 0.5 * (inte1 + inte2));
//                        } else {
//                            Map<Integer, Double> temp = new HashMap<>();
//                            temp.put(n1, 0.5 * (inte1 + inte2));
//                            inEdgeMap.put(n2, temp);
//                        }
//
//                    }
//                }
//                for (int node : nodeMap.keySet()) {
//                    if (!inEdgeMap.containsKey(node)) {
//                        for (int n2 : outEdgeMap.get(node).keySet()) {
//                            outEdgeMap.get(node).put(n2, outEdgeMap.get(node).get(n2) + nodeMap.get(node) / 2);
//                            inEdgeMap.get(n2).put(node, inEdgeMap.get(n2).get(node) + nodeMap.get(node) / 2);
//                        }
//                    }
//                    if (!outEdgeMap.containsKey(node)) {
//                        for (int n1 : inEdgeMap.get(node).keySet()) {
//                            outEdgeMap.get(n1).put(node, outEdgeMap.get(n1).get(node) + nodeMap.get(node) / 2);
//                            inEdgeMap.get(node).put(n1, inEdgeMap.get(node).get(n1) + nodeMap.get(node) / 2);
//                        }
//                    }
//                }
//                // finish prepare graph
//
//                Set<Edge> edgeToDel = new HashSet<>();
//                Map<Integer, Integer> tempToDel = new HashMap<>();
//
//                for (int node : nodeMap.descendingKeySet()) {
//                    if (outEdgeMap.containsKey(node)) {
//                        Map<Integer, Double> tempOutMap = outEdgeMap.get(node);
//                        double maxOut = 0d;
//                        int maxOutPeak = -1;
//                        for (int n2 : tempOutMap.keySet()) {
//                            if (tempOutMap.get(n2) > maxOut) {
//                                maxOut = tempOutMap.get(n2);
//                                maxOutPeak = n2;
//                            }
//                        }
//                        for (int n2 : tempOutMap.keySet()) {
//                            if (n2 != maxOutPeak) {
//                                edgeToDel.add(new Edge(node, n2));
//                                tempToDel.put(node, n2);
//
//                            }
//                        }
//                        if (inEdgeMap.containsKey(node)) {
//
//                            for (int n1 : inEdgeMap.get(node).keySet()) {
//                                outEdgeMap.get(n1).put(node, outEdgeMap.get(n1).get(node) + maxOut);
//                                inEdgeMap.get(node).put(n1, inEdgeMap.get(node).get(n1) + maxOut);
//
//                            }
//                        }
//                    }
//                }
//
//                for (int n1 : tempToDel.keySet()) {
//                    int n2 = tempToDel.get(n1);
//                    inEdgeMap.get(n2).remove(n1);
//                    outEdgeMap.get(n1).remove(n2);
//                }
//                for (int node : nodeMap.keySet()) {
//                    if (inEdgeMap.containsKey(node)) {
//                        Map<Integer, Double> tempInMap = inEdgeMap.get(node);
//                        double maxIn = 0d;
//                        int maxInPeak = -1;
//                        for (int n1 : tempInMap.keySet()) {
//                            if (tempInMap.get(n1) > maxIn) {
//                                maxIn = tempInMap.get(n1);
//                                maxInPeak = n1;
//                            }
//                        }
//                        for (int n1 : tempInMap.keySet()) {
//                            if (n1 != maxInPeak) {
//                                edgeToDel.add(new Edge(n1, node));
//
//                            }
//                        }
//                        if (outEdgeMap.containsKey(node)) {
//
//                            for (int n2 : outEdgeMap.get(node).keySet()) {
//                                outEdgeMap.get(node).put(n2, outEdgeMap.get(node).get(n2) + maxIn);
//                                inEdgeMap.get(n2).put(node, inEdgeMap.get(n2).get(node) + maxIn);
//
//                            }
//                        }
//                    }
//                }
//                Set<Pair<Integer, Integer>> pairToDel = new HashSet<>();
//                for (Edge e : edgeToDel) {
//                    pairToDel.add(new Pair(e.n1, e.n2));
//                }
//                //        List<ThreeExpAA> allLongDenoisedTagList = inferSegment.getLongDenoisedTag(finalPlMap, precursorMass - massTool.H2O + MassTool.PROTON, scanNum, pairToDel);
                FMIndex fmIndex = buildIndex.fmIndex;

                ArrayList<Pair<Integer, Double>> posScoreList = new ArrayList<>();
                Set<String> usedLongTags = new HashSet<>();
                Map<Integer, String> posTagMap = new HashMap<>();
                Queue<Pair<String, Double>> tagQueue = new LinkedList<>();
                Map<String, List<Double>> tagIntesIdRangeMap = new HashMap<>();

                for (ThreeExpAA longTag : allLongTagList.subList(0, Math.min(validTagNum, allLongTagList.size()))) {
                    //            for (ThreeExpAA longTag : allLongTagList) {
                    if (longTag.getTotalIntensity() < longTag.size() * 0.8 + 0.8 && longTag.size() < 10) continue;

                    List<Double> intensityList = new ArrayList<>();
                    for (double inte : longTag.intensityArray) {
                        intensityList.add(inte);
                    }
                    String oriTagStr = longTag.getPtmFreeAAString().replace('#', 'L');
                    tagIntesIdRangeMap.put(oriTagStr, intensityList);
                    tagQueue.offer(new Pair<>(oriTagStr, longTag.getTotalIntensity()));
                }
                //tagQueue.
                Set<String> tagsThisRound = new HashSet<>();
                while (!tagQueue.isEmpty()) {
                    Pair<String, Double> tagPair = tagQueue.poll();
                    int ptnForwardCount = 0;
                    int ptnBackwardCount = 0;
                    SearchInterval searchForward = null;
                    SearchInterval searchBackward = null;
                    String tagStr = tagPair.getFirst();
                    if (!usedLongTags.contains(tagStr)) {
                        usedLongTags.add(tagStr);
                        char[] tagChar = tagStr.toCharArray();
                        searchForward = fmIndex.fmSearch(tagChar);
                        if (searchForward != null) {
                            ptnForwardCount = searchForward.ep - searchForward.sp + 1;
                        }
                    }
                    String revTagStr = new StringBuilder(tagStr).reverse().toString();
                    if (!usedLongTags.contains(revTagStr)) {
                        usedLongTags.add(revTagStr);
                        char[] revTagChar = revTagStr.toCharArray();
                        searchBackward = fmIndex.fmSearch(revTagChar);
                        if (searchBackward != null) {
                            ptnBackwardCount = searchBackward.ep - searchBackward.sp + 1;
                        }
                    }

//                    if (ptnForwardCount > 0 && ptnBackwardCount > 0) {
//                        int a = 1;
//                    }

                    if (ptnForwardCount + ptnBackwardCount > 0) {
                        Set<String> thisRoundProts = new HashSet<>();
                        if (ptnForwardCount > 0) {
                            for (int ii = searchForward.sp; ii <= searchForward.ep; ii++) {
                                int res = Arrays.binarySearch(buildIndex.dotPosArr, fmIndex.SA[ii]);
                                thisRoundProts.add(buildIndex.posProtMap.get(-res - 2));
                            }
                        }
                        if (ptnBackwardCount > 0) {
                            for (int ii = searchBackward.sp; ii <= searchBackward.ep; ii++) {
                                int res = Arrays.binarySearch(buildIndex.dotPosArr, fmIndex.SA[ii]);
                                thisRoundProts.add(buildIndex.posProtMap.get(-res - 2));
                            }
                        }
                        Set<String> group = new HashSet<>();
                        for (String prot1 : thisRoundProts) {
                            String seq1 = protIdSeqMap.get(prot1);
                            Set<String> tempGroup = new HashSet<>();
                            tempGroup.add(prot1);
                            for (String prot2 : thisRoundProts) {
                                if (prot1.contentEquals(prot2)) continue;
                                String seq2 = protIdSeqMap.get(prot2);
                                if (seq1.length() != seq2.length()) continue;
                                int errorNum = 0;
                                for (int iChar = 0; iChar < seq1.length(); iChar++) {
                                    if (seq1.charAt(iChar) != seq2.charAt(iChar)) {
                                        errorNum++;
                                    }
                                }
                                if (errorNum < 0.1 * seq1.length()) {
                                    tempGroup.add(prot2);
                                }
                            }
                            if (tempGroup.size() > 1) {
                                group.addAll(tempGroup);
                                break;
                            }
                        }

                        if (ptnForwardCount > 0) {
                            for (int ii = searchForward.sp; ii <= searchForward.ep; ii++) {
                                posScoreList.add(new Pair<>(fmIndex.SA[ii], tagPair.getSecond() * 1 / (ptnForwardCount + ptnBackwardCount - group.size() + 1)));
                                posTagMap.put(fmIndex.SA[ii], tagStr);
                            }
                        }
                        if (ptnBackwardCount > 0) {
                            for (int ii = searchBackward.sp; ii <= searchBackward.ep; ii++) {
                                posScoreList.add(new Pair<>(fmIndex.SA[ii], tagPair.getSecond() * 1 / (ptnForwardCount + ptnBackwardCount - group.size() + 1)));
                                posTagMap.put(fmIndex.SA[ii], revTagStr);
                            }
                        }
                        break; // if it is processing original tags, dont break
                    }
                    if (tagStr.length() > 6) tagsThisRound.add(tagStr);

                    if (tagQueue.isEmpty()) {
                        List<Pair<String, Double>> tempList = new ArrayList<>();
                        for (String tag : tagsThisRound) {
                            List<Double> thisIntesList = tagIntesIdRangeMap.get(tag);
                            if (!tagIntesIdRangeMap.containsKey(tag)) {
                                int a = 1;
                            }
                            String subTagL = tag.substring(0, tag.length() - 1);
                            String subTagR = tag.substring(1);
                            double scoreL = 0;
                            double scoreR = 0;
                            for (double s : thisIntesList.subList(0, tag.length())) scoreL += s;
                            for (double s : thisIntesList.subList(1, tag.length() + 1)) scoreR += s;

                            if (!usedLongTags.contains(subTagL)) {
                                tempList.add(new Pair<>(subTagL, scoreL));
                                tagIntesIdRangeMap.put(subTagL, thisIntesList.subList(0, tag.length()));
                            }
                            if (!usedLongTags.contains(subTagR)) {
                                tempList.add(new Pair<>(subTagR, scoreR));
                                tagIntesIdRangeMap.put(subTagR, thisIntesList.subList(1, tag.length() + 1));
                            }
                        }
                        tempList.sort(Comparator.comparingDouble(Pair::getSecond));
                        for (int j = tempList.size() - 1; j >= 0; j--) {
                            tagQueue.offer(tempList.get(j));
                        }
                    }
                }

                double highestScore = 0;
                Map<String, Double> protScoreMap = new HashMap<>();
                Map<String, List<Pair<String, Double>>> protTagScoreMap = new HashMap<>();
                Map<String, Double> tagScoreMap = new HashMap<>();
                Map<String, Integer> protTimesMap = new HashMap<>();
                for (Pair<Integer, Double> posScore : posScoreList) {
                    int res = Arrays.binarySearch(buildIndex.dotPosArr, posScore.getFirst());
                    if (res > 0) {
                        System.out.println("located on dot sign," + scanNum);
                    }
                    String prot = buildIndex.posProtMap.get(-res - 2);
                    if (protScoreMap.containsKey(prot)) {
                        protScoreMap.put(prot, protScoreMap.get(prot) + posScore.getSecond());
                        protTimesMap.put(prot, protTimesMap.get(prot) + 1);
                        protTagScoreMap.get(prot).add(new Pair<>(posTagMap.get(posScore.getFirst()), posScore.getSecond()));
                        tagScoreMap.put(posTagMap.get(posScore.getFirst()), posScore.getSecond());
                        if (protScoreMap.get(prot) + posScore.getSecond() > highestScore)
                            highestScore = protScoreMap.get(prot) + posScore.getSecond();
                    } else {
                        protScoreMap.put(prot, posScore.getSecond());
                        protTimesMap.put(prot, 1);
                        List<Pair<String, Double>> temp = new ArrayList<>();
                        temp.add(new Pair<>(posTagMap.get(posScore.getFirst()), posScore.getSecond()));
                        protTagScoreMap.put(prot, temp);
                        tagScoreMap.put(posTagMap.get(posScore.getFirst()), posScore.getSecond());
                        if (posScore.getSecond() > highestScore) highestScore = posScore.getSecond();
                    }
                }


                //            entry.candidateList = tagSeqList;
//                Map<String, Double> top1Pep = new HashMap<>();
//                for (String prot : protScoreMap.keySet()) {
//                    if (protScoreMap.get(prot) > highestScore * 0.9) {
//                        top1Pep.put(prot, protScoreMap.get(prot));
//                    }
//                }
//                entry.protScoreMap = top1Pep;
                entry.protTagScoreMap = protTagScoreMap;
            }

            entry.scanCode = inferSegment.generateSegmentIntensityVector(tag3List);
            entry.scanName = this.scanName;
            return entry;
        } else {
//            System.out.println("nullLsz, "+ scanNum);
            return null;
        }
    }

//    public class Edge {
//        public int n1;
//        public int n2;
//        public double weight = 0d;
//        Edge(int n1, int n2){
//            this.n1 = n1;
//            this.n2 = n2;
//        };
//
//        @Override
//        public int hashCode(){
//            return n1*n2;
//        }
//        @Override
//        public boolean equals(Object other) {
//            if (this.n1 == 31 && this.n2 == 36) {
//                int a = 1;
//            }
//            Edge temp = (Edge) other;
//            return (this.n1 == temp.n1) && (this.n2 == temp.n2);
//        }
//    }
    public class Entry {

        public int scanNum = SpecCoder.this.scanNum;
        public double precursorMass = SpecCoder.this.precursorMass;
        public Map<String, Double> protScoreMap = new HashMap<>();

        public Map<String, List<Pair<String, Double>>> protTagScoreMap = new HashMap<>();
        public String scanName;
        public SparseVector scanCode;
        Entry() {
        }
    }
}