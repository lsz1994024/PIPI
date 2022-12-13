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
import org.apache.commons.math3.util.Pair;
import proteomics.FM.FMIndex;
import proteomics.FM.SearchInterval;
import proteomics.Index.BuildIndex;
import proteomics.Segment.InferSegment;
import proteomics.Spectrum.DatasetReader;
import proteomics.Types.ExpTag;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;

import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.locks.ReentrantLock;

import static proteomics.PIPI.lszDebugScanNum;
import static proteomics.PIPI.minTagLenToReduceProtDb;

public class GetLongTag implements Callable<GetLongTag.Entry> {
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


    public GetLongTag(int scanNum, BuildIndex buildIndex, MassTool massTool, JMzReader spectraParser, double minClear, double maxClear, ReentrantLock lock
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
            rawPLMap = spectraParser.getSpectrumById(scanName.split("\\.")[1]).getPeakList();//fileId.scanId.scanNum
        } finally {
            lock.unlock();
        }
        TreeMap<Double, Double> plMap = specProcessor.preSpectrumTopNStyleWithChargeLimit(rawPLMap, precursorMass, precursorCharge, minClear, maxClear, DatasetReader.topN);

        if (plMap.isEmpty()) {
            return null;
        }

        // Coding
        InferSegment inferSegment = buildIndex.getInferSegment();
        TreeMap<Double, Double> finalPlMap = inferSegment.addVirtualPeaks(precursorMass, plMap);
        Map<String, String> protIdSeqMap = buildIndex.protSeqMap;
        if (lszDebugScanNum.contains(this.scanNum)) {
            int a = 1;
        }
//        System.out.println("in, "+ scanNum);
        List<ExpTag> allLongTagList = inferSegment.getLongTag(finalPlMap, precursorMass - massTool.H2O + MassTool.PROTON, scanNum, minTagLenToReduceProtDb);
//        System.out.println("out, "+ scanNum);

        if (!allLongTagList.isEmpty()) {
            Entry entry = new Entry();
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

            FMIndex fmIndex = buildIndex.fmIndexFull;

            ArrayList<Pair<Integer, Double>> posScoreList = new ArrayList<>();
            Set<String> usedLongTags = new HashSet<>();
            Map<Integer, String> posTagMap = new HashMap<>();
            Queue<Pair<String, Double>> tagQueue = new LinkedList<>();
            Map<String, List<Double>> tagIntesIdRangeMap = new HashMap<>();

            for (ExpTag longTag : allLongTagList.subList(0, Math.min(validTagNum, allLongTagList.size()))) {
                //            for (ThreeExpAA longTag : allLongTagList) {
                if (longTag.getTotalIntensity() < longTag.size() * 0.8 + 0.8 && longTag.size() < 10) continue;

                List<Double> intensityList = new ArrayList<>();
                for (double inte : longTag.intensityArray) {
                    intensityList.add(inte);
                }
                String oriTagStr = longTag.getFreeAaString().replace('#', 'L');
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

                if (ptnForwardCount + ptnBackwardCount > 0) {
                    Set<String> thisRoundProts = new HashSet<>();
                    if (ptnForwardCount > 0) {
                        for (int ii = searchForward.sp; ii <= searchForward.ep; ii++) {
                            int res = Arrays.binarySearch(buildIndex.dotPosArrFull, fmIndex.SA[ii]);
                            thisRoundProts.add(buildIndex.posProtMapFull.get(-res - 2));
                        }
                    }
                    if (ptnBackwardCount > 0) {
                        for (int ii = searchBackward.sp; ii <= searchBackward.ep; ii++) {
                            int res = Arrays.binarySearch(buildIndex.dotPosArrFull, fmIndex.SA[ii]);
                            thisRoundProts.add(buildIndex.posProtMapFull.get(-res - 2));
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
                int res = Arrays.binarySearch(buildIndex.dotPosArrFull, posScore.getFirst());
                if (res > 0) {
                    System.out.println("located on dot sign," + scanNum);
                }
                String prot = buildIndex.posProtMapFull.get(-res - 2);
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

            entry.protTagScoreMap = protTagScoreMap;

            entry.scanName = this.scanName;
            return entry;
        } else {
            return null;
        }
    }

    public class Entry {
        public Map<String, List<Pair<String, Double>>> protTagScoreMap = new HashMap<>();
        public String scanName;
        Entry() {
        }
    }
}
