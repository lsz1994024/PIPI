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
import ProteomicsLibrary.Types.SparseBooleanVector;
import ProteomicsLibrary.Types.SparseVector;
import org.apache.commons.math3.util.Pair;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.FM.FMIndex;
import proteomics.FM.SearchInterval;
import proteomics.Index.BuildIndex;
import proteomics.PTM.InferPTM;
import proteomics.Segment.InferSegment;
import proteomics.Spectrum.DatasetReader;
import proteomics.Types.*;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;

import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.locks.ReentrantLock;

import static proteomics.PIPI.*;
import static proteomics.PTM.InferPTM.*;
import static proteomics.PTM.InferPTM.N_PART;
import static proteomics.Segment.InferSegment.*;

public class PreSearch implements Callable<PreSearch.Entry> {
    private static final Logger logger = LoggerFactory.getLogger(PreSearch.class);

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
    private final double ms2Tolerance;
    private final InferPTM inferPTM;
    private final int minPepLen;
    private final int maxPepLen;
    public PreSearch(int scanNum, BuildIndex buildIndex, MassTool massTool, double ms2Tolerance, double ms1Tolerance, double leftInverseMs1Tolerance, double rightInverseMs1Tolerance
            , int ms1ToleranceUnit, double minPtmMass, double maxPtmMass, int localMaxMs2Charge
            , JMzReader spectraParser, double minClear, double maxClear, ReentrantLock lock, String scanName, int precursorCharge, double precursorMass
            , SpecProcessor specProcessor, String truth, int minPepLen, int maxPepLen) {

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
        this.ms2Tolerance = ms2Tolerance;
        this.inferPTM = buildIndex.getInferPTM();
        this.minPepLen = minPepLen;
        this.maxPepLen = maxPepLen;
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

        double ms1TolAbs = Double.parseDouble(InferPTM.df3.format(precursorMass*ms1Tolerance/1000000));
        if (lszDebugScanNum.contains(this.scanNum)) {
            System.out.println(scanNum + ", entered");
            int a = 1;
        }
        TreeMap<Double, Double> plMap = specProcessor.preSpectrumTopNStyleWithChargeLimit(rawPLMap, precursorMass, precursorCharge, minClear, maxClear, DatasetReader.topN, ms2Tolerance);

        if (plMap.isEmpty()) return null;
        // Coding
        SparseVector expProcessedPL;
        if (PIPI.useXcorr) {
            expProcessedPL = specProcessor.prepareXCorr(plMap, false);
        } else {
            expProcessedPL = specProcessor.digitizePL(plMap);
        }
        InferSegment inferSegment = buildIndex.getInferSegment();
        TreeMap<Double, Double> finalPlMap = inferSegment.addVirtualPeaks(precursorMass, plMap);
        Map<String, Pair<Integer,Integer>> varTagLenParaMap = new HashMap<>();
        varTagLenParaMap.put("3", new Pair<>(3,3));
        varTagLenParaMap.put("4", new Pair<>(4,4));
        varTagLenParaMap.put("5", new Pair<>(5,5));
        varTagLenParaMap.put("6", new Pair<>(6,6));
        varTagLenParaMap.put("7", new Pair<>(7,7));
        varTagLenParaMap.put("8", new Pair<>(8,8));
        varTagLenParaMap.put("9", new Pair<>(9,9));
        varTagLenParaMap.put("V", new Pair<>(3,9));
        double totalMass = precursorMass + 2 * MassTool.PROTON;
        FMIndex fmIndex = buildIndex.fmIndexFull;

        if (lszDebugScanNum.contains(this.scanNum)) {
            System.out.println(scanNum + ", entered");
            int a = 1;
        }
        Entry entry = new Entry(scanNum,scanName);
        for (String lenStr : varTagLenParaMap.keySet()){
            int myminTagLenToExtract = varTagLenParaMap.get(lenStr).getFirst();
            int mymaxTagLenToExtract = varTagLenParaMap.get(lenStr).getSecond();

            List<ExpTag> allLongTagList = inferSegment.getLongTag(finalPlMap, precursorMass - massTool.H2O + MassTool.PROTON, scanNum, myminTagLenToExtract,mymaxTagLenToExtract);
//            List<ExpTag> cleanedAllLongTagList = inferSegment.cleanAbundantTagsPrefix(allLongTagList, myminTagLenToExtract);
//            allLongTagList = cleanedAllLongTagList;
//            for (int i = 0; i < allLongTagList.size();i++){
//                String tagSeq = allLongTagList.get(i).getFreeAaString();
//
//            }
            Set<String> searchedTagStrSet = new HashSet<>();
            int minTagLen = 4;
            Map<String, Double> pepScoreMap = new HashMap<>();
            for (ExpTag tagInfo : allLongTagList.subList(0, Math.min(200, allLongTagList.size()))){
                minTagLen = tagInfo.size() > 4 ? 5 : 4;
                if (buildIndex.posProtMapFull.size() < 5000) { // todo this is for synthetic only
                    minTagLen = 3;
                }
                String tagStr = tagInfo.getFreeAaString();
                String revTagStr = new StringBuilder(tagStr).reverse().toString();

                if (tagInfo.isNorC == N_TAG) { //n tag
                    String tagStrMzStr = tagStr + df3.format(tagInfo.getHeadLocation());
                    if (!searchedTagStrSet.contains(tagStrMzStr)) {
                        char[] tagChar = tagStr.toCharArray();

                        Set<String> protIdSetByThisTag = new HashSet<>();
                        int n_res = searchAndSaveFuzzy(scanNum, tagInfo, lenStr, pepScoreMap, protIdSetByThisTag, fmIndex, tagChar, minTagLen, expProcessedPL,finalPlMap, true);
                        searchedTagStrSet.add(tagStrMzStr);
                        if (tagStr.length() > minTagLen+1 && lenStr.contentEquals("V")){ // if the tag was already long i.e. there is space to sub
                            for (int n_cTermAaToCut = 1; n_cTermAaToCut <= Math.min((tagChar.length)/2,tagChar.length-minTagLen); n_cTermAaToCut++){
                                String subTagStr = tagStr.substring(0, tagStr.length()-n_cTermAaToCut);
                                char[] subTagChar = subTagStr.toCharArray();
                                ExpTag subTagInfo = tagInfo.subTag(0,tagStr.length()-n_cTermAaToCut);
                                if (!searchedTagStrSet.contains(tagStr)) {
                                    int numResSub = searchAndSaveFuzzy(scanNum ,subTagInfo, lenStr, pepScoreMap, protIdSetByThisTag, fmIndex, subTagChar, minTagLen, expProcessedPL, finalPlMap, false);
                                    searchedTagStrSet.add(tagStr);
                                    if (numResSub > 1) {
                                        break;
                                    }
                                }
                            }
                        }
                    }
                } else if (tagInfo.isNorC == C_TAG) { // c tag
                    char[] revTagChar = revTagStr.toCharArray();
                    ExpTag revTagInfo = tagInfo.revTag(totalMass);
                    String revTagStrMzStr = revTagStr + df3.format(revTagInfo.getHeadLocation());
                    if (!searchedTagStrSet.contains(revTagStrMzStr)) {
                        Set<String> protIdSetByThisTag = new HashSet<>();
                        int n_res = searchAndSaveFuzzy(scanNum ,revTagInfo, lenStr, pepScoreMap, protIdSetByThisTag, fmIndex, revTagChar, minTagLen, expProcessedPL, finalPlMap, true);
                        searchedTagStrSet.add(revTagStrMzStr);
                        if (revTagStr.length() > minTagLen+1 && lenStr.contentEquals("V")){ // if the tag was already long i.e. there is space to sub
                            for (int n_cTermAaToCut = 1; n_cTermAaToCut <= Math.min((revTagChar.length)/2,revTagChar.length-minTagLen); n_cTermAaToCut++){
                                String subRevTagStr = revTagStr.substring(0, tagStr.length()-n_cTermAaToCut);
                                char[] subRevTagChar = subRevTagStr.toCharArray();
                                ExpTag subRevTagInfo = revTagInfo.subTag(0,tagStr.length()-n_cTermAaToCut);
                                if (!searchedTagStrSet.contains(subRevTagStr)) {
                                    subRevTagInfo.isNorC = NON_NC_TAG; // if c Tag is cut, it should wont be cTag
                                    int numResSub = searchAndSaveFuzzy(scanNum ,subRevTagInfo, lenStr, pepScoreMap, protIdSetByThisTag, fmIndex, subRevTagChar, minTagLen, expProcessedPL, finalPlMap, false);
                                    searchedTagStrSet.add(subRevTagStr);
                                    if (numResSub > 1) {
                                        break;
                                    }
                                }
                            }
                        }
                    }
                } else { // non-nc tag
                    char[] tagChar = tagStr.toCharArray();
                    String tagStrMzStr = tagStr + df3.format(tagInfo.getHeadLocation());
                    if (!searchedTagStrSet.contains(tagStrMzStr)) {
                        Set<String> protIdSetByThisTag = new HashSet<>();
                        int n_res = searchAndSaveFuzzy(scanNum ,tagInfo, lenStr, pepScoreMap, protIdSetByThisTag, fmIndex, tagChar, minTagLen, expProcessedPL,finalPlMap, true);
                        searchedTagStrSet.add(tagStrMzStr);
//                        if (n_res < 100 ) {
                            if (tagStr.length() > minTagLen+1 && lenStr.contentEquals("V")){ // if the tag was already long i.e. there is space to sub
                                for (int n_cTermAaToCut = 1; n_cTermAaToCut <= Math.min((tagChar.length)/2,tagChar.length-minTagLen); n_cTermAaToCut++){
                                    String subTagStr = tagStr.substring(0, tagStr.length()-n_cTermAaToCut);
                                    ExpTag subTagInfo = tagInfo.subTag(0,tagStr.length()-n_cTermAaToCut);
                                    //sub forward
                                    char[] subTagChar = subTagStr.toCharArray();
                                    int numResSub1 = 0;
                                    if (!searchedTagStrSet.contains(subTagStr)) {
                                        numResSub1 = searchAndSaveFuzzy(scanNum ,subTagInfo, lenStr, pepScoreMap, protIdSetByThisTag, fmIndex, subTagChar, minTagLen, expProcessedPL,finalPlMap, false);
                                        searchedTagStrSet.add(subTagStr);
                                        if (numResSub1 > 1) {
                                            break;
                                        }
                                    }
                                }
                            }
//                        }
                    }
                    char[] revTagChar = revTagStr.toCharArray();
                    ExpTag revTagInfo = tagInfo.revTag(totalMass);
                    String revTagStrMzStr = revTagStr + df3.format(revTagInfo.getHeadLocation());
                    if (!searchedTagStrSet.contains(revTagStrMzStr)) {
                        Set<String> protIdSetByThisTag = new HashSet<>();
                        int n_res = searchAndSaveFuzzy(scanNum ,revTagInfo, lenStr, pepScoreMap, protIdSetByThisTag, fmIndex, revTagChar, minTagLen, expProcessedPL, finalPlMap, true);
                        searchedTagStrSet.add(revTagStrMzStr);
//                        if (n_res < 100 ) {
                            if (revTagStr.length() > minTagLen+1 && lenStr.contentEquals("V")){ // if the tag was already long i.e. there is space to sub
                                for (int n_cTermAaToCut = 1; n_cTermAaToCut <= Math.min((revTagChar.length)/2,revTagChar.length-minTagLen); n_cTermAaToCut++){
                                    String subTagStr = revTagStr.substring(0, revTagStr.length()-n_cTermAaToCut);
                                    ExpTag subTagInfo = tagInfo.revTag(totalMass).subTag(0,revTagStr.length()-n_cTermAaToCut);

                                    //sub forward
                                    char[] subTagChar = subTagStr.toCharArray();
                                    int numResSub1 = 0;
                                    if (!searchedTagStrSet.contains(subTagStr)) {
                                        numResSub1 = searchAndSaveFuzzy(scanNum ,subTagInfo, lenStr, pepScoreMap, protIdSetByThisTag, fmIndex, subTagChar, minTagLen, expProcessedPL, finalPlMap, false);
                                        searchedTagStrSet.add(subTagStr);
                                        if (numResSub1 > 1) {
                                            break;
                                        }
                                    }
                                }
                            }
//                        }
                    }
                }
            } // finish collecting pepScore Map for current lenstr
            List<String> pepStrList = new ArrayList<>(pepScoreMap.keySet());
            Collections.sort(pepStrList, Comparator.comparing(o->pepScoreMap.get(o), Comparator.reverseOrder()));
            double truthGoodness = 0;  // how many candidates has higher score than the truth

            if (allLongTagList.isEmpty()) {
                truthGoodness = 9999; // error   no tags
            } else {
                if (pepScoreMap.containsKey(truth)) {
                    truthGoodness = Math.min(8000, pepStrList.indexOf(truth));
                } else {
                    truthGoodness = 9999; // error  has tags but truth not find
                }
            }

            entry.varLenTruthGoodnessMap.put(lenStr, truthGoodness);
        }

        return entry;
    }


    private int searchAndSaveFuzzy(int scanNum, ExpTag tagInfo, String varLenStr, Map<String, Double> pepScoreMap
            , Set<String> protIdSetByThisTag, FMIndex fmIndex, char[] tagChar, int minTagLen, SparseVector expProcessedPL,TreeMap<Double, Double> plMap, boolean isFirstTime) throws CloneNotSupportedException {
        SearchInterval searchRes = fmIndex.fmSearchFuzzy(tagChar);
        if (searchRes == null) {
            return 0;
        }
        int numRes = searchRes.ep-searchRes.sp+1;
        int solCount = 0;
        if (searchRes.settled) {
            for (int ii = searchRes.sp; ii <= searchRes.ep; ii++) {
                int absTagPos = fmIndex.SA[ii];
                int dotIndex = -2-Arrays.binarySearch(buildIndex.dotPosArrFull, absTagPos);

                String protId = buildIndex.posProtMapFull.get(dotIndex); //protId is the same as pepseq
                String protSeq = buildIndex.protSeqMap.get(protId);
//                if (lszDebugScanNum.contains(this.scanNum) && protSeq.contentEquals(truth)) {
//                    int a = 1;
//                }
                if (protIdSetByThisTag.contains(protId)){
                    continue;
                }
                solCount++;
                protIdSetByThisTag.add(protId);

                double tagScore = tagInfo.getTotalIntensity() / numRes;
                if (pepScoreMap.containsKey(protSeq)) {
                    pepScoreMap.put(protSeq, pepScoreMap.get(protSeq) + tagScore);
                } else {
                    pepScoreMap.put(protSeq, tagScore);
                }
            }
        } else {
            if (!isFirstTime || !varLenStr.contentEquals("V")) return 0;
            int matchedPos = searchRes.matchedPos;
            if (tagInfo.size()-matchedPos < minTagLen) return 0;
            for (int ii = searchRes.sp; ii <= searchRes.ep; ii++) {
                int absTagPos = fmIndex.SA[ii];
                int dotIndex = -2-Arrays.binarySearch(buildIndex.dotPosArrFull, absTagPos);
                String protId = buildIndex.posProtMapFull.get(dotIndex);
                String protSeq = buildIndex.protSeqMap.get(protId);
                if (lszDebugScanNum.contains(this.scanNum) && protSeq.contentEquals(truth)) {
                    int a = 1;
                }
                double tagScore = tagInfo.subTag(matchedPos, tagChar.length).getTotalIntensity() / numRes;

                if (pepScoreMap.containsKey(protSeq)) {
                    pepScoreMap.put(protSeq, pepScoreMap.get(protSeq) + tagScore);
                } else {
                    pepScoreMap.put(protSeq, tagScore);
                }
                solCount++;
            }

        }
        return solCount;

    }

    public class Entry { //copied from PtmSearch
        final int scanNum;
        final String scanName;
        final Map<String, Double> varLenTruthGoodnessMap = new HashMap<>();
        Entry(int scanNum, String scanName) {
            this.scanNum = scanNum;
            this.scanName = scanName;
        }
    }
}
