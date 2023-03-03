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
    public PreSearch(int scanNum, BuildIndex buildIndex, MassTool massTool, double ms2Tolerance, double ms1Tolerance, double leftInverseMs1Tolerance, double rightInverseMs1Tolerance
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
        this.ms2Tolerance = ms2Tolerance;
        this.inferPTM = buildIndex.getInferPTM();
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


        double minPcMass = -1 * ms1Tolerance;
        double maxPcMass = ms1Tolerance;
        if (ms1ToleranceUnit == 1) {
            minPcMass = (precursorMass * leftInverseMs1Tolerance);
            maxPcMass = (precursorMass * rightInverseMs1Tolerance);
        }



        double ms1TolAbs = Double.parseDouble(InferPTM.df.format(precursorMass*ms1Tolerance/1000000));

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

//        List<ExpTag> allLongTagList = inferSegment.getLongTag(finalPlMap, precursorMass - massTool.H2O + MassTool.PROTON, scanNum, minTagLenToExtract,maxTagLenToExtract);
        List<ExpTag> allLongTagList = inferSegment.getLongTag(finalPlMap, precursorMass - massTool.H2O + MassTool.PROTON, scanNum, minTagLenToExtract,maxTagLenToExtract);
        if (lszDebugScanNum.contains(this.scanNum)) {
            System.out.println(scanNum + ", entered");
            int a = 1;
        }

        List<ExpTag> cleanedAllLongTagList = inferSegment.cleanAbundantTagsPrefix(allLongTagList, minTagLenToExtract);

        allLongTagList = cleanedAllLongTagList;
        if (allLongTagList.isEmpty())  return null;

        double totalMass = precursorMass + 2 * MassTool.PROTON;

        FMIndex fmIndex = buildIndex.fmIndexReduced;

        Map<String, List<Peptide>> resPeptideListMap = new HashMap<>();
        Map<String, PeptideInfo> peptideInfoMap = new HashMap<>(50000);

        Set<String> searchedTagStrSet = new HashSet<>();
        int minTagLen = 4;

        int n_tags = 0;
        for (ExpTag tagInfo : allLongTagList.subList(0, Math.min(10, allLongTagList.size()))){

            minTagLen = tagInfo.size() > 4 ? 5 : 4;
            if (buildIndex.posProtMapReduced.size() < 5000) { // todo this is for synthetic only
                minTagLen = 3;
            }
            String tagStr = tagInfo.getFreeAaString();
            String revTagStr = new StringBuilder(tagStr).reverse().toString();

            if (tagInfo.isNorC == N_TAG) { //n tag
                if (!searchedTagStrSet.contains(tagStr)) {
                    char[] tagChar = tagStr.toCharArray();

                    Set<String> protIdSetByThisTag = new HashSet<>();
                    int n_res = searchAndSaveFuzzy(scanNum, tagInfo, ms1TolAbs, resPeptideListMap, peptideInfoMap, protIdSetByThisTag, fmIndex, tagChar, minTagLen, expProcessedPL,finalPlMap, true);
                    searchedTagStrSet.add(tagStr);
                    if (n_res > 100) {
                        continue;
                    }
//                    if (n_tags > 5) {
//                        continue;
//                    }
                    if (tagStr.length() > minTagLen+1){ // if the tag was already long i.e. there is space to sub
                        for (int n_cTermAaToCut = 1; n_cTermAaToCut <= Math.min(2,tagChar.length-minTagLen); n_cTermAaToCut++){
                            char[] subTagChar = tagStr.substring(0, tagStr.length()-n_cTermAaToCut).toCharArray();
                            ExpTag subTagInfo = tagInfo.subTag(0,tagStr.length()-n_cTermAaToCut);
                            if (!searchedTagStrSet.contains(tagStr)) {
                                int numResSub = searchAndSaveFuzzy(scanNum ,subTagInfo, ms1TolAbs, resPeptideListMap, peptideInfoMap, protIdSetByThisTag, fmIndex, subTagChar, minTagLen, expProcessedPL, finalPlMap, false);
                                searchedTagStrSet.add(tagStr);
                                if (numResSub > 100) {
                                    break;
                                }
                            }
                        }
                    }
                }
            } else if (tagInfo.isNorC == C_TAG) { // c tag
                char[] revTagChar = revTagStr.toCharArray();
                ExpTag revTagInfo = tagInfo.revTag(totalMass);
                if (!searchedTagStrSet.contains(revTagStr)) {
                    Set<String> protIdSetByThisTag = new HashSet<>();
                    int n_res = searchAndSaveFuzzy(scanNum ,revTagInfo, ms1TolAbs, resPeptideListMap, peptideInfoMap, protIdSetByThisTag, fmIndex, revTagChar, minTagLen, expProcessedPL, finalPlMap, true);
                    searchedTagStrSet.add(revTagStr);
                    if (n_res > 100) {
                        continue;
                    }
//                    if (n_tags > 5) {
//                        continue;
//                    }
                    if (revTagStr.length() > minTagLen+1){ // if the tag was already long i.e. there is space to sub
                        for (int n_cTermAaToCut = 1; n_cTermAaToCut <= Math.min(2,revTagChar.length-minTagLen); n_cTermAaToCut++){
                            String subRevTagStr = revTagStr.substring(0, tagStr.length()-n_cTermAaToCut);
                            char[] subRevTagChar = subRevTagStr.toCharArray();
                            ExpTag subRevTagInfo = revTagInfo.subTag(0,tagStr.length()-n_cTermAaToCut);
                            if (!searchedTagStrSet.contains(subRevTagStr)) {
                                subRevTagInfo.isNorC = NON_NC_TAG; // if c Tag is cut, it should wont be cTag
                                int numResSub = searchAndSaveFuzzy(scanNum ,subRevTagInfo, ms1TolAbs, resPeptideListMap, peptideInfoMap, protIdSetByThisTag, fmIndex, subRevTagChar, minTagLen, expProcessedPL, finalPlMap, false);
                                searchedTagStrSet.add(subRevTagStr);
                                if (numResSub > 100) {
                                    break;
                                }
                            }
                        }
                    }
                }
            } else { // non-nc tag
                char[] tagChar = tagStr.toCharArray();
                if (!searchedTagStrSet.contains(tagStr)) {
                    Set<String> protIdSetByThisTag = new HashSet<>();
                    int n_res = searchAndSaveFuzzy(scanNum ,tagInfo, ms1TolAbs, resPeptideListMap, peptideInfoMap, protIdSetByThisTag, fmIndex, tagChar, minTagLen, expProcessedPL,finalPlMap, true);
                    searchedTagStrSet.add(tagStr);
                    if (n_res < 100 ) {
                        if (tagStr.length() > minTagLen+1){ // if the tag was already long i.e. there is space to sub
                            for (int n_cTermAaToCut = 1; n_cTermAaToCut <= Math.min(2,tagChar.length-minTagLen); n_cTermAaToCut++){
                                String subTagStr = tagStr.substring(0, tagStr.length()-n_cTermAaToCut);
                                ExpTag subTagInfo = tagInfo.subTag(0,tagStr.length()-n_cTermAaToCut);
                                //sub forward
                                char[] subTagChar = subTagStr.toCharArray();
                                int numResSub1 = 0;
                                if (!searchedTagStrSet.contains(subTagStr)) {
                                    numResSub1 = searchAndSaveFuzzy(scanNum ,subTagInfo, ms1TolAbs, resPeptideListMap, peptideInfoMap, protIdSetByThisTag, fmIndex, subTagChar, minTagLen, expProcessedPL,finalPlMap, false);
                                    searchedTagStrSet.add(subTagStr);
                                    if (numResSub1 > 100) {
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
                char[] revTagChar = revTagStr.toCharArray();
                ExpTag revTagInfo = tagInfo.revTag(totalMass);
                if (!searchedTagStrSet.contains(revTagStr)) {
                    Set<String> protIdSetByThisTag = new HashSet<>();
                    int n_res = searchAndSaveFuzzy(scanNum ,revTagInfo, ms1TolAbs, resPeptideListMap, peptideInfoMap, protIdSetByThisTag, fmIndex, revTagChar, minTagLen, expProcessedPL, finalPlMap, true);
                    searchedTagStrSet.add(revTagStr);
                    if (n_res < 100 ) {
                        if (revTagStr.length() > minTagLen+1){ // if the tag was already long i.e. there is space to sub
                            for (int n_cTermAaToCut = 1; n_cTermAaToCut <= Math.min(2,revTagChar.length-minTagLen); n_cTermAaToCut++){
                                String subTagStr = revTagStr.substring(0, revTagStr.length()-n_cTermAaToCut);
                                ExpTag subTagInfo = tagInfo.revTag(totalMass).subTag(0,revTagStr.length()-n_cTermAaToCut);

                                //sub forward
                                char[] subTagChar = subTagStr.toCharArray();
                                int numResSub1 = 0;
                                if (!searchedTagStrSet.contains(subTagStr)) {
                                    numResSub1 = searchAndSaveFuzzy(scanNum ,subTagInfo, ms1TolAbs, resPeptideListMap, peptideInfoMap, protIdSetByThisTag, fmIndex, subTagChar, minTagLen, expProcessedPL, finalPlMap, false);
                                    searchedTagStrSet.add(subTagStr);
                                    if (numResSub1 > 100) {
                                        break;
                                    }
                                }
                            }
                        }
                    }
                }
            }
            n_tags++;
        }

        TreeSet<Peptide> peptideSet = new TreeSet<>(Collections.reverseOrder());
        if (!resPeptideListMap.isEmpty()) {
            for (String pepSeq : resPeptideListMap.keySet()) {
                Peptide bestPep = resPeptideListMap.get(pepSeq).get(0);
                if (bestPep.getScore() > 0) {
                    if (peptideSet.size() < candisNum) {
                        peptideSet.add(bestPep);
                    } else if (bestPep.getScore() > peptideSet.last().getScore()) {
                        peptideSet.pollLast();
                        peptideSet.add(bestPep);
                    }
                }
            }
        }
        if (peptideSet.isEmpty()) {
            return null;
        }
        List<Peptide> pepList = new ArrayList<>(peptideSet);
        Set<Integer> pepIdsToRemove = new HashSet<>();
        for (int j = pepList.size()-1; j >= 1; j--) {
            if (pepIdsToRemove.contains(j)) continue;
            for (int i = 0; i < j; i++) {
                if (pepIdsToRemove.contains(i)) continue;
                if (isHomo(pepList.get(i), pepList.get(j), peptideInfoMap) || pepList.get(j).getScore() > 0.6*pepList.get(i).getScore()) {
                    int iPriority = pepList.get(i).getPriority();
                    int jPriority = pepList.get(j).getPriority();
                    if (iPriority < jPriority) {
                        pepIdsToRemove.add(i);
                    } else  {
                        pepIdsToRemove.add(j);
                    }
                }
            }
            if (pepList.get(j).getPriority() < 0 && isPtmSimuTest ) { // only simu test todo
                pepIdsToRemove.add(j);
            }
        }
        if (pepList.get(0).getPriority() < 0 && isPtmSimuTest ) {// only simu test  todo
            pepIdsToRemove.add(0);
        }
        List<Peptide> newPepList = new ArrayList<>();
        for (int id = 0; id < pepList.size(); id++){
            if (pepIdsToRemove.contains(id)) continue;
            newPepList.add(pepList.get(id));
        }
        if (newPepList.isEmpty()) {
            return null;
        }

        Peptide[] peptideArray = newPepList.toArray(new Peptide[0]);
        Peptide topPep = peptideArray[0];

        Entry entry;
        String pepSetString = "";
        for (Peptide peptide : peptideArray){
            PeptideInfo peptideInfo = peptideInfoMap.get(peptide.getFreeSeq());
            pepSetString += peptide.getVarPtmContainingSeqNow() + "," + peptide.getScore() + "," + String.join("_", peptideInfo.protIdSet) +",";
        }

        boolean shouldPtm = Math.abs(precursorMass-massTool.calResidueMass(topPep.getFreeSeq()) - massTool.H2O) > 0.01;
        boolean hasPTM = topPep.hasVarPTM();
        int ptmNum = 0;
        boolean isSettled = true;
        double totalPtmMass = 0;
        if (hasPTM) {
            ptmNum = topPep.getVarPTMs().size();
            for (double mass : topPep.getVarPTMs().values()){
                totalPtmMass += mass;
            }
        }
        isSettled = Math.abs(totalPtmMass-(precursorMass-massTool.calResidueMass(topPep.getFreeSeq()) - massTool.H2O)) <= 0.01;

        double deltaLCn = 1; // L means the last?
        if (peptideArray.length > candisNum - 1) {
            deltaLCn = (peptideArray[0].getScore() - peptideArray[candisNum - 1].getScore()) / peptideArray[0].getScore();
        }
        double deltaCn = 1;
        if (peptideArray.length > 1) {
            for(int i = 0; i < peptideArray.length; i++) {
                if (peptideArray[i].getScore() != peptideArray[0].getScore()){
                    deltaCn = (peptideArray[0].getScore() - peptideArray[i].getScore()) / peptideArray[0].getScore();
                    break;
                }
            }
        }
        String otherPtmPatterns = "-";
        entry = new PreSearch.Entry(
                scanNum, scanName, shouldPtm ? 1 : 0, hasPTM ? 1 : 0, ptmNum, isSettled ? 1 : 0
                , precursorCharge, precursorMass, buildIndex.getLabelling(), topPep.getPtmContainingSeq(buildIndex.returnFixModMap())
                , topPep.getTheoMass(), topPep.isDecoy() ? 1 : 0, topPep.getGlobalRank()
                , topPep.getTagVecScore(), topPep.getScore(), deltaLCn, deltaCn
                , topPep.getMatchedPeakNum(), topPep.getIonFrac(), topPep.getMatchedHighestIntensityFrac()
                , topPep.getExplainedAaFrac(), otherPtmPatterns, topPep.getaScore(), ""
                , pepSetString.substring(0, pepSetString.length()-1)
        );
        entry.topPeptide = topPep;
        entry.topPeptide.precursorMass = precursorMass;
        for (Peptide peptide : peptideArray) {
            entry.peptideInfoMapForRef.put(peptide.getFreeSeq(), peptideInfoMap.get(peptide.getFreeSeq()));
        }
//        System.out.println(scanNum + ","+entry.peptideInfoMapForRef.size() + "," + peptideArray.length);
        if (lszDebugScanNum.contains(scanNum)){
            if (peptideInfoMap.containsKey(truth)) {
                int a = 1;
            }
        }
        int c = 1;


        return entry;
    }
    private boolean isHomo(Peptide p1, Peptide p2, Map<String, PeptideInfo> peptideInfoMap) {

        Set<String> temp1 = peptideInfoMap.get(p1.getFreeSeq()).protIdSet;
        Set<String> temp2 = peptideInfoMap.get(p2.getFreeSeq()).protIdSet;

        Set<String> set = new HashSet<>(temp1);
        set.retainAll(temp2);
//        if (set.isEmpty()) return false;
        SparseBooleanVector sbv1 = buildIndex.inferSegment.generateSegmentBooleanVector(p1.getFreeSeq());
        SparseBooleanVector sbv2 = buildIndex.inferSegment.generateSegmentBooleanVector(p2.getFreeSeq());
        return sbv1.dot(sbv2) > 0.3*Math.min(p1.getFreeSeq().length(), p2.getFreeSeq().length());
    }
    private boolean alreadySearched(ExpTag tagInfo, Set<ExpTag> usedTagInfoSet) {
        boolean res = false;
        for (ExpTag tmpTag : usedTagInfoSet) {
            if (tmpTag.size() < tagInfo.size()) continue;

            int alignedPos = tmpTag.getFreeAaString().indexOf(tagInfo.getFreeAaString());
            int endPos = alignedPos + tagInfo.size()-1;

            if (alignedPos >= 0
                    && Math.abs(tmpTag.expAaList.get(alignedPos).getHeadLocation() - tagInfo.getHeadLocation()) < 0.02){
                res = true;
                break;
            } else {
                int n_sharedAa = LCS(tagInfo.getFreeAaString(), tmpTag.getFreeAaString()).length();
                if (n_sharedAa > tagInfo.size()*0.7){
                    res = true;
                    break;
                }
            }
        }
        return res;
    }

//    private boolean sameTagStrUsed(String tagStr, Set<String> usedTagStrSet) {
//        return res;
//    }

    private static String LCS(String a, String b) {
        char[] arrayA = a.toCharArray();
        char[] arrayB = b.toCharArray();
        String maxString = "";
        for (int i = 0; i < arrayA.length; i++) {
            for (int j = 0; j < arrayB.length; j++) {
                if (arrayA[i]==arrayB[j]){
                    // i为起点
                    Integer i1 = i;
                    Integer j1 = j;
                    // 找到相同的字母后，以当前j为起点，持续比较下去，看相同部分到哪里结束
                    for (; i1 <arrayA.length && j1 < arrayB.length && arrayA[i1]==arrayB[j1]; i1++,j1++) {}
                    // i1为终点
                    if (i!=i1-1){
                        // 不断用common共同子串与max最长子串比较长度，只保留最长那个子串
                        String common = a.substring(i,i1);
                        if (maxString.length()<common.length()){
                            maxString = common;
                        }
                    }
                }
            }
        }
        return maxString;
    }
//    private int searchAndSave(ExpTag tagInfo, double minPcMass, double maxPcMass, Map<String, List<Peptide>> ptmPeptideListMap, Map<String,
//            List<Peptide>> freePeptideListMap, Map<String, PeptideInfo> peptideInfoMap, Set<String> peptidesFoundByThisTag, FMIndex fmIndex, char[] tagChar){
//        int numRes = 0;
//        SearchInterval searchRes = fmIndex.fmSearch(tagChar);
//        if (searchRes != null) {
//            numRes = searchRes.ep-searchRes.sp+1;
//            for (int ii = searchRes.sp; ii <= searchRes.ep; ii++) {
//                int absTagPos = fmIndex.SA[ii];
//                int dotIndex = -2-Arrays.binarySearch(buildIndex.dotPosArrReduced, absTagPos);
//                String protId = buildIndex.posProtMapReduced.get(dotIndex);
//                int relPos = absTagPos - buildIndex.dotPosArrReduced[dotIndex] - 1;
//
////                if (searchRes.settled) {
//                updateCandiList(protId, relPos, tagInfo, minPcMass, maxPcMass, ptmPeptideListMap, freePeptideListMap, peptideInfoMap, peptidesFoundByThisTag);
////                }
//            }
//        }
//        return numRes;
//    }

    private int searchAndSaveFuzzy(int scanNum, ExpTag tagInfo, double ms1TolAbs, Map<String, List<Peptide>> resPeptideListMap, Map<String, PeptideInfo> peptideInfoMap
            , Set<String> protIdSetByThisTag, FMIndex fmIndex, char[] tagChar, int minTagLen, SparseVector expProcessedPL,TreeMap<Double, Double> plMap, boolean isFirstTime) throws CloneNotSupportedException {

        int numRes = 0;
        SearchInterval searchRes = fmIndex.fmSearchFuzzy(tagChar);
        numRes = searchRes.ep-searchRes.sp+1;
//        if (numRes > 50) {
//            return numRes;
//        }
        int solCount = 0;
        if (searchRes.settled) {
            for (int ii = searchRes.sp; ii <= searchRes.ep; ii++) {
                if (solCount > 100) break; //At most take 300 sols
                int absTagPos = fmIndex.SA[ii];
                int dotIndex = -2-Arrays.binarySearch(buildIndex.dotPosArrReduced, absTagPos);

                String protId = buildIndex.posProtMapReduced.get(dotIndex);
                if (protIdSetByThisTag.contains(protId)){
                    continue;
                }
                solCount++;
                protIdSetByThisTag.add(protId);
                int relPos = absTagPos - buildIndex.dotPosArrReduced[dotIndex] - 1;
                updateCandiList(scanNum, protId, relPos, tagInfo, ms1TolAbs, resPeptideListMap, peptideInfoMap, expProcessedPL, plMap);
            }
        } else {
            if (!isFirstTime) return 0;
            int matchedPos = searchRes.matchedPos;
            if (tagInfo.size()-matchedPos < minTagLen) return 0;
            for (int ii = searchRes.sp; ii <= searchRes.ep; ii++) {
                if (solCount > 100) break;
                int absTagPos = fmIndex.SA[ii];
                int dotIndex = -2-Arrays.binarySearch(buildIndex.dotPosArrReduced, absTagPos);
                String protId = buildIndex.posProtMapReduced.get(dotIndex);
                String protSeq = buildIndex.protSeqMap.get(protId);
                int relPos = absTagPos - buildIndex.dotPosArrReduced[dotIndex] - 1;
                solCount++;
                updateCandiList(scanNum, protId, relPos, tagInfo.subTag(matchedPos, tagChar.length), ms1TolAbs, resPeptideListMap, peptideInfoMap, expProcessedPL, plMap);
            }

        }
//        System.out.println(scanNum + "," + tagInfo.getFreeAaString() + "," + numRes);

        return solCount;

    }

    private boolean isCMassValid(double mass) {
        double deltaMass = mass - precursorMass;
        if (deltaMass > minPtmMass - 530 && mass < maxPtmMass) {
            return true;
        }
        return false;
    }

    private boolean isNMassValid(double mass) {
        double deltaMass = mass - precursorMass;
        if (deltaMass > minPtmMass - 530 && mass < maxPtmMass) {
            return true;
        }
        return false;
    }
    private void updateCandiList(int scanNum, String protId, int tagPosInProt, ExpTag finderTag, double ms1TolAbs, Map<String, List<Peptide>> resPeptideListMap
            , Map<String, PeptideInfo> peptideInfoMap, SparseVector expProcessedPL, TreeMap<Double, Double> plMap) throws CloneNotSupportedException {

        double tagCMass = finderTag.getTailLocation() + massTool.H2O-MassTool.PROTON; // +massTool.H2O-MassTool.PROTON  is to mimic the mass of the real neutral precursor mass
        //tagCMass is correct even the tag is fuzzy. because FM starts searching from C
        String protSeq = buildIndex.protSeqMap.get(protId);
        Map<Integer, Double> cPoscMassMap = new HashMap<>();
        if ( (tagPosInProt+finderTag.size() == protSeq.length() && finderTag.isNorC != C_TAG)  // not C tag but found at prot Cterm, impossible
                || (tagPosInProt == 0 && finderTag.isNorC != N_TAG) // not N tag but found at prot Cterm, impossible
                || tagPosInProt < 0
                || tagPosInProt+finderTag.size() > protSeq.length()
        ) {
            return;
        }

        int missCleav = getNumOfMissCleavSite(finderTag.getFreeAaString());
        if (finderTag.isNorC == C_TAG || Math.abs(tagCMass - precursorMass) < ms1TolAbs){ // C tag, the end of tag must be KR and the tagCMass must be exactly settle. Otherwise no valid cPos will be, then no valid n pos.
            if (isKR(protSeq.charAt(tagPosInProt+finderTag.size()-1)) && Math.abs(tagCMass - precursorMass) < ms1TolAbs) { //fixme, should use MS1 tolerance not MS2 tolerance
                cPoscMassMap.put(tagPosInProt+finderTag.size()-1, tagCMass); // amino acid at cPos is also counted
            }
        } else { // only when oriTag is not C oriTag can it extend to c
            for (int i = tagPosInProt+finderTag.size(); i < protSeq.length(); i++) {  //
                //the if-block is in order, dont change the order.
                if (isX(protSeq.charAt(i))) break;
                tagCMass += massTool.getMassTable().get(protSeq.charAt(i)); //even when tag is fuzzy, tagcmass wont be disturbed
                if (tagCMass > precursorMass+maxPtmMass) break; // the total max minus ptm is -250. although the max minus ptm for single ptm is -156

                if (isKR(protSeq.charAt(i))) {
                    missCleav++; //current num of missed cleavage
                }
                if (tagCMass >= precursorMass-maxPtmMass && isKR(protSeq.charAt(i))) { // the total max plus ptm is 600. though the max plus ptm for single ptm is 527
                    cPoscMassMap.put(i, tagCMass);
                }
                if (missCleav > massTool.missedCleavage  && !isPtmSimuTest) { // if current num of KR is max, dont need to extend to c because it is impossible to take one more KR
                    break;         // stop extend to C term
                }
            }
        }

        int n_extension = 0;
        for (int cPos : cPoscMassMap.keySet()) {
            //check and find valid solution for cDeltaMass
            double cDeltaMass = precursorMass - cPoscMassMap.get(cPos);
            TreeMap<Integer, VarPtm> cPosVarPtmResMap = new TreeMap<>();
            char rightFlank;
            byte isProtNorC_Term = NON_TERM_PROT;
            if (cPos == protSeq.length()-1) {
                rightFlank = '-';
                isProtNorC_Term = C_TERM_PROT;
            } else {
                rightFlank = protSeq.charAt(cPos+1);
            }
            String cPartSeq = protSeq.substring(tagPosInProt+finderTag.size(), cPos+1);
            double cCutMass = finderTag.getTailLocation() - MassTool.PROTON;
            boolean isCTermFree = false;
            if (Math.abs(cDeltaMass) > 0.1) { // not ptm free in c part
                Set<Integer> fixModIdxes = inferPTM.getFixModIdxes(cPartSeq);
                Map<Integer, VarPtm[]> partPosVarModArrayMap = inferPTM.getIdxVarModMapNew(cPartSeq, fixModIdxes, C_PART, isProtNorC_Term); //todo no need to generate var mod list for aa again and again, make it stored.
                ModPepPool cPartPepsWithPtm = inferPTM.settlePtmOnSide(scanNum, expProcessedPL, plMap, precursorMass, cPartSeq, false,
                        partPosVarModArrayMap, cCutMass, cDeltaMass, precursorCharge, C_PART,ms1TolAbs);

                if (cPartPepsWithPtm.peptideTreeSet.isEmpty()) {
                    continue;
                }
                cPosVarPtmResMap = cPartPepsWithPtm.getTopPepPtn().posVarPtmResMap;
                if (cPosVarPtmResMap.isEmpty()) {
                    continue;  // c part should has PTM but unsettled.
                }
            } else {
                isCTermFree = true;
            }


            //extend and settle N part
            double nDeltaMass = finderTag.getHeadLocation() -MassTool.PROTON; //n term delta mass should be independent of cDeltaMass
            missCleav = getNumOfMissCleavSite(protSeq.substring(tagPosInProt, cPos)) ; //dont contain the c-term K, it does not count as miss cleav

            int min_nPos = (finderTag.isNorC == N_TAG && Math.abs(finderTag.getHeadLocation()-MassTool.PROTON) <= 0.02 ) ? tagPosInProt : 0;  //if it is from N tag and (when fuzzy) it is still N tag
            int max_nPos = (finderTag.isNorC == N_TAG && Math.abs(finderTag.getHeadLocation()-MassTool.PROTON) <= 0.02 ) ? tagPosInProt : tagPosInProt-1;
            for (int nPos = max_nPos; nPos >= min_nPos; nPos--) {
                if (isX(protSeq.charAt(nPos))) break;
                if (nPos < tagPosInProt){
                    nDeltaMass -= massTool.getMassTable().get(protSeq.charAt(nPos));
                }
                if (nDeltaMass < minPtmMass) break;
                if (nDeltaMass > maxPtmMass) continue;// the total max plus ptm is 600. though the max plus ptm for single ptm is 527

                if (nTermSpecific) {
                    if (nPos != 0 && !isKR(protSeq.charAt(nPos-1))) {// n term must be specific
                        continue;
                    }
                }

                if (isKR(protSeq.charAt(nPos))) {
                    missCleav++; //current num of missed cleavage
                }
                if (missCleav > massTool.missedCleavage  && !isPtmSimuTest) {
                    break;         // stop extend to n term
                }
                if (cPos+1-nPos < 6) { // min length of pep
                    continue;
                }

                // check N part
                String nPartSeq = protSeq.substring(nPos, tagPosInProt+1);
                char leftFlank;
                isProtNorC_Term = NON_TERM_PROT;
                if (nPos == 0 || (nPos == 1 && protSeq.charAt(0) == 'M')){
                    leftFlank = '-';
                    isProtNorC_Term = N_TERM_PROT;
                } else {
                    leftFlank = protSeq.charAt(nPos-1);
                }
                double nCutMass = finderTag.getTailLocation() - MassTool.PROTON;
                TreeMap<Integer, VarPtm> nPosVarPtmResMap = new TreeMap<>();
                boolean isNTermFree = false;
                if (Math.abs(nDeltaMass) > 0.1){
                    Set<Integer> fixModIdxes = inferPTM.getFixModIdxes(nPartSeq);
                    Map<Integer, VarPtm[]> partPosVarModArrayMap = inferPTM.getIdxVarModMapNew(nPartSeq, fixModIdxes, N_PART, isProtNorC_Term); //todo no need to generate var mod list for aa again and again, make it stored.
                    ModPepPool nPartModPepsSettled = inferPTM.settlePtmOnSide(scanNum, expProcessedPL, plMap, precursorMass, nPartSeq, false,
                            partPosVarModArrayMap, nCutMass, nDeltaMass, precursorCharge, N_PART,ms1TolAbs);

                    if (nPartModPepsSettled.peptideTreeSet.isEmpty()) {
//                        System.out.println(scanNum + "," + nPartSeq);
                        continue; // how could it be return!!!
                    }
                    nPosVarPtmResMap = nPartModPepsSettled.getTopPepPtn().posVarPtmResMap;
                    if (nPosVarPtmResMap.isEmpty()) {
                        continue;  // n part should has PTM but unsettled.
                    }
                } else {
                    isNTermFree = true;
                }
                n_extension++;
                int tagPosInPep = tagPosInProt-nPos;
                // store this peptide
                String freePepSeq = protSeq.substring(nPos, cPos+1);
                StringBuilder ptmPepSeqSB = new StringBuilder(freePepSeq);
                ptmPepSeqSB.replace(tagPosInPep, tagPosInPep+finderTag.size(), finderTag.getPtmAaString());

                Peptide peptide = new Peptide(freePepSeq, true, massTool, 1, 0.999, 0);// these paras are dummy answer will be deleted
                peptide.tagPosInPep = tagPosInPep;
                peptide.ptmSeq = ptmPepSeqSB.toString();
                peptide.finderTag = finderTag;
                peptide.cDeltaMass = cDeltaMass;
                peptide.nDeltaMass = nDeltaMass;

                PosMassMap fullPosMassMap = new PosMassMap(freePepSeq.length());
                TreeMap<Integer, VarPtm> posVarPtmResMap = new TreeMap<>();
                for (int pos : nPosVarPtmResMap.keySet()) { //copy the top 1 ptm pattern in n part // whhat if also choose the largeset priority
                    fullPosMassMap.put(pos, nPosVarPtmResMap.get(pos).mass); // copy the ptms from partModPepsUnsettled
                    posVarPtmResMap.put(pos, nPosVarPtmResMap.get(pos));
                }
                for (int pos : cPosVarPtmResMap.keySet()) { //copy the top 1 ptm pattern in n part // whhat if also choose the largeset priority
                    fullPosMassMap.put(pos+tagPosInPep+finderTag.size(), cPosVarPtmResMap.get(pos).mass); // copy the ptms from partModPepsUnsettled
                    posVarPtmResMap.put(pos+tagPosInPep+finderTag.size(), cPosVarPtmResMap.get(pos));
                }
                int idOfAa = -1;
                for (char aaChar : finderTag.getPtmAaString().toCharArray()) {
                    if (Character.isUpperCase(aaChar)) {
                        idOfAa += 1;
                    } else {
                        fullPosMassMap.put(tagPosInPep+idOfAa, massTool.labelVarPtmMap.get(aaChar).mass);
                        posVarPtmResMap.put(tagPosInPep+idOfAa, massTool.labelVarPtmMap.get(aaChar));
                    }
                }

                if (! fullPosMassMap.isEmpty()) { // has ptm in n , c or on the tag, it is a ptm-containing peptide
                    peptide.setVarPTM(fullPosMassMap);
                    peptide.posVarPtmResMap = posVarPtmResMap;
                }
                peptide.setScore(massTool.buildVectorAndCalXCorr(peptide.getIonMatrixNow(), 1, expProcessedPL, peptide.matchedBions,  peptide.matchedYions));//todo decide the penalty
//                peptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, fullPeptide.getIonMatrixNow(), ms2Tolerance));

                List<Peptide> peptideList = resPeptideListMap.get(freePepSeq);
                if (peptideList == null){
                    peptideList = new ArrayList<>();
                    peptideList.add(peptide);
                    resPeptideListMap.put(freePepSeq, peptideList);
                } else {
                    peptideList.add(peptide);
                }

                //update peptideInfoMapForRef,  just the flanks and proteins, these are standard info should be independent of peptide candidates.
                if (peptideInfoMap.containsKey(freePepSeq)) {
                    PeptideInfo peptideInfo = peptideInfoMap.get(freePepSeq);
                    if (!peptideInfo.protIdSet.contains(protId)) { //if this pep with prot is not recorded, add this prot
                        if (peptideInfo.leftFlank != '-' && peptideInfo.rightFlank != '-') {
                            if (rightFlank == '-' || leftFlank == '-') {
                                peptideInfo.leftFlank = leftFlank;
                                peptideInfo.rightFlank = rightFlank;
                            }
                        }
                        peptideInfo.protIdSet.add(protId);
                        if (!protId.startsWith("DECOY_")) {
                            peptideInfo.isDecoy = false;
                        }
                    }
                } else {
                    PeptideInfo peptideInfo = new PeptideInfo(freePepSeq, protId.startsWith("DECOY_"), leftFlank, rightFlank);
                    peptideInfo.protIdSet.add(protId);
                    peptideInfoMap.put(freePepSeq, peptideInfo);
                }

            }
        }
//        System.out.println(scanNum + ", n_extention," + finderTag.getFreeAaString() +","+ n_extension);

        return;
    }

    private int getNumOfMissCleavSite(String seq) {
        String str1 = seq.substring(0,seq.length());
        String str2 = str1.replaceAll("[KR]","");
        return str1.length()-str2.length();
    }
    private boolean isKR(char aa){
        return aa == 'K' || aa == 'R';
    }

    private boolean isX(char aa){
        return aa == 'X';
    }
    public class PresearchEntry {

        public Map<String, PeptideInfo> peptideInfoMapForRef = new HashMap<>();
        public double precursorMass = PreSearch.this.precursorMass;
        public int precursorCharge = PreSearch.this.precursorCharge;
        public List<Peptide> ptmCandiList = new ArrayList<>();
        public List<Peptide> freeCandiList = new ArrayList<>();

        public String scanName;
        PresearchEntry() {
        }
    }
    public class Entry { //copied from PtmSearch
        public Map<String, PeptideInfo> peptideInfoMapForRef = new HashMap<>();
        public Peptide topPeptide = null;
        final int scanNum;
        final String scanName;
        final int precursorCharge;
        final double precursorMass;
        final String labelling;
        public final String peptide;
        final double theoMass;
        final int isDecoy;
        final int globalRank;
        final double normalizedCorrelationCoefficient;
        public final double score;
        final double deltaLCn;
        final double deltaCn;
        final int matchedPeakNum;
        final double ionFrac;
        final double matchedHighestIntensityFrac;
        final double explainedAaFrac;
        final String otherPtmPatterns; // It has 4 decimal because it is write the the result file for checking. It is not used in scoring or other purpose.
        final String aScore;
        final String candidates;
        final String peptideSet;
        final int hasPTM;
        final int ptmNum;
        final int isSettled;
        final int shouldPtm;
        Entry(int scanNum, String scanName, int shouldPtm, int hasPTM, int ptmNum, int isSetteld, int precursorCharge, double precursorMass
                ,String labelling, String peptide, double theoMass, int isDecoy, int globalRank, double normalizedCorrelationCoefficient
                , double score, double deltaLCn, double deltaCn, int matchedPeakNum, double ionFrac, double matchedHighestIntensityFrac
                , double explainedAaFrac, String otherPtmPatterns, String aScore, String candidates, String peptideSet) {
            this.scanNum = scanNum;
            this.scanName = scanName;
            this.shouldPtm = shouldPtm;
            this.hasPTM = hasPTM;
            this.ptmNum = ptmNum;
            this.isSettled = isSetteld;
            this.precursorCharge = precursorCharge;
            this.precursorMass = precursorMass;
            this.labelling = labelling;
            this.peptide = peptide;
            this.theoMass = theoMass;
            this.isDecoy = isDecoy;
            this.globalRank = globalRank;
            this.normalizedCorrelationCoefficient = normalizedCorrelationCoefficient;
            this.score = score;
            this.deltaLCn = deltaLCn;
            this.deltaCn = deltaCn;
            this.matchedPeakNum = matchedPeakNum;
            this.ionFrac = ionFrac;
            this.matchedHighestIntensityFrac = matchedHighestIntensityFrac;
            this.explainedAaFrac = explainedAaFrac;
            this.otherPtmPatterns = otherPtmPatterns;
            this.aScore = aScore;
            this.candidates = candidates;
            this.peptideSet = peptideSet;
        }
    }
}
