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
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.FM.FMIndex;
import proteomics.FM.SearchInterval;
import proteomics.Index.BuildIndex;
import proteomics.Segment.InferSegment;
import proteomics.Spectrum.DatasetReader;
import proteomics.Types.*;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;

import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.locks.ReentrantLock;

import static proteomics.PIPI.*;
import static proteomics.Segment.InferSegment.C_TAG;
import static proteomics.Segment.InferSegment.N_TAG;

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
        if (lszDebugScanNum.contains(this.scanNum)) {
            System.out.println(scanNum + ", entered");
            int a = 1;
        }

        double minPcMass = -1 * ms1Tolerance;
        double maxPcMass = ms1Tolerance;
        if (ms1ToleranceUnit == 1) {
            minPcMass = (precursorMass * leftInverseMs1Tolerance);
            maxPcMass = (precursorMass * rightInverseMs1Tolerance);
        }

        TreeMap<Double, Double> plMap = specProcessor.preSpectrumTopNStyleWithChargeLimit(rawPLMap, precursorMass, precursorCharge, minClear, maxClear, DatasetReader.topN, ms2Tolerance);

        if (plMap.isEmpty()) return null;
        // Coding
        InferSegment inferSegment = buildIndex.getInferSegment();
        TreeMap<Double, Double> finalPlMap = inferSegment.addVirtualPeaks(precursorMass, plMap);

        List<ExpTag> allLongTagList = inferSegment.getLongTag(finalPlMap, precursorMass - massTool.H2O + MassTool.PROTON, scanNum, minTagLenToExtract,maxTagLenToExtract);

        if (allLongTagList.isEmpty())  return null;
        Entry entry = new Entry();


        double totalMass = precursorMass + 2 * MassTool.PROTON;

        FMIndex fmIndex = buildIndex.fmIndexReduced;

        Map<String, List<Peptide>> ptmPeptideListMap = new HashMap<>();
        Map<String, List<Peptide>> freePeptideListMap = new HashMap<>();
        Map<String, PeptideInfo> peptideInfoMap = new HashMap<>(50000);

        Set<String> searchedTagStrSet = new HashSet<>();
        for (ExpTag tagInfo : allLongTagList){
            String tagStr = tagInfo.getFreeAaString();
            String revTagStr = new StringBuilder(tagStr).reverse().toString();
            Set<String> pepStringFoundByThisTag = new HashSet<>();

            if (tagInfo.isNorC == N_TAG) { //n tag
                if (!searchedTagStrSet.contains(tagStr)) {
                    char[] tagChar = tagStr.toCharArray();
                    int numRes = searchAndSaveFuzzy(tagInfo, minPcMass, maxPcMass, ptmPeptideListMap, freePeptideListMap, peptideInfoMap, pepStringFoundByThisTag, fmIndex, tagChar);
                    searchedTagStrSet.add(tagStr);
                    if (numRes == 0 && tagStr.length() > 6){ // if the tag was already long
                        for (int subTime = 1; subTime <= 3; subTime++){
                            char[] subTagChar = tagStr.substring(0, tagStr.length()-subTime).toCharArray();
                            ExpTag subTagInfo = tagInfo.subTag(0,tagStr.length()-subTime);
                            if (!searchedTagStrSet.contains(tagStr)) {
                                int numResSub = searchAndSaveFuzzy(subTagInfo, minPcMass, maxPcMass, ptmPeptideListMap, freePeptideListMap, peptideInfoMap, pepStringFoundByThisTag, fmIndex, subTagChar);
                                searchedTagStrSet.add(tagStr);
                                if (numResSub > 0) break;
                            }
                        }
                    }
                }
            } else if (tagInfo.isNorC == C_TAG) { // c tag
                char[] revTagChar = revTagStr.toCharArray();
                ExpTag revTagInfo = tagInfo.revTag(totalMass);
                if (!searchedTagStrSet.contains(revTagStr)) {
                    int numRes = searchAndSaveFuzzy(revTagInfo, minPcMass, maxPcMass, ptmPeptideListMap, freePeptideListMap, peptideInfoMap, pepStringFoundByThisTag, fmIndex, revTagChar);
                    searchedTagStrSet.add(revTagStr);
                    if (numRes == 0 && tagStr.length() > 6){ // if the tag was already long
                        for (int subTime = 1; subTime <= 3; subTime++){
                            String subRevTagStr = revTagStr.substring(0, tagStr.length()-subTime);
                            char[] subRevTagChar = subRevTagStr.toCharArray();
                            ExpTag subRevTagInfo = revTagInfo.subTag(0,tagStr.length()-subTime);
                            if (!searchedTagStrSet.contains(subRevTagStr)) {
                                int numResSub = searchAndSaveFuzzy(subRevTagInfo, minPcMass, maxPcMass, ptmPeptideListMap, freePeptideListMap, peptideInfoMap, pepStringFoundByThisTag, fmIndex, subRevTagChar);
                                searchedTagStrSet.add(subRevTagStr);
                                if (numResSub > 0) break;
                            }
                        }
                    }
                }
            } else { // non-nc tag
                char[] tagChar = tagStr.toCharArray();
                if (!searchedTagStrSet.contains(tagStr)) {
                    int numResForward = searchAndSaveFuzzy(tagInfo, minPcMass, maxPcMass, ptmPeptideListMap, freePeptideListMap, peptideInfoMap, pepStringFoundByThisTag, fmIndex, tagChar);
                    searchedTagStrSet.add(tagStr);
                    if (numResForward == 0 && tagStr.length() > 6){ // if the tag was already long
                        for (int subTime = 1; subTime <= 3; subTime++){
                            String subTagStr = tagStr.substring(0, tagStr.length()-subTime);
                            ExpTag subTagInfo = tagInfo.subTag(0,tagStr.length()-subTime);

                            String revSubTagStr = (new StringBuilder(subTagStr)).reverse().toString();
                            //sub forward
                            char[] subTagChar = subTagStr.toCharArray();
                            int numResSub1 = 0;
                            if (!searchedTagStrSet.contains(subTagStr)) {
                                numResSub1 = searchAndSaveFuzzy(subTagInfo, minPcMass, maxPcMass, ptmPeptideListMap, freePeptideListMap, peptideInfoMap, pepStringFoundByThisTag, fmIndex, subTagChar);
                                searchedTagStrSet.add(subTagStr);
                            }

                            //sub backward
                            char[] revSubTagChar = revSubTagStr.toCharArray();
                            int numResSub2 = 0;
                            ExpTag revSubTagInfo = subTagInfo.revTag(totalMass);
                            if (!searchedTagStrSet.contains(revSubTagStr)) {
                                numResSub2 = searchAndSaveFuzzy(revSubTagInfo, minPcMass, maxPcMass, ptmPeptideListMap, freePeptideListMap, peptideInfoMap, pepStringFoundByThisTag, fmIndex, revSubTagChar);
                                searchedTagStrSet.add(revSubTagStr);
                            }
                            if (numResSub1 + numResSub2 > 0) break;
                        }
                    }
                }
                char[] revTagChar = revTagStr.toCharArray();
                ExpTag revTagInfo = tagInfo.revTag(totalMass);
                if (!searchedTagStrSet.contains(revTagStr)) {
                    int numResBackward = searchAndSaveFuzzy(revTagInfo, minPcMass, maxPcMass, ptmPeptideListMap, freePeptideListMap, peptideInfoMap, pepStringFoundByThisTag, fmIndex, revTagChar);
                    searchedTagStrSet.add(revTagStr);
                    if (numResBackward == 0 && revTagStr.length() > 6){ // if the tag was already long
                        for (int subTime = 1; subTime <= 3; subTime++){
                            String subTagStr = revTagStr.substring(0, revTagStr.length()-subTime);
                            ExpTag subTagInfo = tagInfo.revTag(totalMass).subTag(0,revTagStr.length()-subTime);

                            String revSubTagStr = (new StringBuilder(subTagStr)).reverse().toString();
                            //sub forward
                            char[] subTagChar = subTagStr.toCharArray();
                            int numResSub1 = 0;
                            if (!searchedTagStrSet.contains(subTagStr)) {
                                numResSub1 = searchAndSaveFuzzy(subTagInfo, minPcMass, maxPcMass, ptmPeptideListMap, freePeptideListMap, peptideInfoMap, pepStringFoundByThisTag, fmIndex, subTagChar);
                                searchedTagStrSet.add(subTagStr);
                            }
                            //sub backward
                            char[] revSubTagChar = revSubTagStr.toCharArray();
                            ExpTag revSubTagInfo = subTagInfo.revTag(totalMass);
                            int numResSub2 = 0;
                            if (!searchedTagStrSet.contains(revSubTagStr)) {
                                numResSub2 = searchAndSaveFuzzy(revSubTagInfo, minPcMass, maxPcMass, ptmPeptideListMap, freePeptideListMap, peptideInfoMap, pepStringFoundByThisTag, fmIndex, revSubTagChar);
                                searchedTagStrSet.add(revSubTagStr);
                            }
                            if (numResSub1 + numResSub2 > 0) break;
                        }
                    }
                }
            }
        }

        entry.scanName = this.scanName;
//        List<ExpTag> allLongTagList = inferSegment.getLongTag(finalPlMap, precursorMass - massTool.H2O + MassTool.PROTON, scanNum, 4);

        Map<String, Double> scanTagStrMap = inferSegment.getTagStrMap(allLongTagList);
//        Search search = new Search(entry, scanNum, buildIndex.inferSegment, precursorMass, scanTagStrMap, massTool, localMaxMs2Charge, peptideInfoMap);


        List<Pair<String, Double>> ptmPeptideTotalScoreList = new LinkedList<>();
        List<Pair<String, Double>> freePeptideTotalScoreList = new LinkedList<>();
        for (String seq : ptmPeptideListMap.keySet()) {
            double maxScore = 0;
            for (Peptide peptide : ptmPeptideListMap.get(seq)) {
                maxScore = Math.max(peptide.finderTag.getTotalIntensity(), maxScore);
            }
            ptmPeptideTotalScoreList.add(new Pair<>(seq, maxScore));
        }
        for (String seq : freePeptideListMap.keySet()) {
            double maxScore = 0;
            for (Peptide peptide : freePeptideListMap.get(seq)) {
                maxScore = Math.max(peptide.finderTag.getTotalIntensity(), maxScore);
            }
            freePeptideTotalScoreList.add(new Pair<>(seq, maxScore));
        }
        Collections.sort(ptmPeptideTotalScoreList, Comparator.comparing(o -> o.getValue(), Comparator.reverseOrder()));
        Collections.sort(freePeptideTotalScoreList, Comparator.comparing(o -> o.getValue(), Comparator.reverseOrder()));

//        cleanRedundantPeps = 0
        for (int i = 0; i < Math.min(ptmPeptideTotalScoreList.size(), 10); i++) {
            String seq = ptmPeptideTotalScoreList.get(i).getKey();
            List<Peptide> ptmPeptideList = ptmPeptideListMap.get(seq);
            Collections.sort(ptmPeptideList, Comparator.comparing(o -> o.finderTag.size(), Comparator.reverseOrder())); // make sure that ABCDE is in front of BCDE, to remove redundant BCDE
            Set<ExpTag> localUseTagInfoSet = new HashSet<>();
            for (Peptide peptide : ptmPeptideList) {
                if (!alreadySearched(peptide.finderTag, localUseTagInfoSet)) {
                    peptide.isDecoy = peptideInfoMap.get(seq).isDecoy; //override the fake isDecoy with true one
                    entry.peptideInfoMapForRef.put(seq, peptideInfoMap.get(seq));
                    entry.ptmCandiList.add(peptide);
                    localUseTagInfoSet.add(peptide.finderTag);
                }
            }
        }
        for (int i = 0; i < Math.min(freePeptideTotalScoreList.size(), 10); i++) {
            String seq = freePeptideTotalScoreList.get(i).getKey();
            List<Peptide> freePeptideList = freePeptideListMap.get(seq);
            Collections.sort(freePeptideList, Comparator.comparing(o -> o.finderTag.size(), Comparator.reverseOrder())); // make sure that ABCDE is in front of BCDE, to remove redundant BCDE
            Set<ExpTag> localUseTagInfoSet = new HashSet<>();
            for (Peptide peptide : freePeptideListMap.get(seq)) {
                if (!alreadySearched(peptide.finderTag, localUseTagInfoSet)) {
                    peptide.isDecoy = peptideInfoMap.get(seq).isDecoy; //override the fake isDecoy with true one
                    entry.peptideInfoMapForRef.put(seq, peptideInfoMap.get(seq));
                    entry.freeCandiList.add(peptide);
                    localUseTagInfoSet.add(peptide.finderTag);
                }
            }
        }
        if (lszDebugScanNum.contains(scanNum)){
            for (Peptide pep : entry.ptmCandiList){
                System.out.println(scanNum + ","+scanName + "," + pep.getFreeSeq()+ "," +  pep.finderTag.getFreeAaString()+ "," + pep.nDeltaMass+ "," + pep.cDeltaMass);
            }
            if (peptideInfoMap.containsKey(truth)) {
                int a = 1;
            }
        }
        int c = 1;
        return entry;
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
    private int searchAndSave(ExpTag tagInfo, double minPcMass, double maxPcMass, Map<String, List<Peptide>> ptmPeptideListMap, Map<String,
            List<Peptide>> freePeptideListMap, Map<String, PeptideInfo> peptideInfoMap, Set<String> peptidesFoundByThisTag, FMIndex fmIndex, char[] tagChar){
        int numRes = 0;
        SearchInterval searchRes = fmIndex.fmSearch(tagChar);
        if (searchRes != null) {
            numRes = searchRes.ep-searchRes.sp+1;
            for (int ii = searchRes.sp; ii <= searchRes.ep; ii++) {
                int absTagPos = fmIndex.SA[ii];
                int dotIndex = -2-Arrays.binarySearch(buildIndex.dotPosArrReduced, absTagPos);
                String protId = buildIndex.posProtMapReduced.get(dotIndex);
                int relPos = absTagPos - buildIndex.dotPosArrReduced[dotIndex] - 1;

//                if (searchRes.settled) {
                updateCandiList(protId, relPos, tagInfo, minPcMass, maxPcMass, ptmPeptideListMap, freePeptideListMap, peptideInfoMap, peptidesFoundByThisTag);
//                }
            }
        }
        return numRes;
    }

    private int searchAndSaveFuzzy(ExpTag tagInfo, double minPcMass, double maxPcMass, Map<String, List<Peptide>> ptmPeptideListMap, Map<String,
            List<Peptide>> freePeptideListMap, Map<String, PeptideInfo> peptideInfoMap, Set<String> peptidesFoundByThisTag, FMIndex fmIndex, char[] tagChar){

//        if (tagInfo.getFreeAaString().contentEquals("HFAHAAAELR") || tagInfo.getFreeAaString().contentEquals("AHAAAELR")) {
//            int a = 1;
//        }
        int numRes = 0;
        SearchInterval searchRes = fmIndex.fmSearchFuzzy(tagChar);
        if (searchRes.settled) {
            numRes = searchRes.ep-searchRes.sp+1;
            for (int ii = searchRes.sp; ii <= searchRes.ep; ii++) {
                int absTagPos = fmIndex.SA[ii];
                int dotIndex = -2-Arrays.binarySearch(buildIndex.dotPosArrReduced, absTagPos);
                String protId = buildIndex.posProtMapReduced.get(dotIndex);
                int relPos = absTagPos - buildIndex.dotPosArrReduced[dotIndex] - 1;
                updateCandiList(protId, relPos, tagInfo, minPcMass, maxPcMass, ptmPeptideListMap, freePeptideListMap, peptideInfoMap, peptidesFoundByThisTag);
//                usedTagInfoSet.add(tagInfo);
            }
        } else {
            if (tagInfo.size() < 6) {
                return numRes; // if an unsettled tag is shorter than 6, no meaning for fuzzy
            }
            int matchedPos = searchRes.matchedPos;
            for (int ii = searchRes.sp; ii <= searchRes.ep; ii++) {
                int absTagPos = fmIndex.SA[ii];
                int dotIndex = -2-Arrays.binarySearch(buildIndex.dotPosArrReduced, absTagPos);
                String protId = buildIndex.posProtMapReduced.get(dotIndex);
                String protSeq = buildIndex.protSeqMap.get(protId);
                int relPos = absTagPos - buildIndex.dotPosArrReduced[dotIndex] - 1;
                if (relPos < matchedPos) {
                    continue; // there is not enough aa in the front part of the prot seq, impossible to match the tag
                }

                int n_ExtraAa = 0;
                int startPos = matchedPos;
                if (matchedPos >= 2) {  //when matchedPos == 1, there is 1 wrong aa and 0 missed aa, when matchedPos == 0, searchInterval will be settled.
                    for (int missedAaPos = matchedPos-2; missedAaPos >= 0; missedAaPos--) {
                        if (tagChar[missedAaPos] != protSeq.charAt(relPos-matchedPos+missedAaPos)) {
                            startPos = missedAaPos+1; //last pos
                            break;
                        }
                        n_ExtraAa++;
                    }
                }

                int correctedRelPos;
                if (n_ExtraAa == 0) {
                    startPos = matchedPos; // if none extra aa is extended, reset startPos to matchedPos, just ignore the wrong and missed part
                    correctedRelPos = relPos;
                } else {
                    startPos = matchedPos-n_ExtraAa-1;
                    correctedRelPos = relPos-n_ExtraAa-1;
                }

                if (startPos < tagChar.length/2) {// the unsettled tag must be at least half of the original tag.
                    updateCandiList(protId, correctedRelPos, tagInfo.subTag(startPos, tagChar.length), minPcMass, maxPcMass, ptmPeptideListMap, freePeptideListMap, peptideInfoMap, peptidesFoundByThisTag);
                    numRes++;
//                        usedTagInfoSet.add(tagInfo.subTag(startPos, tagChar.length));
                }
            }
        }
        return numRes;
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
    private void updateCandiList(String protId, int pos, ExpTag tag, double minPcMass, double maxPcMass
            , Map<String, List<Peptide>> ptmPeptideListMap, Map<String, List<Peptide>> freePeptideListMap, Map<String, PeptideInfo> peptideInfoMap, Set<String> peptidesFoundByThisTag) {

//        if (tag.getFreeAaString().contentEquals("LPNLN")
//            && (protId.contentEquals("XP_003519282.1") || protId.contentEquals("NP_001241115.1"))
//            && scanNum == 72611) {
//            int a = 1;
//        }
        double tagCMass = tag.getTailLocation() + massTool.H2O-MassTool.PROTON; // +massTool.H2O-MassTool.PROTON  is to mimic the mass of the real neutral precursor mass
        //tagCMass is correct even the tag is fuzzy. because FM starts searching from C
        String protSeq = buildIndex.protSeqMap.get(protId);
        Map<Integer, Double> cPoscMassMap = new HashMap<>();
        if ( (pos+tag.size() >= protSeq.length() && tag.isNorC != C_TAG)  // not C tag but found at prot Cterm, impossible
                || (pos == 0 && tag.isNorC != N_TAG) // not N tag but found at prot Cterm, impossible
        ) {
            return;
        }

        int missCleav = getNumOfMissCleavSite(tag.getFreeAaString());
        if (tag.isNorC == C_TAG || Math.abs(tagCMass - precursorMass) < 0.02){ // C tag, the end of tag must be KR and the tagCMass must be exactly settle. Otherwise no valid cPos will be, then no valid n pos.
            if (isKR(protSeq.charAt(pos+tag.size()-1)) && Math.abs(tagCMass - precursorMass) < 0.02) { //fixme, should use MS1 tolerance not MS2 tolerance
                cPoscMassMap.put(pos+tag.size()-1, tagCMass); // amino acid at cPos is also counted
            }
        } else { // only when oriTag is not C oriTag can it extend to c
            for (int i = pos+tag.size(); i < protSeq.length(); i++) {  //
                //the if-block is in order, dont change the order.
                if (isX(protSeq.charAt(i))) break;
                tagCMass += massTool.getMassTable().get(protSeq.charAt(i)); //even when tag is fuzzy, tagcmass wont be disturbed
                if (tagCMass > precursorMass+maxPtmMass) break; // the total max minus ptm is -250. although the max minus ptm for single ptm is -156

                if (isKR(protSeq.charAt(i))) {
                    missCleav++; //current num of missed cleavage
                }
                if (tagCMass >= precursorMass-250 && isKR(protSeq.charAt(i))) { // the total max plus ptm is 600. though the max plus ptm for single ptm is 527
                    cPoscMassMap.put(i, tagCMass);
                }
                if (missCleav > massTool.missedCleavage  && !isPtmSimuTest) { // if current num of KR is max, dont need to extend to c because it is impossible to take one more KR
                    break;         // stop extend to C term
                }
            }
        }

        for (int cPos : cPoscMassMap.keySet()) {
            double nDeltaMass = tag.getHeadLocation() -MassTool.PROTON; //n term delta mass should be independent of cDeltaMass
            char rightFlank;
            if (cPos == protSeq.length()-1) {
                rightFlank = '-';
            } else {
                rightFlank = protSeq.charAt(cPos+1);
            }
//                    int numMissCleave = 0; // todo miss cleavage
            missCleav = getNumOfMissCleavSite(protSeq.substring(pos, cPos)) ; //dont contain the c-term K, it does not count as miss cleav

            int min_nPos = (tag.isNorC == N_TAG && Math.abs(tag.getHeadLocation()-MassTool.PROTON) <= 0.02 ) ? pos : 0;  //if it is from N tag and (when fuzzy) it is still N tag
            int max_nPos = (tag.isNorC == N_TAG && Math.abs(tag.getHeadLocation()-MassTool.PROTON) <= 0.02 ) ? pos : pos-1;
            for (int nPos = max_nPos; nPos >= min_nPos; nPos--) {
                if (isX(protSeq.charAt(nPos))) break;
                if (nPos < pos){
                    nDeltaMass -= massTool.getMassTable().get(protSeq.charAt(nPos));
                }
                if (nDeltaMass < -250) break;
                if (nDeltaMass > 250) continue;// the total max plus ptm is 600. though the max plus ptm for single ptm is 527

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
                String freePepSeq = protSeq.substring(nPos, cPos+1);


                StringBuilder ptmPepSeqSB = new StringBuilder(freePepSeq);
                ptmPepSeqSB.replace(pos-nPos, pos-nPos+tag.size(), tag.getPtmAaString());

                char leftFlank;
                if (nPos == 0 || (nPos == 1 && protSeq.charAt(0) == 'M')){
                    leftFlank = '-';
                } else {
                    leftFlank = protSeq.charAt(nPos-1);
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

                double theoTotalMassWithPtm = massTool.calLabeledSeqMass(ptmPepSeqSB.toString()) + massTool.H2O;
                if (peptidesFoundByThisTag.contains(freePepSeq)){
                    continue;
                } else {
                    peptidesFoundByThisTag.add(freePepSeq);
                }

                if (freePepSeq.contentEquals("DRCHFAHGAAELR")){
                    int a = 1;
                }
                Peptide peptide = new Peptide(freePepSeq, true, massTool, 1, 0.999, 0);
                // these paras are dummy answer will be deleted
                peptide.tagPosInPep = pos-nPos;
                peptide.ptmSeq = ptmPepSeqSB.toString();
                peptide.finderTag = tag;
                peptide.cDeltaMass = precursorMass - cPoscMassMap.get(cPos);
                peptide.nDeltaMass = nDeltaMass;
//                if (theoTotalMassWithPtm > minPcMass && theoTotalMassWithPtm < maxPcMass){// fixme, should it be cDeltaMass==0 and nDeltaMass==0 ?
                if (Math.abs(peptide.cDeltaMass) <= 0.1 && Math.abs(peptide.nDeltaMass) <= 0.1){// fixme, should it be cDeltaMass==0 and nDeltaMass==0 ?

                        //free
                    List<Peptide> peptideList = freePeptideListMap.get(freePepSeq);
                    if (peptideList == null){
                        peptideList = new ArrayList<>();
                        peptideList.add(peptide);
                        freePeptideListMap.put(freePepSeq, peptideList);
                    } else {
                        peptideList.add(peptide);
                    }
                } else {
                    //ptm
                    List<Peptide> peptideList = ptmPeptideListMap.get(freePepSeq);
                    if (peptideList == null){
                        peptideList = new ArrayList<>();
                        peptideList.add(peptide);
                        ptmPeptideListMap.put(freePepSeq, peptideList);
                    } else {
                        peptideList.add(peptide);
                    }
                }
            }
        }
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
    public class Entry {

        public Map<String, PeptideInfo> peptideInfoMapForRef = new HashMap<>();
        public double precursorMass = PreSearch.this.precursorMass;
        public int precursorCharge = PreSearch.this.precursorCharge;
        public List<Peptide> ptmCandiList = new ArrayList<>();
        public List<Peptide> freeCandiList = new ArrayList<>();

        public String scanName;
        Entry() {
        }
    }
}
