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
        if (lszDebugScanNum.contains(this.scanNum)) {
            int a = 1;
        }

        double minPcMass = -1 * ms1Tolerance;
        double maxPcMass = ms1Tolerance;
        if (ms1ToleranceUnit == 1) {
            minPcMass = (precursorMass * leftInverseMs1Tolerance);
            maxPcMass = (precursorMass * rightInverseMs1Tolerance);
        }

        TreeMap<Double, Double> plMap = specProcessor.preSpectrumTopNStyleWithChargeLimit(rawPLMap, precursorMass, precursorCharge, minClear, maxClear, DatasetReader.topN);

        if (plMap.isEmpty()) return null;
        // Coding
        InferSegment inferSegment = buildIndex.getInferSegment();
        TreeMap<Double, Double> finalPlMap = inferSegment.addVirtualPeaks(precursorMass, plMap);

        List<ExpTag> allLongTagList = inferSegment.getLongTag(finalPlMap, precursorMass - massTool.H2O + MassTool.PROTON, scanNum, 4);

        if (allLongTagList.isEmpty())  return null;
        Entry entry = new Entry();


        double totalMass = precursorMass + 2 * MassTool.PROTON;

        FMIndex fmIndex = buildIndex.fmIndexReduced;

        Map<String, List<Peptide>> ptmPeptideListMap = new HashMap<>();
        Map<String, List<Peptide>> freePeptideListMap = new HashMap<>();
        Map<String, PeptideInfo> peptideInfoMap = new HashMap<>(50000);

        for (ExpTag tagInfo : allLongTagList){
            String tag = tagInfo.getFreeAaString();
            String revTagStr = new StringBuilder(tag).reverse().toString();

            SearchInterval searchForward = null;
            SearchInterval searchBackward = null;


            Set<String> peptidesFoundByThisTag = new HashSet<>();
            if (tagInfo.isNorC == -1) { //n tag
                char[] tagChar = tag.toCharArray();
                searchForward = fmIndex.fmSearch(tagChar);
                if (searchForward != null) {
//                    ptnForwardCount = searchForward.ep - searchForward.sp + 1;
                    for (int ii = searchForward.sp; ii <= searchForward.ep; ii++) {
                        int absTagPos = fmIndex.SA[ii];
                        int dotIndex = -2-Arrays.binarySearch(buildIndex.dotPosArrReduced, absTagPos);
                        // fmIndex.SA[ii]  the index of tag's first aa in the whole long proteins string
                        // -res-2  the dot sign index which the tag index is closed to
                        String protId = buildIndex.posProtMapReduced.get(dotIndex);
                        int relPos = absTagPos - buildIndex.dotPosArrReduced[dotIndex] - 1;

                        updateCandiList(protId, relPos, tagInfo, minPcMass, maxPcMass, ptmPeptideListMap, freePeptideListMap, peptideInfoMap, peptidesFoundByThisTag);
                    }
                }
            } else if (tagInfo.isNorC == 1) { // c tag
                char[] revTagChar = revTagStr.toCharArray();
                searchBackward = fmIndex.fmSearch(revTagChar);
                if (searchBackward != null) {
//                    ptnBackwardCount = searchBackward.ep - searchBackward.sp + 1;
                    for (int ii = searchBackward.sp; ii <= searchBackward.ep; ii++) {
                        int absTagPos = fmIndex.SA[ii];
                        int dotIndex = -2-Arrays.binarySearch(buildIndex.dotPosArrReduced, absTagPos);
                        String protId = buildIndex.posProtMapReduced.get(dotIndex);
                        int relPos = absTagPos - buildIndex.dotPosArrReduced[dotIndex] - 1;
                        if (relPos == 21478) {
                            int a = 1;
                        }
                        updateCandiList(protId, relPos, tagInfo.revTag(totalMass), minPcMass, maxPcMass, ptmPeptideListMap, freePeptideListMap, peptideInfoMap, peptidesFoundByThisTag); // put the reversed tag
                    }
                }
            } else { // non-nc tag
                char[] tagChar = tag.toCharArray();
                searchForward = fmIndex.fmSearch(tagChar);
                if (searchForward != null) {
                    for (int ii = searchForward.sp; ii <= searchForward.ep; ii++) {
                        int absTagPos = fmIndex.SA[ii];
                        int dotIndex = -2-Arrays.binarySearch(buildIndex.dotPosArrReduced, absTagPos);
                        String protId = buildIndex.posProtMapReduced.get(dotIndex);
                        int relPos = absTagPos - buildIndex.dotPosArrReduced[dotIndex] - 1;

                        updateCandiList(protId, relPos, tagInfo, minPcMass, maxPcMass, ptmPeptideListMap, freePeptideListMap, peptideInfoMap, peptidesFoundByThisTag);
                    }
                }

                char[] revTagChar = revTagStr.toCharArray();
                searchBackward = fmIndex.fmSearch(revTagChar);
                if (searchBackward != null) {
                    for (int ii = searchBackward.sp; ii <= searchBackward.ep; ii++) {
                        int absTagPos = fmIndex.SA[ii];
                        int dotIndex = -2-Arrays.binarySearch(buildIndex.dotPosArrReduced, absTagPos);
                        String protId = buildIndex.posProtMapReduced.get(dotIndex);
                        int relPos = absTagPos - buildIndex.dotPosArrReduced[dotIndex] - 1;

                        updateCandiList(protId, relPos, tagInfo.revTag(totalMass), minPcMass, maxPcMass, ptmPeptideListMap, freePeptideListMap, peptideInfoMap, peptidesFoundByThisTag); // put the reversed tag
                    }
                }
            }
        }

        entry.scanName = this.scanName;
//        List<ExpTag> allLongTagList = inferSegment.getLongTag(finalPlMap, precursorMass - massTool.H2O + MassTool.PROTON, scanNum, 4);

        Map<String, Double> scanTagStrMap = inferSegment.getTagStrMap(allLongTagList);
//        Search search = new Search(entry, scanNum, buildIndex.inferSegment, precursorMass, scanTagStrMap, massTool, localMaxMs2Charge, peptideInfoMap);
        if (lszDebugScanNum.contains(scanNum)){
            if (peptideInfoMap.containsKey(truth)) {
                int a = 1;
            }
        }

        List<Pair<String, Double>> ptmPeptideTotalScoreList = new LinkedList<>();
        List<Pair<String, Double>> freePeptideTotalScoreList = new LinkedList<>();
        for (String seq : ptmPeptideListMap.keySet()) {
            double totalScore = 0;
            for (Peptide peptide : ptmPeptideListMap.get(seq)) {
                totalScore += peptide.finderTag.getTotalIntensity();
            }
            ptmPeptideTotalScoreList.add(new Pair<>(seq, totalScore));
        }
        for (String seq : freePeptideListMap.keySet()) {
            double totalScore = 0;
            for (Peptide peptide : freePeptideListMap.get(seq)) {
                totalScore += peptide.finderTag.getTotalIntensity();
            }
            freePeptideTotalScoreList.add(new Pair<>(seq, totalScore));
        }
        Collections.sort(ptmPeptideTotalScoreList, Comparator.comparing(o -> o.getValue(), Comparator.reverseOrder()));
        Collections.sort(freePeptideTotalScoreList, Comparator.comparing(o -> o.getValue(), Comparator.reverseOrder()));
        for (int i = 0; i < Math.min(ptmPeptideTotalScoreList.size(), 10); i++) {
            String seq = ptmPeptideTotalScoreList.get(i).getKey();
            for (Peptide peptide : ptmPeptideListMap.get(seq)) {
                peptide.isDecoy = peptideInfoMap.get(seq).isDecoy; //override the fake isDecoy with true one
                entry.peptideInfoMapForRef.put(seq, peptideInfoMap.get(seq));
                entry.ptmCandiList.add(peptide);
            }
        }
        for (int i = 0; i < Math.min(freePeptideTotalScoreList.size(), 10); i++) {
            String seq = freePeptideTotalScoreList.get(i).getKey();
            for (Peptide peptide : freePeptideListMap.get(seq)) {
                peptide.isDecoy = peptideInfoMap.get(seq).isDecoy; //override the fake isDecoy with true one
                entry.peptideInfoMapForRef.put(seq, peptideInfoMap.get(seq));
                entry.freeCandiList.add(peptide);
            }
        }
        int c = 1;
        return entry;
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


        double tagCMass = tag.getTailLocation() + massTool.H2O-MassTool.PROTON; // +massTool.H2O-MassTool.PROTON  is to mimic the mass of the real neutral precursor mass
        String protSeq = buildIndex.protSeqMap.get(protId);
        Map<Integer, Double> cPoscMassMap = new HashMap<>();
        if (pos+tag.size() >= protSeq.length()) {
            int a = 1;
        }
        if (isKR(protSeq.charAt(pos+tag.size()-1)) && Math.abs(tagCMass - precursorMass) < 250) {
            cPoscMassMap.put(pos+tag.size()-1, tagCMass); // amino acid at cPos is also counted
        }

        int missCleav = getNumOfMissCleavSite(tag.getFreeAaString());
        if (tag.isNorC != 1) { // only when oriTag is not C oriTag can it extend to c
            for (int i = pos+tag.size(); i < protSeq.length(); i++) {  //

                //the if-block is in order, dont change the order.
                if (isX(protSeq.charAt(i))) break;
                tagCMass += massTool.getMassTable().get(protSeq.charAt(i));
                if (tagCMass > precursorMass+maxPtmMass) break; // the total max minus ptm is -250. although the max minus ptm for single ptm is -156

                if (isKR(protSeq.charAt(i))) {
                    missCleav++; //current num of missed cleavage
                }
                if (tagCMass >= precursorMass-600 && isKR(protSeq.charAt(i))) { // the total max plus ptm is 600. though the max plus ptm for single ptm is 527
                    cPoscMassMap.put(i, tagCMass);
                }
                if (missCleav > maxMissCleav) { // if current num of KR is max, dont need to extend to c because it is impossible to take one more KR
                    break;         // stop extend to C term
                }
            }
        }

        for (int cPos : cPoscMassMap.keySet()) {
            double nDeltaMass = precursorMass - cPoscMassMap.get(cPos) + tag.getHeadLocation() -MassTool.PROTON;
            char rightFlank;
            if (cPos == protSeq.length()-1) {
                rightFlank = '-';
            } else {
                rightFlank = protSeq.charAt(cPos+1);
            }
//                    int numMissCleave = 0; // todo miss cleavage
            missCleav = getNumOfMissCleavSite(protSeq.substring(pos, cPos)) ; //dont contain the c-term K, it does not count as miss cleav

            int min_nPos = (tag.isNorC == -1) ? pos : 0;
            int max_nPos = (tag.isNorC == -1) ? pos : pos-1;
            for (int nPos = max_nPos; nPos >= min_nPos; nPos--) {
                if (isX(protSeq.charAt(nPos))) break;
                if (nPos < pos){
                    nDeltaMass -= massTool.getMassTable().get(protSeq.charAt(nPos));
                }
                if (nDeltaMass < -250) break;
                if (nDeltaMass > 600) continue;// the total max plus ptm is 600. though the max plus ptm for single ptm is 527

                if (nTermSpecific) {
                    if (nPos != 0 && !isKR(protSeq.charAt(nPos-1))) {// n term must be specific
                        continue;
                    }
                }

                if (isKR(protSeq.charAt(nPos))) {
                    missCleav++; //current num of missed cleavage
                }
                if (missCleav > maxMissCleav) {
                    break;         // stop extend to n term
                }
                if (cPos+1-nPos < 6) {
                    continue;
                }
                String freePepSeq = protSeq.substring(nPos, cPos+1);

                if (freePepSeq.contentEquals("TNVNPSEVGDLVVGSVLAPGAQR")) {
                    int a = 1;
                }
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

                Peptide peptide = new Peptide(freePepSeq, true, massTool, 1, 0.999, 0);
                // these paras are dummy answer will be deleted
                peptide.tagPosInPep = pos-nPos;
                peptide.ptmSeq = ptmPepSeqSB.toString();
                peptide.finderTag = tag;
                peptide.cDeltaMass = precursorMass - cPoscMassMap.get(cPos);
                peptide.nDeltaMass = nDeltaMass;
                if (theoTotalMassWithPtm > minPcMass && theoTotalMassWithPtm < maxPcMass){
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
