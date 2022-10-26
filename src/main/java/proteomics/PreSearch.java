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
import proteomics.Index.BuildIndex;
import proteomics.Segment.InferSegment;
import proteomics.Spectrum.DatasetReader;
import proteomics.Types.*;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;

import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.locks.ReentrantLock;

import static proteomics.PIPI.lszDebugScanNum;

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
        if (lszDebugScanNum.contains(this.scanNum)) {
            int a = 1;
        }
        TreeMap<Double, Double> plMap = specProcessor.preSpectrumTopNStyleWithChargeLimit(rawPLMap, precursorMass, precursorCharge, minClear, maxClear, DatasetReader.topN);

        if (plMap.isEmpty()) {
            return null;
        }

        // Coding
        InferSegment inferSegment = buildIndex.getInferSegment();
        TreeMap<Double, Double> finalPlMap = inferSegment.addVirtualPeaks(precursorMass, plMap);

        List<ThreeExpAA> tag4List = inferSegment.getAllTag5(precursorMass, finalPlMap, scanNum);
        Collections.sort(tag4List, Comparator.comparingDouble(ThreeExpAA::getTotalIntensity));
        Map<String, Set<Pair<String, Integer>>> tagProtPosMap = buildIndex.tagProtPosMap;
        Map<String, PeptideInfo> peptideInfoMap = new HashMap<>();

        if (tag4List.isEmpty())  return null;
        Entry entry = new Entry();

        double totalMass = precursorMass + 2 * MassTool.PROTON;
        List<ThreeExpAA> compTag4List = new ArrayList<>(2* tag4List.size());
        compTag4List.addAll(tag4List);
        for (ThreeExpAA tag4 : tag4List) {
            compTag4List.add(tag4.revTag(totalMass));
        }   // duplicate reverse tags for every tag, then dont need to care the both direction of a tag because already contained in this compTag4List

        double pcMassL = precursorMass - 250;
        double pcMassR = precursorMass + 250;
        for (ThreeExpAA tagInfo : compTag4List) {
            String tag = tagInfo.getPtmFreeAAString();

            if (!tagProtPosMap.containsKey(tagInfo.getPtmFreeAAString())) {
                int a = 1;
                continue;
            }
            Set<Pair<String, Integer>> protPosPairs = tagProtPosMap.get(tagInfo.getPtmFreeAAString());
            for (Pair<String, Integer> protPos : protPosPairs){
                String protId = protPos.getFirst();
                if (protId.contentEquals("sp|P35637|FUS_HUMAN") && lszDebugScanNum.contains(scanNum)) {
                    int a = 1;
                }
                int pos = protPos.getSecond();
                String protSeq = buildIndex.protSeqMap.get(protId);
                Set<Integer> cPosSet = new HashSet<>();

                double tagNMass = tagInfo.getHeadLocation();
                double tagCMass = tagInfo.getTailLocation();
                if ("KR".contains(protSeq.substring(pos+3, pos+4)) && tagCMass <= pcMassR && tagCMass  >= pcMassL) {
                    cPosSet.add(pos+3);
                }
                for (int i = pos+tag.length(); i < protSeq.length(); i++) {  //
                    tagCMass += massTool.calResidueMass(protSeq.substring(i, i+1));
                    if (tagCMass < pcMassL) continue;
                    if (tagCMass > pcMassR) break;
                    if ("KR".contains(protSeq.substring(i, i+1))) {
                        cPosSet.add(i);
                    }
                }
                for (int cPos : cPosSet) {
                    double curMass = massTool.calResidueMass(protSeq.substring(pos, cPos+1)) + massTool.H2O;
                    int numMissCleave = 0;
                    for (int nPos = pos-1; nPos > 0; nPos--) {
                        curMass += massTool.calResidueMass(protSeq.substring(nPos, nPos+1));
                        if (curMass < pcMassL) continue;
                        if (curMass > pcMassR) break;
                        String pepSeq = protSeq.substring(nPos, cPos+1);
                        if (peptideInfoMap.containsKey(pepSeq)) {
                            PeptideInfo pepInfo = peptideInfoMap.get(pepSeq);
                            pepInfo.protIdSet.add(protId);
                            if (!protId.startsWith("DECOY_")) {
                                pepInfo.isTarget = true;
                            }
                        } else {
                            char leftFlank = '-';
                            char rightFlank = '-';
                            if (nPos-1 >= 0) leftFlank = protSeq.charAt(nPos-1);
                            if (cPos+1 <= protSeq.length()-1) rightFlank = protSeq.charAt(cPos+1);

                            PeptideInfo pepInfo = new PeptideInfo(pepSeq, !protId.startsWith("DECOY_"), leftFlank, rightFlank);
                            pepInfo.protIdSet.add(protId);
                            peptideInfoMap.put(pepSeq, pepInfo);
                        }
                    }
                }
            }
        }
        System.out.println(scanNum + ","+ compTag4List.size()+"," + peptideInfoMap.size());

        entry.peptideInfoSet = new HashSet<>(peptideInfoMap.values());
        entry.scanName = this.scanName;

        return entry;
    }
    public class Entry {

        public Set<PeptideInfo> peptideInfoSet;
        public double precursorMass = PreSearch.this.precursorMass;
        public int precursorCharge = PreSearch.this.precursorCharge;
        public List<Peptide> ptmOnlyList = new ArrayList<>();
        public List<Peptide> ptmFreeList = new ArrayList<>();

        public String scanName;
        Entry() {
        }
    }
}
