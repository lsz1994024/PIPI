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

import ProteomicsLibrary.Score;
//import gurobi.GRB;
//import gurobi.GRBEnv;
//import gurobi.GRBModel;
import proteomics.Index.BuildIndex;
import proteomics.PTM.InferPTM;
import ProteomicsLibrary.Binomial;
import proteomics.Search.CalSubscores;
import ProteomicsLibrary.SpecProcessor;
import ProteomicsLibrary.MassTool;
import ProteomicsLibrary.Types.*;
import proteomics.Spectrum.DatasetReader;
import proteomics.Types.*;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;

import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.locks.ReentrantLock;

import static proteomics.PIPI.lszDebugScanNum;

public class PtmSearch implements Callable<PtmSearch.Entry> {
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
    private final Map<String, PepInfo> peptide0Map;
    private final JMzReader spectraParser;
    private final double minClear;
    private final double maxClear;
    private final ReentrantLock lock;
    private final String scanName;
    private final int precursorCharge;
    private final double precursorMass;
    private final InferPTM inferPTM;
    private final SpecProcessor preSpectrum;
    private final Binomial binomial;
    private final int scanNum;
    private final Set<Peptide> ptmOnlyList;
    private final Set<Peptide> ptmFreeList;
    private Map<String, PeptideInfo> peptideInfoMap;

    private boolean isHomo(Peptide p1, Peptide p2) {

        Set<String> temp1 = peptideInfoMap.get(p1.getPTMFreePeptide()).protIdSet;
        Set<String> temp2 = peptideInfoMap.get(p2.getPTMFreePeptide()).protIdSet;

        Set<String> set = new HashSet<>(temp1);
        set.retainAll(temp2);
        if (set.isEmpty()) return false;
//        peptide0Map.get(p1.getPTMFreePeptide()).code.
        SparseBooleanVector sbv1 = buildIndex.inferSegment.generateSegmentBooleanVector(p1.getPTMFreePeptide());
        SparseBooleanVector sbv2 = buildIndex.inferSegment.generateSegmentBooleanVector(p2.getPTMFreePeptide());
        return sbv1.dot(sbv2) > 0.3*Math.min(p1.getPTMFreePeptide().length(), p2.getPTMFreePeptide().length());
    }

    public PtmSearch(int scanNum, BuildIndex buildIndex, MassTool massTool, double ms1Tolerance, double leftInverseMs1Tolerance, double rightInverseMs1Tolerance, int ms1ToleranceUnit, double ms2Tolerance
            , double minPtmMass, double maxPtmMass, int localMaxMs2Charge, JMzReader spectraParser, double minClear, double maxClear, ReentrantLock lock, String scanName, int precursorCharge
            , double precursorMass, InferPTM inferPTM, SpecProcessor preSpectrum, String sqlPath, Binomial binomial, int precursorScanNo, Set<Peptide> ptmOnlyList, Set<Peptide> ptmFreeList, Map<String, PeptideInfo> peptideInfoMap)  {
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
        this.scanName = scanName;
        this.precursorCharge = precursorCharge;
        this.precursorMass = precursorMass;
        this.inferPTM = inferPTM;
        this.preSpectrum = preSpectrum;
        this.binomial = binomial;
        peptide0Map = buildIndex.getPepInfoMap();
        this.scanNum = scanNum;
        this.ptmOnlyList = ptmOnlyList;
        this.ptmFreeList = ptmFreeList;
        this.peptideInfoMap = peptideInfoMap;
    }

    @Override
    public Entry call() throws Exception {
        Map<Double, Double> rawPLMap;
        try {
            lock.lock();
            // Reading peak list.
            rawPLMap = spectraParser.getSpectrumById(scanName.split("\\.")[1]).getPeakList();
        } finally {
            lock.unlock();
        }

        // preprocess peak list
        TreeMap<Double, Double> plMap = preSpectrum.preSpectrumTopNStyleWithChargeLimit(rawPLMap, precursorMass, precursorCharge, minClear, maxClear, DatasetReader.topN);

        if (plMap.isEmpty()) {
            return null;
        }

        if (lszDebugScanNum.contains(scanNum)) {
            int a = 1;
        }
        // Coding
        if (true) {
            SparseVector expProcessedPL;
            if (PIPI.useXcorr) {
                expProcessedPL = preSpectrum.prepareXCorr(plMap, false);
            } else {
                expProcessedPL = preSpectrum.digitizePL(plMap);
            }

            double localMS1ToleranceL = -1 * ms1Tolerance;
            double localMS1ToleranceR = ms1Tolerance;
            if (ms1ToleranceUnit == 1) {
                localMS1ToleranceL = (precursorMass * leftInverseMs1Tolerance) - precursorMass;
                localMS1ToleranceR = (precursorMass * rightInverseMs1Tolerance) - precursorMass;
            }

            // infer PTM using the new approach// my first modification
            TreeSet<Peptide> peptideSet = new TreeSet<>(Collections.reverseOrder());
            Map<String, TreeSet<Peptide>> modSequences = new TreeMap<>();
            int whereIsTopCand = 0; // 0 for still top, -1 for no PTM pattern, -2 for PTM free but score < 0, other number is the final ranking
            PeptideInfo pepInfo;
            for (Peptide peptide : ptmOnlyList) {
                pepInfo = peptideInfoMap.get(peptide.getPTMFreePeptide());

                PeptidePTMPattern peptidePTMPattern = inferPTM.findPTM(scanNum, expProcessedPL, plMap, precursorMass, peptide.getPTMFreePeptide(), peptide.isDecoy(), peptide.getNormalizedCrossCorr()
                        , pepInfo.leftFlank, pepInfo.rightFlank, peptide.getGlobalRank(), precursorCharge, localMaxMs2Charge, localMS1ToleranceL, localMS1ToleranceR);

                if (!peptidePTMPattern.getPeptideTreeSet().isEmpty()) {
                    for (Peptide tempPeptide : peptidePTMPattern.getPeptideTreeSet()) {
                        tempPeptide.bestPep = peptidePTMPattern.bestPep;

//                        if (Math.abs(massTool.calResidueMass(tempPeptide.getVarPtmContainingSeqNow())+ massTool.H2O - precursorMass) > 1) {
//                            continue;// skip the unsettle peptides
//                        }
                        if (tempPeptide.getScore() > 0) {
                            if (peptideSet.size() < candisNum) {
                                peptideSet.add(tempPeptide);
                            } else if (tempPeptide.getScore() > peptideSet.last().getScore()) {
                                peptideSet.pollLast();
                                peptideSet.add(tempPeptide);
                            }
                        }
                    }
                    // record scores with different PTM patterns for calculating PTM delta score.
                    modSequences.put(peptidePTMPattern.ptmFreePeptide, peptidePTMPattern.getPeptideTreeSet());
                }
            }
            // Calculate Score for PTM free peptide  for PTM free score is calXcorr score, for PTM only score is PTM score
            for (Peptide peptide : ptmFreeList) {
                double score = massTool.buildVectorAndCalXCorr(peptide.getIonMatrix(), precursorCharge, expProcessedPL);
                if (score > 0) {
                    peptide.setScore(score + peptide.getNormalizedCrossCorr());
                    double[][] temp =peptide.getIonMatrix();
                    if (temp.length == 1) {
                        int a = 1;
                    }
                    peptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, temp, ms2Tolerance));
                    if (peptideSet.size() < candisNum) {
                        peptideSet.add(peptide);
                    } else if (peptide.getScore() > peptideSet.last().getScore()) {
                        peptideSet.pollLast();
                        peptideSet.add(peptide);
                    }
                }
            }
            // Calculate Score for PTM free peptide  for PTM free score is calXcorr score, for PTM only score is PTM score

            if (peptideSet.isEmpty()) {
                return null;
            }

            List<Peptide> pepList = new ArrayList<>(peptideSet);

            for (int j = pepList.size()-1; j >= 1; j--) {
                for (int i = 0; i < j; i++) {
                    if (isHomo(pepList.get(i), pepList.get(j))) {    // to clean homo pep candidates from the list
                        boolean isIGood = Math.abs(massTool.calResidueMass(pepList.get(i).getVarPtmContainingSeqNow()) + massTool.H2O - precursorMass) < 0.5;
                        boolean isJGood = Math.abs(massTool.calResidueMass(pepList.get(j).getVarPtmContainingSeqNow()) + massTool.H2O - precursorMass) < 0.5;
                        if (!isIGood && isJGood) {
                            pepList.remove(i);
                        } else {
                            pepList.remove(j);
                        }
                        break;
                    }
                }
            }


            Peptide[] peptideArray = pepList.toArray(new Peptide[0]);

            Peptide topPep = peptideArray[0];
            Entry entry;
            List<ExtraEntry> extraEntryList = new ArrayList<>();
            TreeSet<Peptide> ptmPatterns = null;
            if (topPep.hasVarPTM()) {
                ptmPatterns = modSequences.get(topPep.getPTMFreePeptide());
            }
            new CalSubscores(topPep, ms2Tolerance, plMap, precursorCharge, ptmPatterns, binomial);


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
            if (ptmPatterns != null) {
                List<String> tempList = new LinkedList<>();
                Iterator<Peptide> ptmPatternsIterator = ptmPatterns.iterator();
                ptmPatternsIterator.next();
                while (ptmPatternsIterator.hasNext()) {
                    Peptide temp = ptmPatternsIterator.next();
                    tempList.add(String.format(Locale.US, "%s-%.4f", temp.getPtmContainingSeq(buildIndex.returnFixModMap()), temp.getScore())); // Using 4 decimal here because it is write the the result file for checking. It is not used in scoring or other purpose.
                }
                otherPtmPatterns = String.join(";", tempList);
            }
            String pepSetString = "";
            for (Peptide peptide : peptideArray){
                pepInfo = peptideInfoMap.get(peptide.getPTMFreePeptide());
                pepSetString += peptide.getVarPtmContainingSeqNow() + "," + peptide.getScore() + "," + String.join("_", pepInfo.protIdSet) +",";
            }    // pepSetString saves ptmContaining seq

            boolean shouldPtm = Math.abs(precursorMass-massTool.calResidueMass(topPep.getPTMFreePeptide()) - massTool.H2O) > 0.01;
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
            isSettled = Math.abs(totalPtmMass-(precursorMass-massTool.calResidueMass(topPep.getPTMFreePeptide()) - massTool.H2O)) <= 0.01;
            entry = new Entry(
                    scanNum, scanName, shouldPtm ? 1 : 0, hasPTM ? 1 : 0, ptmNum, isSettled ? 1 : 0
                    , precursorCharge, precursorMass, buildIndex.getLabelling(), topPep.getPtmContainingSeq(buildIndex.returnFixModMap())
                    , topPep.getTheoMass(), topPep.isDecoy() ? 1 : 0, topPep.getGlobalRank()
                    , topPep.getNormalizedCrossCorr(), topPep.getScore(), deltaLCn, deltaCn
                    , topPep.getMatchedPeakNum(), topPep.getIonFrac(), topPep.getMatchedHighestIntensityFrac()
                    , topPep.getExplainedAaFrac(), otherPtmPatterns, topPep.getaScore(), ""
                    , pepSetString.substring(0, pepSetString.length()-1), whereIsTopCand, extraEntryList
                    );

            return entry;
        } else {
            return null;
        }
    }
    public class ExtraEntry {

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
        final int whereIsTopCand;
        final int hasPTM;
        final int ptmNum;
        final int isSettled;
        final int shouldPtm;

        ExtraEntry(int scanNum, String scanName, int shouldPtm, int hasPTM, int ptmNum, int isSetteld, int precursorCharge, double precursorMass,String labelling, String peptide, double theoMass, int isDecoy, int globalRank, double normalizedCorrelationCoefficient, double score, double deltaLCn, double deltaCn, int matchedPeakNum, double ionFrac, double matchedHighestIntensityFrac, double explainedAaFrac, String otherPtmPatterns, String aScore, String candidates, String peptideSet, int whereIsTopCand) {
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
            this.whereIsTopCand = whereIsTopCand;
        }
    }
    public class Entry {

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
        final int whereIsTopCand;
        final int hasPTM;
        final int ptmNum;
        final int isSettled;
        final int shouldPtm;
        final List<ExtraEntry> extraEntryList;
        Entry(int scanNum, String scanName, int shouldPtm, int hasPTM, int ptmNum, int isSetteld, int precursorCharge, double precursorMass
                ,String labelling, String peptide, double theoMass, int isDecoy, int globalRank, double normalizedCorrelationCoefficient
                , double score, double deltaLCn, double deltaCn, int matchedPeakNum, double ionFrac, double matchedHighestIntensityFrac
                , double explainedAaFrac, String otherPtmPatterns, String aScore, String candidates, String peptideSet, int whereIsTopCand, List<ExtraEntry> extraEntryList) {
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
            this.whereIsTopCand = whereIsTopCand;
            this.extraEntryList = extraEntryList;
        }
    }
}
