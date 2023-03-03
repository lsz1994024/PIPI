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

import static proteomics.PIPI.isPtmSimuTest;
import static proteomics.PIPI.lszDebugScanNum;
import static proteomics.PTM.InferPTM.*;

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
    private final PriorityQueue<Peptide> peptidePriorityQueue;
    private Map<String, PeptideInfo> peptideInfoMap;

    private boolean isHomo(Peptide p1, Peptide p2) {

        Set<String> temp1 = peptideInfoMap.get(p1.getFreeSeq()).protIdSet;
        Set<String> temp2 = peptideInfoMap.get(p2.getFreeSeq()).protIdSet;

        Set<String> set = new HashSet<>(temp1);
        set.retainAll(temp2);
//        if (set.isEmpty()) return false;
        SparseBooleanVector sbv1 = buildIndex.inferSegment.generateSegmentBooleanVector(p1.getFreeSeq());
        SparseBooleanVector sbv2 = buildIndex.inferSegment.generateSegmentBooleanVector(p2.getFreeSeq());
        return sbv1.dot(sbv2) > 0.3*Math.min(p1.getFreeSeq().length(), p2.getFreeSeq().length());
    }

    public PtmSearch(int scanNum, BuildIndex buildIndex, MassTool massTool, double ms1Tolerance, double leftInverseMs1Tolerance, double rightInverseMs1Tolerance, int ms1ToleranceUnit, double ms2Tolerance
            , double minPtmMass, double maxPtmMass, int localMaxMs2Charge, JMzReader spectraParser, double minClear, double maxClear, ReentrantLock lock, String scanName, int precursorCharge
            , double precursorMass, InferPTM inferPTM, SpecProcessor preSpectrum, Binomial binomial, PriorityQueue<Peptide> peptidePriorityQueue, Map<String, PeptideInfo> peptideInfoMap)  {
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
        this.peptidePriorityQueue = peptidePriorityQueue;
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
        TreeMap<Double, Double> plMap = preSpectrum.preSpectrumTopNStyleWithChargeLimit(rawPLMap, precursorMass, precursorCharge, minClear, maxClear, DatasetReader.topN,ms2Tolerance);

        if (plMap.isEmpty()) {
            return null;
        }

        if (lszDebugScanNum.contains(scanNum)) {
            int a = 1;
        }
        double ms1TolAbs = Double.parseDouble(InferPTM.df.format(precursorMass*ms1Tolerance/1000000));

        // Coding
        SparseVector expProcessedPL;
        if (PIPI.useXcorr) {
            expProcessedPL = preSpectrum.prepareXCorr(plMap, false);
        } else {
            expProcessedPL = preSpectrum.digitizePL(plMap);
        }


        // infer PTM using the new approach// my first modification
        TreeSet<Peptide> peptideSet = new TreeSet<>(Collections.reverseOrder());
        Map<String, TreeSet<Peptide>> modSequences = new TreeMap<>();
        int whereIsTopCand = 0; // 0 for still top, -1 for no PTM pattern, -2 for PTM free but score < 0, other number is the final ranking

        for (Peptide peptide : peptidePriorityQueue) {
            PeptideInfo peptideInfo = peptideInfoMap.get(peptide.getFreeSeq());
            String pepSeq = peptide.getFreeSeq();
            Map<Integer, VarPtm> refVarPtmMap = peptide.posVarPtmResMap;

            byte isProtNorC_Term = NON_TERM_PROT;
            if (peptideInfo.leftFlank == '-' || peptideInfo.leftFlank == 'M'){
                isProtNorC_Term = N_TERM_PROT;
            } else if (peptideInfo.rightFlank == '-') {
                isProtNorC_Term = C_TERM_PROT;
            }

            double deltaMass = precursorMass - (massTool.calResidueMass(peptide.getFreeSeq()) + massTool.H2O);
            if (!peptide.posVarPtmResMap.isEmpty()) { // ptm containing todo this criteria might be wrong
                Set<Integer> fixModIdxes = inferPTM.getFixModIdxes(pepSeq);
                Map<Integer, VarPtm[]> posVarModArrayMap = inferPTM.getIdxVarModMapNewComple(pepSeq, fixModIdxes, isProtNorC_Term); //todo no need to generate var mod list for aa again and again, make it stored.
                ModPepPool modPepsSettled = inferPTM.settlePtmOnSideComple(scanNum, expProcessedPL, plMap, precursorMass, pepSeq, false,
                        posVarModArrayMap, deltaMass, precursorCharge,ms1TolAbs, refVarPtmMap);

                if (modPepsSettled.peptideTreeSet.isEmpty()) {
                    continue; // how could it be return!!!
                }
                for (Peptide tempPeptide : modPepsSettled.getPeptideTreeSet()) {
                    tempPeptide.bestPep = modPepsSettled.bestPep;
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
                modSequences.put(modPepsSettled.freeSeq, modPepsSettled.getPeptideTreeSet());
            } else {  // the comple pep does not have PTM
                double score = massTool.buildVectorAndCalXCorr(peptide.getIonMatrixNow(), 1, expProcessedPL);
                if (score > 0) {
                    peptide.setScore(score + 0);
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
        }
        // Calculate Score for PTM free peptide  for PTM free score is calXcorr score, for PTM only score is PTM score
//        for (Peptide peptide : ptmFreeList) {
//            ExpTag tag = peptide.finderTag;
//
//            if ( tag != null && tag.getPtmAaString().length() != tag.getFreeAaString().length()) { // means the tag is with label var mod, need to put this mod on the final pep
//                PosMassMap newPtmPtn = new PosMassMap(peptide.getFreeSeq().length());
//                int idOfAa = -1;
//                for (char letter : tag.getPtmAaString().toCharArray()) {
//                    if (Character.isUpperCase(letter)) {
//                        idOfAa += 1;
//                    } else {
//                        newPtmPtn.put(peptide.tagPosInPep+idOfAa, massTool.labelVarPtmMap.get(letter).mass);
//                    }
//                }
//                peptide.setVarPTM(newPtmPtn);
//                modSequences.put(peptide.getFreeSeq(), new TreeSet<>(Arrays.asList(peptide))); // just a dummy
//            }
//            double score = massTool.buildVectorAndCalXCorr(peptide.getIonMatrixNow(), 1, expProcessedPL);
//            if (score > 0) {
//                peptide.setScore(score + 0);
//                double[][] temp =peptide.getIonMatrix();
//                if (temp.length == 1) {
//                    int a = 1;
//                }
//                peptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, temp, ms2Tolerance));
//                if (peptideSet.size() < candisNum) {
//                    peptideSet.add(peptide);
//                } else if (peptide.getScore() > peptideSet.last().getScore()) {
//                    peptideSet.pollLast();
//                    peptideSet.add(peptide);
//                }
//            }
//        }
        // Calculate Score for PTM free peptide  for PTM free score is calXcorr score, for PTM only score is PTM score

        if (peptideSet.isEmpty()) {
            return null;
        }

        List<Peptide> pepList = new ArrayList<>(peptideSet);
        Set<Integer> pepIdsToRemove = new HashSet<>();
        for (int j = pepList.size()-1; j >= 1; j--) {
            if (pepIdsToRemove.contains(j)) continue;

            for (int i = 0; i < j; i++) {
                if (pepIdsToRemove.contains(i)) continue;

//                    if (isHomo(pepList.get(i), pepList.get(j))) {    // to clean homo pep candidates from the list
                if (isHomo(pepList.get(i), pepList.get(j)) || pepList.get(j).getScore() > 0.6*pepList.get(i).getScore()) {
                    int iPriority = pepList.get(i).getPriority();
                    int jPriority = pepList.get(j).getPriority();
                    if (iPriority < jPriority) {
                        pepIdsToRemove.add(i);
                    } else  {
                        pepIdsToRemove.add(j);
                    }
//                        break;
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

        for (Peptide pep : newPepList){
            if (pep.getPriority() < 0 && isPtmSimuTest) {
                System.out.println(scanNum + " wtf ? " + pep.getVarPtmContainingSeqNow());
            }
        }

        Peptide[] peptideArray = newPepList.toArray(new Peptide[0]);

        Peptide topPep = peptideArray[0];
        Entry entry;
        TreeSet<Peptide> ptmPatterns = null;
        if (topPep.hasVarPTM()) {
            ptmPatterns = modSequences.get(topPep.getFreeSeq());
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
            PeptideInfo peptideInfo = peptideInfoMap.get(peptide.getFreeSeq());
            pepSetString += peptide.getVarPtmContainingSeqNow() + "," + peptide.getScore() + "," + String.join("_", peptideInfo.protIdSet) +",";
        }    // pepSetString saves ptmContaining seq

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
        entry = new Entry(
                scanNum, scanName, shouldPtm ? 1 : 0, hasPTM ? 1 : 0, ptmNum, isSettled ? 1 : 0
                , precursorCharge, precursorMass, buildIndex.getLabelling(), topPep.getPtmContainingSeq(buildIndex.returnFixModMap())
                , topPep.getTheoMass(), topPep.isDecoy() ? 1 : 0, topPep.getGlobalRank()
                , topPep.getTagVecScore(), topPep.getScore(), deltaLCn, deltaCn
                , topPep.getMatchedPeakNum(), topPep.getIonFrac(), topPep.getMatchedHighestIntensityFrac()
                , topPep.getExplainedAaFrac(), otherPtmPatterns, topPep.getaScore(), ""
                , pepSetString.substring(0, pepSetString.length()-1), whereIsTopCand
                );

        if (lszDebugScanNum.contains(scanNum)){
            int a = 1;
        }
        return entry;
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
        Entry(int scanNum, String scanName, int shouldPtm, int hasPTM, int ptmNum, int isSetteld, int precursorCharge, double precursorMass
                ,String labelling, String peptide, double theoMass, int isDecoy, int globalRank, double normalizedCorrelationCoefficient
                , double score, double deltaLCn, double deltaCn, int matchedPeakNum, double ionFrac, double matchedHighestIntensityFrac
                , double explainedAaFrac, String otherPtmPatterns, String aScore, String candidates, String peptideSet, int whereIsTopCand ) {
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
}
