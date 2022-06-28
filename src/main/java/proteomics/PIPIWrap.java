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
import proteomics.Search.Search;
import proteomics.Segment.InferSegment;
import ProteomicsLibrary.PrepareSpectrum;
import ProteomicsLibrary.MassTool;
import ProteomicsLibrary.Types.*;
import proteomics.Spectrum.PreSpectra;
import proteomics.Types.*;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;

import java.sql.*;
import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.locks.ReentrantLock;

public class PIPIWrap implements Callable<PIPIWrap.Entry> {
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
    private final Binomial binomial;
    private final int scanNum;
    private final int precursorScanNo;
    private final Set<Peptide> ptmOnlyList;
    private final Set<Peptide> ptmFreeList;



    public PIPIWrap(int scanNum, BuildIndex buildIndex, MassTool massTool, double ms1Tolerance, double leftInverseMs1Tolerance, double rightInverseMs1Tolerance, int ms1ToleranceUnit, double ms2Tolerance
            , double minPtmMass, double maxPtmMass, int localMaxMs2Charge, JMzReader spectraParser, double minClear, double maxClear, ReentrantLock lock, String scanId, int precursorCharge
            , double precursorMass, InferPTM inferPTM, PrepareSpectrum preSpectrum, String sqlPath, Binomial binomial, int precursorScanNo,  Set<Peptide> ptmOnlyList,  Set<Peptide> ptmFreeList) {
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
        this.binomial = binomial;
        peptide0Map = buildIndex.getPeptide0Map();
        this.scanNum = scanNum;
        this.precursorScanNo = precursorScanNo;
        this.ptmOnlyList = ptmOnlyList;
        this.ptmFreeList = ptmFreeList;
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
        TreeMap<Double, Double> plMap = preSpectrum.
                preSpectrumTopNStyle(rawPLMap, precursorMass, precursorCharge, minClear, maxClear, PreSpectra.topN);

        if (plMap.isEmpty()) {
            return null;
        }

        // Coding
        InferSegment inferSegment = buildIndex.getInferSegment();
        List<ThreeExpAA> expAaLists = inferSegment.inferSegmentLocationFromSpectrum(precursorMass, plMap, scanNum);
//        for(ThreeExpAA tag : expAaLists) {
//            System.out.println(tag.getPtmFreeAAString());
//        }
        if (true) {
            SparseVector scanCode = inferSegment.generateSegmentIntensityVector(expAaLists);

            // Begin search.
//            Search search = new Search(scanNum, buildIndex, precursorMass, scanCode, massTool, ms1Tolerance, leftInverseMs1Tolerance, rightInverseMs1Tolerance, ms1ToleranceUnit, minPtmMass, maxPtmMass, localMaxMs2Charge);
            // prepare the spectrum

//            String pepHighestSimiScore = "";
//            double highestScore = -1;
//            for (Peptide peptide : search.getPTMOnlyResult()) {
//                if (peptide.getNormalizedCrossCorr() >= highestScore) {
//                    pepHighestSimiScore = peptide.getPTMFreePeptide();
//                    highestScore = peptide.getNormalizedCrossCorr();
//                }
//            }
//
//            for (Peptide peptide : search.getPTMFreeResult()) {
//                if (peptide.getNormalizedCrossCorr() >= highestScore) {
//                    pepHighestSimiScore = peptide.getPTMFreePeptide();
//                    highestScore = peptide.getNormalizedCrossCorr();
//                }
//            }


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
//            GRBEnv env = new GRBEnv(true);
//            env.set(GRB.IntParam.OutputFlag,0);
//            env.start();
            int indexOfPtmContainCandidates = 1;
            for (Peptide peptide : ptmOnlyList) {
                Peptide0 peptide0 = peptide0Map.get(peptide.getPTMFreePeptide());

                PeptidePTMPattern peptidePTMPattern = inferPTM.findPTM(scanNum, expProcessedPL, plMap, precursorMass, peptide.getPTMFreePeptide(), peptide.isDecoy(), peptide.getNormalizedCrossCorr(), peptide0.leftFlank, peptide0.rightFlank, peptide.getGlobalRank(), precursorCharge, localMaxMs2Charge, localMS1ToleranceL, localMS1ToleranceR, expAaLists);

                if (!peptidePTMPattern.getPeptideTreeSet().isEmpty()) {
                    for (Peptide tempPeptide : peptidePTMPattern.getPeptideTreeSet()) {
                        tempPeptide.bestPep = peptidePTMPattern.bestPep;
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
//            env.dispose();
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

            if (!peptideSet.isEmpty()) {
                Peptide[] peptideArray = peptideSet.toArray(new Peptide[0]);
                Peptide topPeptide = peptideArray[0];
//                if (topPeptide.shouldPTM && !topPeptide.hasVarPTM()){
//
//                    if (topPeptide.bestPep != null && topPeptide.bestPep.getScore() > topPeptide.getScore()) {
//                        topPeptide = topPeptide.bestPep;
//                    }
//                }
                TreeSet<Peptide> ptmPatterns = null;
                if (topPeptide.hasVarPTM()) {
                    ptmPatterns = modSequences.get(topPeptide.getPTMFreePeptide());
                }
                new CalSubscores(topPeptide, ms2Tolerance, plMap, precursorCharge, ptmPatterns, binomial);

                Connection sqlConnection = DriverManager.getConnection(sqlPath);
                Statement sqlStatement = sqlConnection.createStatement();
                ResultSet sqlResultSet = sqlStatement.executeQuery(String.format(Locale.US, "SELECT scanNum, scanId, precursorCharge, precursorMass, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient, isDecoy, score FROM spectraTable WHERE scanId='%s'", scanId));
                if (sqlResultSet.next()) {
                    int scanNum = sqlResultSet.getInt("scanNum");
                    String scanId = sqlResultSet.getString("scanId");
                    int precursorCharge = sqlResultSet.getInt("precursorCharge");
                    double precursorMass = sqlResultSet.getDouble("precursorMass");
                    String mgfTitle = sqlResultSet.getString("mgfTitle");
                    int isotopeCorrectionNum = sqlResultSet.getInt("isotopeCorrectionNum");
                    double ms1PearsonCorrelationCoefficient = sqlResultSet.getDouble("ms1PearsonCorrelationCoefficient");

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
//                        deltaCn = (peptideArray[0].getScore() - peptideArray[1].getScore()) / peptideArray[0].getScore();
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
                        Peptide0 pep0 = peptide0Map.get(peptide.getPTMFreePeptide());
                        pepSetString += peptide.getPTMFreePeptide() + "," + peptide.getScore() + "," + peptide.isDecoy() + "," + peptide.hasVarPTM() + "," + String.join("_", pep0.proteins) +",";
                    }
                    List<PepWithScore> candidatesList = new ArrayList<>();

//                    Collections.sort(candidatesList);
                    String candidatesString = "";
                    for (PepWithScore pS : candidatesList) {
                        candidatesString += pS.pepSeq + "," + pS.score + "," + pS.isDecoy + "," + pS.hasPTM + "," + pS.proteins + ",";
                    }
                    if (scanNum == 2327) {
                        System.out.println("lsz");
                    }
                    boolean shouldPtm = Math.abs(precursorMass-massTool.calResidueMass(topPeptide.getPTMFreePeptide()) - massTool.H2O) > 0.01;
                    boolean hasPTM = topPeptide.hasVarPTM();
                    int ptmNum = 0;
                    boolean isSettled = true;
                    double totalPtmMass = 0;
                    if (hasPTM) {
                        ptmNum = topPeptide.getVarPTMs().size();
                        for (double mass : topPeptide.getVarPTMs().values()){
                            totalPtmMass += mass;
                        }
                    }
                    isSettled = Math.abs(totalPtmMass-(precursorMass-massTool.calResidueMass(topPeptide.getPTMFreePeptide()) - massTool.H2O)) <= 0.01;

                    Entry entry = new Entry(scanNum, scanId, shouldPtm ? 1 : 0, hasPTM ? 1 : 0, ptmNum, isSettled ? 1 : 0, precursorCharge, precursorMass, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient, buildIndex.getLabelling(), topPeptide.getPtmContainingSeq(buildIndex.returnFixModMap()), topPeptide.getTheoMass(), topPeptide.isDecoy() ? 1 : 0, topPeptide.getGlobalRank(), topPeptide.getNormalizedCrossCorr(), topPeptide.getScore(), deltaLCn, deltaCn, topPeptide.getMatchedPeakNum(), topPeptide.getIonFrac(), topPeptide.getMatchedHighestIntensityFrac(), topPeptide.getExplainedAaFrac(), otherPtmPatterns, topPeptide.getaScore(), candidatesString, pepSetString.substring(0,pepSetString.length()-1), whereIsTopCand);

                    sqlResultSet.close();
                    sqlStatement.close();
                    return entry;
                } else {
                    throw new NullPointerException(String.format(Locale.US, "There is no record %s in the spectraTable.", scanId));
                }
            } else {
                return null;
            }
        } else {
            return null;
        }
    }


    public class Entry {

        final int scanNum;
        final String scanId;
        final int precursorCharge;
        final double precursorMass;
        final String mgfTitle;
        final int isotopeCorrectionNum;
        final double ms1PearsonCorrelationCoefficient;
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

        Entry(int scanNum, String scanId, int shouldPtm, int hasPTM, int ptmNum, int isSetteld, int precursorCharge, double precursorMass, String mgfTitle, int isotopeCorrectionNum, double ms1PearsonCorrelationCoefficient, String labelling, String peptide, double theoMass, int isDecoy, int globalRank, double normalizedCorrelationCoefficient, double score, double deltaLCn, double deltaCn, int matchedPeakNum, double ionFrac, double matchedHighestIntensityFrac, double explainedAaFrac, String otherPtmPatterns, String aScore, String candidates, String peptideSet, int whereIsTopCand) {
            this.scanNum = scanNum;
            this.scanId = scanId;
            this.shouldPtm = shouldPtm;
            this.hasPTM = hasPTM;
            this.ptmNum = ptmNum;
            this.isSettled = isSetteld;
            this.precursorCharge = precursorCharge;
            this.precursorMass = precursorMass;
            this.mgfTitle = mgfTitle;
            this.isotopeCorrectionNum = isotopeCorrectionNum;
            this.ms1PearsonCorrelationCoefficient = ms1PearsonCorrelationCoefficient;
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
