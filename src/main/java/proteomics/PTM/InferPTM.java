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

package proteomics.PTM;

import ProteomicsLibrary.*;
import ProteomicsLibrary.MassTool;
import ProteomicsLibrary.Types.*;
//import proteomics.OutputPeff;
import gurobi.*;
import org.apache.commons.math3.util.Pair;
import proteomics.Types.*;

import java.io.*;
import java.text.DecimalFormat;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import org.dom4j.Document;
import org.dom4j.Element;
import org.dom4j.io.SAXReader;

import static proteomics.PIPI.lszDebugScanNum;


public class InferPTM {

    private static final Pattern pattern = Pattern.compile("([0-9A-Za-z]+)(\\(([0-9\\-]+)\\))?");
    private static final double ptmMassTolerance = 0.1;
    public static final byte N_PART = 0;
    public static final byte C_PART = 1;
    public static final byte N_TERM_PROT = -1;
    public static final byte NON_TERM_PROT = 0;
    public static final byte C_TERM_PROT = 1;
    public final static DecimalFormat df3 = new DecimalFormat("0.000");
    public final static DecimalFormat df2 = new DecimalFormat("0.00");
    private final MassTool massTool;
    private final Map<String, Double> elementTable;
    private final Map<Character, Double> massTable;
    private final Map<Character, Double> fixModMap;
    private Set<VarPtm> varPtmSet = new HashSet<>();
    private final double minPtmMass;
    private final double maxPtmMass;
    public double minUserPtmMass = 0;
    public double maxUserPtmMass = 530;
    private final double ms2Tolerance;
    private Map<Character, List<VarPtm>> finalPtmMap = new HashMap<>();
    private final Set<Character> aaCharSet = new HashSet<>(Arrays.asList('A','C','D','E','F','G','H','I','K','L','M','N','O','P','Q','R','S','T','U','V','W','Y'));

    private Set<Character> aaWithFixModSet = new HashSet<>();
    public InferPTM(MassTool massTool, Map<Character, Double> fixModMap, Map<String, String> parameterMap) throws Exception{
        this.massTool = massTool;
        elementTable = massTool.getElementTable();
        massTable = massTool.getMassTable();
        this.fixModMap = fixModMap;
        for (Character c : fixModMap.keySet()){
            if (Math.abs(fixModMap.get(c)) > 0.02) {
                aaWithFixModSet.add(c);
            }
        }
//        this.minPtmMass = Math.min(Double.valueOf(parameterMap.get("min_ptm_mass")), -600);// todo, this is correct way to handle user specified big mass PTM
        this.minPtmMass = Double.valueOf(parameterMap.get("min_ptm_mass"));

        this.maxPtmMass = Double.valueOf(parameterMap.get("max_ptm_mass"));
        this.ms2Tolerance = Double.valueOf(parameterMap.get("ms2_tolerance"));

        char[] aaArray = new char[]{'G', 'A', 'S', 'P', 'V', 'T', 'C', 'L', 'N', 'D', 'Q', 'K', 'E', 'M', 'H', 'F', 'R', 'Y', 'W'};
        int n_varPtm = 0;
        for (String k : parameterMap.keySet()) {
            if (!k.startsWith("mod")) continue;
            String v = parameterMap.get(k);
            if (v.startsWith("0.0")) break;
            n_varPtm++ ;
        }
        for (String k : parameterMap.keySet()) {
            if (!k.startsWith("mod")) continue;

            String v = parameterMap.get(k);
            if (v.startsWith("0.0")) break;

            //15.994915,M,0,Oxidation
            String[] modStr = v.split(",");
            double modMass = Double.valueOf(modStr[0]);
            char modSite = modStr[1].charAt(0);
            int modPosition = Integer.valueOf(modStr[2]);
            int priority = 1;
            if (modSite == 'M' && modStr[3].contentEquals("Oxidation") && n_varPtm != 1) {
                priority = 0; // oxidation on M is a common phonomenon but not an enriched modification
            }
            if (modPosition == 4) {//  position anywhere, highest prority
                //must happen in one aa, X (any aa) not allowed
                if (Math.abs(fixModMap.get(modSite)) < 0.1) {
                    varPtmSet.add(new VarPtm(modMass, modSite, modPosition, modStr[3], "ByUser", priority)); // var mods from the parameter file have the highest priority, those PTM can exist in peptide terminal.
                }
            } else {// position N C term, middle  prority // when find tags we can't differ pepN from protN
                if (modSite == 'X') { //on any aa
                    for (char oriAa : aaArray){
                        if (Math.abs(fixModMap.get(modSite)) < 0.1) {
                            varPtmSet.add(new VarPtm(modMass, oriAa, modPosition, modStr[3], "ByUser", priority)); // var mods from the parameter file have the highest priority, those PTM can exist in peptide terminal.
                        }
                    }
                } else { //on single aa
                    if (Math.abs(fixModMap.get(modSite)) < 0.1) {
                        varPtmSet.add(new VarPtm(modMass, modSite, modPosition, modStr[3], "ByUser", priority)); // var mods from the parameter file have the highest priority, those PTM can exist in peptide terminal.
                    }
                }
            }
        }

        // Multimap<Character, ModEntry>
        // Reading Mod table...
        readModFromUnimod();

        // update ptm table with the high priority mods in the parameter file.
        for (VarPtm varPtm : varPtmSet) {
            if (finalPtmMap.containsKey(varPtm.site)) {
                finalPtmMap.get(varPtm.site).add(varPtm);
            } else {
                List<VarPtm> tempList = new LinkedList<>();
                tempList.add(varPtm);
                finalPtmMap.put(varPtm.site, tempList);
            }
        }
    }

    public ModPepPool findPtmNew1(int scanNum, SparseVector expProcessedPL, TreeMap<Double, Double> plMap, double precursorMass, Peptide candiPep, PeptideInfo peptideInfo
            , int precursorCharge, int localMaxMS2Charge, double localMS1ToleranceL, double localMS1ToleranceR, Set<VarPtm> refVarPtmList) {
//        double ptmFreeMass = massTool.calResidueMass(ptmFreePeptide) + massTool.H2O;
        // use the old way of findPtm but blocking tag region in the way as I block region with fixed mod
        String freeSeq = candiPep.getFreeSeq();
        boolean isDecoy = candiPep.isDecoy;
        char leftFlank = peptideInfo.leftFlank;
        char rightFlank = peptideInfo.rightFlank;

//        double nDeltaMass = candiPep.nDeltaMass;
//        double cDeltaMass = candiPep.cDeltaMass;
        int tagPosInPep = candiPep.tagPosInPep;
        int tagLen = 0;

        double score = 0;
        if (candiPep.finderTag != null) {
            score = candiPep.finderTag.getTotalIntensity();
            tagLen = candiPep.finderTag.size();
        }

        ModPepPool modPepPool = new ModPepPool(freeSeq);
        Peptide peptide = new Peptide(freeSeq, isDecoy, massTool);

        double totalDeltaMass = precursorMass - peptide.getTheoMass();
        peptide.absDeltaMass = totalDeltaMass;
        if (Math.abs(totalDeltaMass) < 0.01 ) {
            return modPepPool;
        }
        Set<Integer> fixModIdxes = getFixModIdxes(freeSeq);  // positions that has fixed mod on it. Those postions should not bear var mod then.

//        String ptmFreePeptideOrdinary = freeSeq.replaceAll("I","L");

        int pepPosInProt = 0; // none of the terms is protein term
        if (leftFlank == '-') {
            pepPosInProt = -1; //nterm is protein N term
        } else if (rightFlank == '-') {
            pepPosInProt = 1;  //cterm is protein C term
        }

        Map<Integer, Set<VarPtm>> idxVarModMap = getIdxVarModMapNewOld(freeSeq, fixModIdxes, pepPosInProt, tagPosInPep, tagLen); //todo no need to generate var mod list for aa again and again, make it stored.
        Map<Integer, VarPtm[]> idxVarModArrayMap = new HashMap<>();
        for (int id : idxVarModMap.keySet()){
            VarPtm[] modArray = new VarPtm[idxVarModMap.get(id).size()];
            idxVarModMap.get(id).toArray(modArray);
            Arrays.sort(modArray, Comparator.comparingDouble(VarPtm::getMass));
            idxVarModArrayMap.put(id, modArray);
        }

        Set<Integer> modifiedZone = new HashSet<>(idxVarModArrayMap.keySet());// e.g. nABCDEFc, modifiedZone={1,2,3,4,5,6}

        if (pepPosInProt != 0){
            int a = 1;
        }


        if (lszDebugScanNum.contains(scanNum)) {
            int a = 1;
        }
        if (modifiedZone.size() == 0) {
//            System.out.println(scanNum + " is empty modifiedZone after tag 2");
            return modPepPool; //Some scans are not valid Scans. Will be deleted soon.
        }

        TreeMap<Double, Double> unUsedPlMap = new TreeMap<>(plMap);

        ModPepPool allPtmPattern = new ModPepPool(freeSeq,1);
        ModPepPool allPtmPatternBad = new ModPepPool(freeSeq,1);
        modifiedZone = IntStream.rangeClosed(0, freeSeq.length()-1).boxed().collect(Collectors.toSet());// e.g. nABCDEFc, modifiedZone={1,2,3,4,5,6}

        Peptide cleanPep = new Peptide(freeSeq, isDecoy, massTool);

        int lb = 0;  //lb included
        int rb = freeSeq.length() - 1;//rb included
        Map<Integer, Double> matchedBions = new HashMap<>();
        Map<Integer, Double> matchedYions = new HashMap<>();
        double[][] temp1 = cleanPep.getIonMatrixNow();
        Set<Integer> jRange = IntStream.rangeClosed(0, freeSeq.length()-1).boxed().collect(Collectors.toSet());
        double[][] cleanPepIonMatrix = cleanPep.getIonMatrixNow();
        double cleanScore = massTool.buildVectorAndCalXCorr(cleanPepIonMatrix, 1, expProcessedPL, matchedBions, matchedYions, jRange) ;
        cleanPep.setScore(cleanScore);
        cleanPep.setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, cleanPep.getIonMatrix(), ms2Tolerance));
        cleanPep.matchedBions.putAll(matchedBions);
        cleanPep.matchedYions.putAll(matchedYions);

        List<Integer> matchedBIndex = new ArrayList<>(matchedBions.keySet());
        List<Integer> matchedYIndex = new ArrayList<>(matchedYions.keySet());
        Collections.sort(matchedBIndex);
        Collections.sort(matchedYIndex);
        if (matchedBIndex.size() >= 2) {
            if (matchedBIndex.get(matchedBIndex.size()-1) + (1-matchedBions.get(matchedBIndex.get(matchedBIndex.size()-1))) > matchedBIndex.get(matchedBIndex.size()-2)+4) {
                matchedBions.remove(matchedBIndex.get(matchedBIndex.size()-1));
            }
        }
        if (matchedYIndex.size() >= 2) {
            if (matchedYIndex.get(0) - (1-matchedYions.get(matchedYIndex.get(0))) < matchedYIndex.get(1)-4) {
                matchedYions.remove(matchedYIndex.get(0));
            }
        }

        if (!matchedBions.isEmpty()) {
            lb = Collections.max(matchedBions.keySet()) + 1;
        }
        if (!matchedYions.isEmpty()) {
            rb = Collections.min(matchedYions.keySet()) - 1;
        }

        if (rb - lb <= 0) {
            double bSumIntens = 0;
            for (double intes : matchedBions.values()) bSumIntens += intes;
            double ySumIntens = 0;
            for (double intes : matchedYions.values()) ySumIntens += intes;
            if (bSumIntens > ySumIntens) {
                rb = freeSeq.length() - 1;
            } else {
                lb = 0;
            }
        }
        modifiedZone = IntStream.rangeClosed(lb, rb).boxed().collect(Collectors.toSet());

        allPtmPatternBad.push(cleanPep);
        allPtmPattern.push(cleanPep);
        allPtmPatternBad.bestPep = allPtmPatternBad.getTopPepPtn();
        if (modifiedZone.isEmpty()) {
            return allPtmPatternBad;
        }

        ModPepPool ptmInitialTemp = new ModPepPool(freeSeq, 1);
        Peptide cleanPeptide = new Peptide(freeSeq, isDecoy, massTool);
        cleanPeptide.setVarPTM(new PosMassMap(freeSeq.length()));
        ptmInitialTemp.push(cleanPeptide);
//        long t1 = System.currentTimeMillis();

        // 1st
        ModPepPool ptmN1Good = new ModPepPool(freeSeq, 1);
        ModPepPool ptmN1Bad = new ModPepPool(freeSeq, 1);

        DividedZone z1Z2Res = dividePep(scanNum, modifiedZone, ptmN1Good, ptmN1Bad, ptmInitialTemp, idxVarModArrayMap, totalDeltaMass, freeSeq, isDecoy, expProcessedPL, plMap, precursorCharge, localMaxMS2Charge, unUsedPlMap);
        if (!ptmN1Good.getPeptideTreeSet().isEmpty()) {
            ptmN1Good.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN1Good.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 1*0.1);
            ptmN1Good.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, ptmN1Good.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));

            allPtmPattern.push(ptmN1Good.getTopPepPtn());
        }
        if (!ptmN1Bad.getPeptideTreeSet().isEmpty()) {
            ptmN1Bad.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN1Bad.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 1*0.1);
            ptmN1Bad.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, ptmN1Bad.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
            allPtmPatternBad.push(ptmN1Bad.getTopPepPtn());
        }
        Set<Integer> zone1 = z1Z2Res.settledZone;
        Set<Integer> zone2 = z1Z2Res.toModZone;
        if (ptmN1Bad.peptideTreeSet.isEmpty() || zone2.isEmpty()) {
            if (allPtmPattern.peptideTreeSet.isEmpty()) {
                return allPtmPatternBad;
            }
            allPtmPattern.bestPep = allPtmPatternBad.getTopPepPtn();
            return allPtmPattern;
        }
        Map<Integer, Double> matchedBionsBad1 = new HashMap<>();
        Map<Integer, Double> matchedYionsBad1 = new HashMap<>();
//        ptmN1Bad.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN1Bad.getTopPepPtn().getIonMatrixNow(), 1,  expProcessedPL, matchedBionsBad1, matchedYionsBad1, IntStream.range(0, ptmFreePeptide.length()-2).boxed().collect(Collectors.toSet())));

//        ptmN1Bad.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN1Bad.getTopPepPtn().getIonMatrixNow(), 1,  expProcessedPL) - 1*0.1);
//        ptmN1Bad.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, localMaxMS2Charge, ptmN1Bad.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
//        allPtmPatternBad.push(ptmN1Bad.getTopPepPtn());

        double zone1Mass = ptmN1Bad.getTopPepPtn().getVarPTMs().values().iterator().next();
        double zone2Mass = totalDeltaMass - zone1Mass;

        // 2nd
        ModPepPool ptmN2Good = new ModPepPool(freeSeq, 1);
        ModPepPool ptmN2Bad = new ModPepPool(freeSeq, 1);

        DividedZone z3Z4Res = dividePep(scanNum, zone2, ptmN2Good, ptmN2Bad ,ptmN1Bad, idxVarModArrayMap, zone2Mass, freeSeq, isDecoy, expProcessedPL, plMap, precursorCharge, localMaxMS2Charge, unUsedPlMap);
        if (!ptmN2Good.getPeptideTreeSet().isEmpty()) {
            ptmN2Good.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN2Good.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 2*0.1);
            ptmN2Good.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, ptmN2Good.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));

            allPtmPattern.push(ptmN2Good.getTopPepPtn());
        }
        if (!ptmN2Bad.getPeptideTreeSet().isEmpty()) {
            ptmN2Bad.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN2Bad.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 2*0.1);
            ptmN2Bad.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, ptmN2Bad.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
            allPtmPatternBad.push(ptmN2Bad.getTopPepPtn());
        }
        Set<Integer> zone3 = z3Z4Res.settledZone;
        Set<Integer> zone4 = z3Z4Res.toModZone;
        if (ptmN2Bad.peptideTreeSet.isEmpty() || zone4.isEmpty()) {
            if (allPtmPattern.peptideTreeSet.isEmpty()) {
                return allPtmPatternBad;
            }
            allPtmPattern.bestPep = allPtmPatternBad.getTopPepPtn();
            return allPtmPattern;
        }
//        ptmN2Bad.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN2Bad.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 2*0.1);
//        ptmN2Bad.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, localMaxMS2Charge, ptmN2Bad.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
//        allPtmPatternBad.push(ptmN2Bad.getTopPepPtn());

        double zone3Mass = z3Z4Res.ptmMass;
        double zone4Mass = zone2Mass - zone3Mass;

        // 3rd
        ModPepPool ptmN3Good = new ModPepPool(freeSeq, 1);
        ModPepPool ptmN3Bad = new ModPepPool(freeSeq, 1);

        DividedZone z5Z6Res = dividePep(scanNum, zone4, ptmN3Good, ptmN3Bad, ptmN2Bad, idxVarModArrayMap, zone4Mass, freeSeq, isDecoy, expProcessedPL, plMap, precursorCharge, localMaxMS2Charge, unUsedPlMap);
        if (!ptmN3Good.peptideTreeSet.isEmpty()) {
            ptmN3Good.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN3Good.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 3*0.1);
            ptmN3Good.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, ptmN3Good.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));

            allPtmPattern.push(ptmN3Good.getTopPepPtn());
        }
        if (!ptmN3Bad.getPeptideTreeSet().isEmpty()) {
            ptmN3Bad.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN3Bad.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 3*0.1);
            ptmN3Bad.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, ptmN3Bad.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
            allPtmPatternBad.push(ptmN3Bad.getTopPepPtn());
        }
        Set<Integer> zone5 = z5Z6Res.settledZone;
        Set<Integer> zone6 = z5Z6Res.toModZone;
        if (ptmN3Bad.peptideTreeSet.isEmpty() || zone6.isEmpty()) {
            if (allPtmPattern.peptideTreeSet.isEmpty()) {
                return allPtmPatternBad;
            }
            allPtmPattern.bestPep = allPtmPatternBad.getTopPepPtn();
            return allPtmPattern;
        }
//        ptmN3Bad.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN3Bad.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 3*0.1);
//        ptmN3Bad.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, localMaxMS2Charge, ptmN3Bad.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
//        allPtmPatternBad.push(ptmN3Bad.getTopPepPtn());

        double zone5Mass = z5Z6Res.ptmMass;
        double zone6Mass = zone4Mass - zone5Mass;

        ModPepPool ptmN4Good = new ModPepPool(freeSeq, 1);
        ModPepPool ptmN4Bad = new ModPepPool(freeSeq, 1);
        DividedZone z7Z8Res = dividePep(scanNum, zone6, ptmN4Good, ptmN4Bad ,ptmN3Bad, idxVarModArrayMap, zone6Mass, freeSeq, isDecoy, expProcessedPL, plMap, precursorCharge, localMaxMS2Charge, unUsedPlMap);
        if (!ptmN4Good.getPeptideTreeSet().isEmpty()) {
            ptmN4Good.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN4Good.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 4*0.1);
            ptmN4Good.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, ptmN4Good.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
            allPtmPattern.push(ptmN4Good.getTopPepPtn());
        }
        if (!ptmN4Bad.getPeptideTreeSet().isEmpty()) {
            ptmN4Bad.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN4Bad.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 4*0.1);
            ptmN4Bad.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, ptmN4Bad.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));

            allPtmPatternBad.push(ptmN4Bad.getTopPepPtn());
        }



        long end = System.currentTimeMillis();

//        System.out.println("N1 time, "+ (end - start));
        if (allPtmPattern.peptideTreeSet.isEmpty()) {
            return allPtmPatternBad;
        }
        allPtmPattern.bestPep = allPtmPatternBad.getTopPepPtn();
        return allPtmPattern;
    }

    public List<Map<Double, Integer>> findBestPtmMIP(int scanNum, GRBEnv env, Map<Double, Integer> allMassMaxTimeMap, double totalDeltaMass, String partSeq
            , double cutMass, byte isNorC_Part, SparseVector expProcessedPL, double ms1TolAbs) { // Sometimes, the precursor mass error may affects the digitized spectrum.
        int numTimeout = 0;
        int numOfSols = 0;

//        ModPepPool modPepPool = new ModPepPool(partSeq, 20);
        int numPtmsOnPep = 2;
        long t1 = 0,t2 = 0,t3 = 0,t4 = 0;
        t1=System.currentTimeMillis();
        List<Map<Double, Integer>> massTimeMapList = new LinkedList<>();
        try {
            GRBModel model = new GRBModel(env);
            double t = 0.01;
            //obj function
            GRBLinExpr objFunction = new GRBLinExpr();

            objFunction.addConstant(t);

            //variables
            GRBLinExpr totalNumsOnPepConstr = new GRBLinExpr();
            GRBLinExpr totalMassOnPepConstr = new GRBLinExpr();

            GRBVar fakeProtonVar = model.addVar(-1, 1, 0, GRB.INTEGER, "fakeProton");
//            GRBVar fake001Var = model.addVar(-1, 1, 0, GRB.INTEGER, "fake001");
            totalMassOnPepConstr.addTerm(MassTool.PROTON, fakeProtonVar);
//            totalMassOnPepConstr.addTerm(0.01, fake001Var);

            Map<Double, GRBVar> massVarMap = new HashMap<>();
            massVarMap.put(MassTool.PROTON, fakeProtonVar);
//            massVarMap.put(0.01, fake001Var);
            for (double mass : allMassMaxTimeMap.keySet()) {
                GRBVar tmpVar = model.addVar(0, allMassMaxTimeMap.get(mass), 0, GRB.INTEGER, "isPtmSelected_"+mass);
                massVarMap.put(mass, tmpVar);
                totalMassOnPepConstr.addTerm(mass, tmpVar); // + m_i * x_i
                totalNumsOnPepConstr.addTerm(1,tmpVar); // + 1 * x_i
            }
//            model.addConstr(totalMassOnPepConstr, GRB.EQUAL, Double.parseDouble(df2.format(totalDeltaMass)) , "totalMassOnPepConstr");
            model.addConstr(totalMassOnPepConstr, GRB.GREATER_EQUAL, totalDeltaMass - ms1TolAbs , "constrM1");
            model.addConstr(totalMassOnPepConstr, GRB.LESS_EQUAL, totalDeltaMass +ms1TolAbs, "constrM2"); //or put this to constraints as a model.addRange
            model.addConstr(totalNumsOnPepConstr, GRB.LESS_EQUAL, Math.min(partSeq.length(),4), "totalNumsOnPepConstr"); // this value should not exceed the number of aa in the partSeq

            Set<Double> allMassSet = new HashSet<>(allMassMaxTimeMap.keySet());
            allMassSet.add(MassTool.PROTON);
//            allMassSet.add(0.01);

            int poolSolutions = 1000;
            model.set(GRB.IntParam.MIPFocus, 1);
            model.set(GRB.IntParam.PoolSearchMode, 1);
            model.set(GRB.IntParam.PoolSolutions, poolSolutions);
            model.set(GRB.DoubleParam.TimeLimit, 1); // second

            t2=System.currentTimeMillis();

            model.optimize();
            t3=System.currentTimeMillis();
            switch (model.get(GRB.IntAttr.Status)) {
                case GRB.OPTIMAL:
                    break;
                case GRB.TIME_LIMIT:
                    break;
                default:
                    model.dispose();
                    return new LinkedList<>();
            }
            int solCount = model.get(GRB.IntAttr.SolCount);
            if (solCount == 0) return new LinkedList<>();
            for (int solNum = 0; solNum < solCount; solNum++) {
                model.set(GRB.IntParam.SolutionNumber, solNum);
                Map<Double, Integer> massTimeMap = new HashMap<>();
                for (double mass : allMassSet){
                    GRBVar var = massVarMap.get(mass);
                    int time = (int) Math.round(var.get(GRB.DoubleAttr.Xn));
                    if (time > 0) {
                        massTimeMap.put(mass, time);
                    }
                }
                massTimeMapList.add(massTimeMap);
            }

            for (Map<Double, Integer> map : massTimeMapList) {
                List<Double> massList = new ArrayList<>(map.keySet());
                Collections.sort(massList, Comparator.reverseOrder());
                System.out.println(massList);
            }
            t4=System.currentTimeMillis();
            model.dispose();
        } catch (GRBException e) {
            System.out.println("Error code: " + e.getErrorCode() + ". " + e.getMessage());
        }
        return massTimeMapList;
    }

    private int collectPtmCombs(int lastPeakNum, int scanNum, GRBEnv env, List<Map<Integer, Integer>> allSolsList, double[] bIons,
                                double[] yIons, int numPtmsOnPep, Set<Integer> modifiedZone, Map<Integer, VarPtm[]> idxVarModMap,
                                double totalDeltaMass, Set<Integer> constraintZone, double extraDeltaMass, boolean shouldExtraConstr,
                                double leftMassBound, double rightMassBound, String ptmFreePeptide, boolean isDecoy, double normalizedCrossCorr,
                                int globalRank, ModPepPool peptidePTMPattern, SparseVector expProcessedPL, TreeMap<Double, Double> plMap,
                                int precursorCharge, int localMaxMS2Charge, SortedMap<Double, Double> unUsedPlMap) throws Exception{ // Sometimes, the precursor mass error may affects the digitized spectrum.
        System.out.print("[");
        for (double mz: yIons){
            System.out.print(mz+",");
        }
        System.out.println("]\n");
//        System.out.println(unUsedPlMap);
        double averageMass = totalDeltaMass/numPtmsOnPep;
        double massL = averageMass - 0.01;
        double massR = averageMass + 0.01;
        try {
            GRBModel model = new GRBModel(env);
            double t = 0.01;

            //variables
            Map<Integer, GRBVar[]> posVarsMap = new HashMap<>();
            GRBLinExpr totalPtmsOnPepConstr = new GRBLinExpr();
            GRBLinExpr totalMassOnPepConstr = new GRBLinExpr();
            GRBLinExpr totalMassInConstrZoneConstr = new GRBLinExpr();
            GRBLinExpr massLConstr = new GRBLinExpr();
            GRBLinExpr massRConstr = new GRBLinExpr();
            for (int pos : modifiedZone){
                if (!idxVarModMap.containsKey(pos)) continue;
                int numPtmsOnAA = idxVarModMap.get(pos).length;
                GRBVar[] varsArray = new GRBVar[numPtmsOnAA];
                for (int i = 0; i < numPtmsOnAA; i++){
                    varsArray[i] = model.addVar(0, 1, 0, GRB.BINARY, "B_x_pos"+pos+ptmFreePeptide.charAt(pos)+"_ptmId"+i);
                }
                posVarsMap.put(pos, varsArray);

                GRBLinExpr onePtmOnAaConstr = new GRBLinExpr();
                double[] coeffOneAaArray = new double[numPtmsOnAA];
                Arrays.fill(coeffOneAaArray, 1);
                onePtmOnAaConstr.addTerms(coeffOneAaArray, posVarsMap.get(pos));
                model.addConstr(onePtmOnAaConstr, GRB.LESS_EQUAL, 1, "CON_1AA_Max1Ptm_pos"+pos+ptmFreePeptide.charAt(pos));

                totalPtmsOnPepConstr.addTerms(coeffOneAaArray, posVarsMap.get(pos));

                double[] coeffMassAaArray = new double[numPtmsOnAA];
                for (int i = 0; i < numPtmsOnAA; i++){
                    coeffMassAaArray[i] = idxVarModMap.get(pos)[i].mass;
                    if (coeffMassAaArray[i] < massL){
                        massLConstr.addTerm(1, posVarsMap.get(pos)[i]);
                    }
                    if(coeffMassAaArray[i] > massR) {
                        massRConstr.addTerm(1, posVarsMap.get(pos)[i]);
                    }
                }
                totalMassOnPepConstr.addTerms(coeffMassAaArray, posVarsMap.get(pos));

                if (shouldExtraConstr && constraintZone.contains(pos)) {
                    totalMassInConstrZoneConstr.addTerms(coeffMassAaArray, posVarsMap.get(pos));
                }
            }

            //constraints
//            model.addConstr(totalPtmsOnPepConstr, GRB.GREATER_EQUAL, 1, "CON_totalNumPtmsOnPep");
//            model.addConstr(totalPtmsOnPepConstr, GRB.LESS_EQUAL, 4, "CON_totalNumPtmsOnPep");
            model.addConstr(totalPtmsOnPepConstr, GRB.EQUAL, numPtmsOnPep, "CON_totalNumPtmsOnPep");

//            model.addRange(totalMassOnPepConstr, totalDeltaMass - t, totalDeltaMass + t, "test");
            model.addConstr(totalMassOnPepConstr, GRB.GREATER_EQUAL, totalDeltaMass - t, "CON_totalMassOnPep_GEQ");
            model.addConstr(totalMassOnPepConstr, GRB.LESS_EQUAL, totalDeltaMass + t, "CON_totalMassOnPep_LEQ"); //or put this to constraints as a model.addRange
            if (shouldExtraConstr && extraDeltaMass != -9999d && !constraintZone.isEmpty()) {
                model.addConstr(totalMassInConstrZoneConstr, GRB.GREATER_EQUAL, extraDeltaMass - t, "CON_extraMassInZone_GEQ");
                model.addConstr(totalMassInConstrZoneConstr, GRB.LESS_EQUAL, extraDeltaMass + t, "CON_extraMassInZone_LEQ"); //or put this to constraints as a model.addRange
            }

            model.addConstr(massLConstr, GRB.LESS_EQUAL, numPtmsOnPep-1, "CON_NumPtm_LeftHalf");
            model.addConstr(massRConstr, GRB.LESS_EQUAL, numPtmsOnPep-1, "CON_NumPtm_RightHalf");

            //matching
            double Drec = 1/2500d;
            double tau = 1/(1-0.01*Drec);
            List<Map.Entry<Double, Double>> peakList = new ArrayList<>(unUsedPlMap.entrySet());

            // b ions
            Map<Integer, Map<Integer, GRBVar>> bPosVars = new HashMap<>();
            Map<Integer, List<Integer>> bPosPeakIdMap = new HashMap<>();
            double aaMass = massTool.calResidueMass(ptmFreePeptide.substring(0, Collections.min(modifiedZone))) + massTool.PROTON;
            double lbMass = aaMass;
            double rbMass = aaMass;
            for (int pos : modifiedZone) {
                aaMass += massTool.getMassTable().get(ptmFreePeptide.charAt(pos));
                lbMass += massTool.getMassTable().get(ptmFreePeptide.charAt(pos));
                rbMass += massTool.getMassTable().get(ptmFreePeptide.charAt(pos));
                if (idxVarModMap.containsKey(pos)){
                    if (idxVarModMap.get(pos)[0].mass < 0) lbMass += idxVarModMap.get(pos)[0].mass;
                    if (idxVarModMap.get(pos)[idxVarModMap.get(pos).length-1].mass > 0) rbMass += idxVarModMap.get(pos)[idxVarModMap.get(pos).length-1].mass;
                }
                GRBLinExpr bVarSumAtPos = new GRBLinExpr();
                Map<Integer, GRBVar> varsMap = new HashMap<>();

                List<Integer> peakIdList = new ArrayList<>();
                for (int pId = 0; pId < peakList.size(); pId++){
                    if (peakList.get(pId).getKey() < lbMass - 0.01) continue;
                    if (peakList.get(pId).getKey() > rbMass + 0.01) break;
                    peakIdList.add(pId);
                    GRBVar tempVar = model.addVar(0, 1, 0, GRB.BINARY, "B_b_pos"+pos+ptmFreePeptide.charAt(pos)+"_peak"+pId);
                    varsMap.put(pId, tempVar);
                    bVarSumAtPos.addTerm(1, tempVar);

                    GRBLinExpr bVarMass = new GRBLinExpr();
                    GRBLinExpr bVarMass2 = new GRBLinExpr();
                    bVarMass.addTerm(1, tempVar);
                    bVarMass2.addTerm(1, tempVar);
                    for (int i : modifiedZone) {
                        if (!idxVarModMap.containsKey(i)) continue;
                        if (i > pos) break;
                        for (int ptmId = 0; ptmId < idxVarModMap.get(i).length; ptmId++){
                            bVarMass.addTerm(-Drec*idxVarModMap.get(i)[ptmId].mass, posVarsMap.get(i)[ptmId]);
                            bVarMass2.addTerm(Drec*idxVarModMap.get(i)[ptmId].mass, posVarsMap.get(i)[ptmId]);
                        }
                    }
                    model.addConstr(bVarMass, GRB.LESS_EQUAL, Drec*(bIons[pos]-peakList.get(pId).getKey())+tau, "CON_b_pos"+pos+ptmFreePeptide.charAt(pos)+"_peak"+pId+"_LEQ1");
                    model.addConstr(bVarMass2, GRB.LESS_EQUAL, Drec*(peakList.get(pId).getKey()-bIons[pos])+tau, "CON_b_pos"+pos+ptmFreePeptide.charAt(pos)+"_peak"+pId+"_LEQ2");
                }
                bPosPeakIdMap.put(pos, peakIdList);
                bPosVars.put(pos, varsMap);
                model.addConstr(bVarSumAtPos, GRB.LESS_EQUAL, 1, "CON_b_pos"+pos+ptmFreePeptide.charAt(pos)+"_max1Peak");
            }

            for (int pos1 : modifiedZone) {   //for No cross constraints
                if (pos1 == Collections.max(modifiedZone)) continue;
                int pos2 = pos1 + 1;
                for (int pId1 : bPosPeakIdMap.get(pos1)) {
                    if (bPosPeakIdMap.get(pos2) == null) {
                        int a = 1;
                        System.out.println(scanNum);
                    }
                    if (!bPosPeakIdMap.get(pos2).contains(pId1)) continue;
                    for (int pId2 : bPosPeakIdMap.get(pos2)) {
                        if (peakList.get(pId2).getKey() - peakList.get(pId1).getKey() > 57.05) break;
                        GRBLinExpr cross = new GRBLinExpr();
                        cross.addTerm(1, bPosVars.get(pos1).get(pId1));
                        cross.addTerm(1, bPosVars.get(pos2).get(pId2));
                        model.addConstr(cross, GRB.LESS_EQUAL ,1, "CON_noCross_bpos1_"+pos1+ptmFreePeptide.charAt(pos1)+"toPeak"+pId1+"_bpos2_"+pos2+ptmFreePeptide.charAt(pos2)+"toPeak"+pId2);
                    }
                }
            }
            // y ions
            Map<Integer, Map<Integer, GRBVar>> yPosVars = new HashMap<>();
            Map<Integer, List<Integer>> yPosPeakIdMap = new HashMap<>();
            double aaMassYion = massTool.calResidueMass(ptmFreePeptide.substring(1+Collections.max(modifiedZone)))+ massTool.H2O + massTool.PROTON;
            double lbMassYion = aaMassYion;
            double rbMassYion = aaMassYion;
            for (int pos = Collections.max(modifiedZone); pos >= Collections.min(modifiedZone); pos--) {
                aaMassYion += massTool.getMassTable().get(ptmFreePeptide.charAt(pos));
                lbMassYion += massTool.getMassTable().get(ptmFreePeptide.charAt(pos));
                rbMassYion += massTool.getMassTable().get(ptmFreePeptide.charAt(pos));
                if (idxVarModMap.containsKey(pos)) {
                    if (idxVarModMap.get(pos)[0].mass < 0) lbMassYion += idxVarModMap.get(pos)[0].mass;
                    if (idxVarModMap.get(pos)[idxVarModMap.get(pos).length-1].mass > 0) rbMassYion += idxVarModMap.get(pos)[idxVarModMap.get(pos).length-1].mass;
                }
                GRBLinExpr yVarSumAtPos = new GRBLinExpr();
                Map<Integer, GRBVar> varsMap = new HashMap<>();

                List<Integer> peakIdList = new ArrayList<>();
                for (int pId = 0; pId < peakList.size(); pId++){
                    if (peakList.get(pId).getKey() < lbMassYion - 0.01) continue;
                    if (peakList.get(pId).getKey() > rbMassYion + 0.01) break;
                    peakIdList.add(pId);
                    GRBVar tempVar = model.addVar(0, 1, 0, GRB.BINARY, "B_y_pos"+pos+ptmFreePeptide.charAt(pos)+"_peak"+pId);
                    varsMap.put(pId, tempVar);
                    yVarSumAtPos.addTerm(1, tempVar);

                    GRBLinExpr yVarMass = new GRBLinExpr();
                    GRBLinExpr yVarMass2 = new GRBLinExpr();
                    yVarMass.addTerm(1, tempVar);
                    yVarMass2.addTerm(1, tempVar);
                    for (int i : modifiedZone) {
                        if (i < pos) continue;
                        if (!idxVarModMap.containsKey(i)) continue;
                        for (int ptmId = 0; ptmId < idxVarModMap.get(i).length; ptmId++){
                            yVarMass.addTerm(-Drec*idxVarModMap.get(i)[ptmId].mass, posVarsMap.get(i)[ptmId]);
                            yVarMass2.addTerm(Drec*idxVarModMap.get(i)[ptmId].mass, posVarsMap.get(i)[ptmId]);
                        }
                    }
                    model.addConstr(yVarMass, GRB.LESS_EQUAL, Drec*(yIons[pos-1]-peakList.get(pId).getKey())+tau, "CON_y_pos"+pos+ptmFreePeptide.charAt(pos)+"_peak"+pId+"_LEQ1");
                    model.addConstr(yVarMass2, GRB.LESS_EQUAL, Drec*(peakList.get(pId).getKey()-yIons[pos-1])+tau, "CON_y_pos"+pos+ptmFreePeptide.charAt(pos)+"_peak"+pId+"_LEQ2");
                }
                yPosPeakIdMap.put(pos, peakIdList);
                model.addConstr(yVarSumAtPos, GRB.LESS_EQUAL, 1, "CON_y_pos"+pos+ptmFreePeptide.charAt(pos)+"_max1Peak");
                yPosVars.put(pos, varsMap);
            }
//            System.out.println(yPosPeakIdMap);

            for (int pos1 = Collections.max(modifiedZone); pos1 >= Collections.min(modifiedZone)+1; pos1--) {
                int pos2 = pos1 - 1;
                for (int pId1 : yPosPeakIdMap.get(pos1)) {
                    if (!yPosPeakIdMap.get(pos2).contains(pId1)) continue;
                    for (int pId2 : yPosPeakIdMap.get(pos2)) {
                        if (peakList.get(pId2).getKey() - peakList.get(pId1).getKey() > 57.05) break;
                        GRBLinExpr cross = new GRBLinExpr();
                        cross.addTerm(1, yPosVars.get(pos1).get(pId1));
                        cross.addTerm(1, yPosVars.get(pos2).get(pId2));
                        model.addConstr(cross, GRB.LESS_EQUAL ,1, "CON_noCross_ypos1_"+pos1+ptmFreePeptide.charAt(pos1)+"toPeak"+pId1+"_ypos2_"+pos2+ptmFreePeptide.charAt(pos2)+"toPeak"+pId2);
                    }
                }
            }

            // end of matching

            GRBLinExpr numMatchedPeaks = new GRBLinExpr();
            for (int pos : bPosPeakIdMap.keySet()) {
                for (int pId : bPosPeakIdMap.get(pos)) {
                    numMatchedPeaks.addTerm(1, bPosVars.get(pos).get(pId));
                }
            }
            for (int pos : yPosPeakIdMap.keySet()){
                for (int pId : yPosPeakIdMap.get(pos)){
                    numMatchedPeaks.addTerm(1, yPosVars.get(pos).get(pId));
                }
            }
            model.addConstr(numMatchedPeaks, GRB.GREATER_EQUAL, lastPeakNum+1, "test");

            //obj
//            GRBLinExpr objFunction = new GRBLinExpr();
//            for (int pos = 1; pos <= ptmFreePeptide.length()-3; pos++) {
//                for (int pId : bPosPeakIdMap.get(pos)) {
//                    objFunction.addTerm(peakList.get(pId).getValue(), bPosVars.get(pos).get(pId));
//                }
//                for (int pId : yPosPeakIdMap.get(pos)){
//                    objFunction.addTerm(peakList.get(pId).getValue(), yPosVars.get(pos).get(pId));
//                }
//            }
//            for (int pId : yPosPeakIdMap.get(ptmFreePeptide.length()-2)){
//                objFunction.addTerm(peakList.get(pId).getValue(), yPosVars.get(ptmFreePeptide.length()-2).get(pId));
//            }
////            }
//
//            for (int pos: posVarsMap.keySet()){
//                for (GRBVar var : posVarsMap.get(pos)){
//                    objFunction.addTerm(-0.1, var);
//                }
//            }
            model.setObjective(numMatchedPeaks, GRB.MAXIMIZE);

            int solId = 0;
            for (Map<Integer, Integer> sol : allSolsList){
                GRBLinExpr solConstr = new GRBLinExpr();
                for (int aId : sol.keySet()){
                    solConstr.addTerm(1, posVarsMap.get(aId)[sol.get(aId)]);
                }
                model.addConstr(solConstr, GRB.LESS_EQUAL, numPtmsOnPep-1, "sol_"+solId);
                solId++;
            }
            allSolsList.clear();

            //settings
            model.set(GRB.IntParam.MIPFocus, 1);
//            model.set(GRB.DoubleParam.TimeLimit, 1.5*numPtmsOnPep);
            // solve the model
//            model.write("/home/slaiad/Data/Simulation_Data/test.lp");
//            model.write("/home/slaiad/Data/Simulation_Data/test.mps");
            int numOfSols = 0;
            int maxLastNPeaks = 0;

            whileLoop:
            while ( numOfSols < 1 ) {//10 - 2*numPtmsOnPep
                long start=System.currentTimeMillis();   //获取开始时间
                model.optimize();
                long end=System.currentTimeMillis(); //获取结束时间
                switch (model.get(GRB.IntAttr.Status)) {
                    case GRB.OPTIMAL:
//                        model.write("/home/slaiad/Data/Simulation_Data/test.sol");

                        System.out.println("\nNew sol:"+numOfSols +"," + numPtmsOnPep+ ", "+(end-start)+"ms");
                        break;//                        System.out.println(scanNum + " :TimeOut");
                    default:
                        model.dispose();
                        break whileLoop;
                }

                numOfSols++;
                Map<Integer, Integer> positionOneMap = new HashMap<>();
//                positionOneMap.clear();
                for (int pos : modifiedZone) {
                    if (!posVarsMap.containsKey(pos)) continue;
                    GRBVar[] varsResArray = posVarsMap.get(pos);
                    for (int varId = 0; varId < varsResArray.length; varId++ ) {  //should I be careful about the binary variable to be really 0/1 instead of decimal
                        if (1 == Math.round(varsResArray[varId].get(GRB.DoubleAttr.X))) {
                            positionOneMap.put(pos, varId);
                        }
                    }
                }
                allSolsList.add(positionOneMap);

                //stub test true other
//                positionOneMap.clear();
//                positionOneMap.put(1,5);
//                positionOneMap.put(10,13);
//                positionOneMap.put(12,1);
                //stub end

                PosMassMap positionDeltaMassMap = new PosMassMap(ptmFreePeptide.length());
                Peptide peptide = new Peptide(ptmFreePeptide, isDecoy, massTool);
                for (Map.Entry<Integer, Integer> entry : positionOneMap.entrySet()){
                    positionDeltaMassMap.put(entry.getKey(), idxVarModMap.get(entry.getKey())[entry.getValue()].mass);
                }

//                for (Integer coor : lastPeptide.getVarPTMs().keySet()) {
//                    newPosMassMap.put(coor, lastPeptide.getVarPTMs().get(coor)); // copy the ptms from lastPeptide
//                }

                peptide.setVarPTM(positionDeltaMassMap);

                double[][] temp = peptide.getIonMatrix();
                double score = massTool.buildVectorAndCalXCorr(peptide.getIonMatrix(), precursorCharge, expProcessedPL) - 0.01*(numPtmsOnPep-1);
//                double lastLpScore = 0;
                int currentNPeaks = 0;
                if (score > 0) {
//                    peptide.lpScore = model.get(GRB.DoubleAttr.ObjVal);
//                    System.out.println("numPtmsOnPep "+numPtmsOnPep+" lpScore "+peptide.lpScore + " XCorr " + score);
                    peptide.setScore(score);
                    peptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, localMaxMS2Charge, peptide.getIonMatrix(), ms2Tolerance));
                    peptidePTMPattern.push(peptide);
                    for (int pos : bPosVars.keySet()) {
                        for (int peakId : bPosVars.get(pos).keySet()){
                            if (1 == Math.round(bPosVars.get(pos).get(peakId).get(GRB.DoubleAttr.X))){
//                                lastLpScore += peakList.get(peakId).getValue();
                                currentNPeaks++;
                                System.out.println("bMatch "+bPosVars.get(pos).get(peakId).get(GRB.StringAttr.VarName) + ", "+ peakList.get(peakId).getKey()+","+peakList.get(peakId).getValue());
                            }
                        }
                    }
                    for (int pos : yPosVars.keySet()) {
                        for (int peakId : yPosVars.get(pos).keySet()){
                            if (1 == Math.round(yPosVars.get(pos).get(peakId).get(GRB.DoubleAttr.X))){
//                                lastLpScore += peakList.get(peakId).getValue();
                                currentNPeaks++;
                                System.out.println("yMatch "+yPosVars.get(pos).get(peakId).get(GRB.StringAttr.VarName) + ", "+ peakList.get(peakId).getKey()+","+peakList.get(peakId).getValue());
                            }
                        }
                    }
//                    System.out.println("this  "+ lastLpScore + "," + currentNPeaks + "," +peptide);

                    List<Double> massList = new ArrayList<>();
                    for (double mass : positionDeltaMassMap.values()){
                        massList.add(mass);
                    }
                    Collections.sort(massList);
//                    System.out.println(numOfSols + ", " + score + ", "+peptide.getMatchedPeakNum()+", "+ massList);
                }

                GRBLinExpr forMoreFeasiSolConstr = new GRBLinExpr();
                for (Map.Entry<Integer, Integer> entry : positionOneMap.entrySet()){
                    forMoreFeasiSolConstr.addTerm(1, posVarsMap.get(entry.getKey())[entry.getValue()]);
                }
//                model.addConstr(forMoreFeasiSolConstr, GRB.LESS_EQUAL, numPtmsOnPep - 0.8, "forFeasiSol"); //or put this to constraints as a model.addRange

                GRBLinExpr numMatchedPeaks1 = new GRBLinExpr();
                for (int pos : bPosPeakIdMap.keySet()) {
                    for (int pId : bPosPeakIdMap.get(pos)) {
                        numMatchedPeaks1.addTerm(1, bPosVars.get(pos).get(pId));
                    }
                }
                for (int pos : yPosPeakIdMap.keySet()){
                    for (int pId : yPosPeakIdMap.get(pos)){
                        numMatchedPeaks1.addTerm(1, yPosVars.get(pos).get(pId));
                    }
                }
                model.addConstr(numMatchedPeaks1, GRB.GREATER_EQUAL, currentNPeaks + 1, "test");
                maxLastNPeaks = currentNPeaks;
            }
            model.dispose();
            return maxLastNPeaks;
        } catch (GRBException e) {
            System.out.println("Error code: " + e.getErrorCode() + ". " + e.getMessage());
        }
        return 0;
    }
    private Set<Integer> getInitModZoneComple(String freeSeq, boolean isDecoy, MassTool massTool, int localMaxMS2Charge,
                                        double tagVecScore, int globalRank, SparseVector expProcessedPL, TreeMap<Double, Double> plMap){

        Peptide cleanPep = new Peptide(freeSeq, isDecoy, massTool);
        int lb = 0;  //lb included
        int rb = freeSeq.length() - 1;//rb included
        Map<Integer, Double> matchedBions = new HashMap<>();
        Map<Integer, Double> matchedYions = new HashMap<>();
        double[][] ionMatrix = cleanPep.getIonMatrixNow();

        Set<Integer> jRange = IntStream.rangeClosed(0, freeSeq.length()-1).boxed().collect(Collectors.toSet());

        double cleanScore = massTool.buildVectorAndCalXCorr(ionMatrix, 1, expProcessedPL, matchedBions, matchedYions, jRange) ;
        cleanPep.setScore(cleanScore);
        cleanPep.setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, cleanPep.getIonMatrix(), ms2Tolerance));
        cleanPep.matchedBions.putAll(matchedBions);
        cleanPep.matchedYions.putAll(matchedYions);

        List<Integer> matchedBIndex = new ArrayList<>(matchedBions.keySet());
        List<Integer> matchedYIndex = new ArrayList<>(matchedYions.keySet());
        Collections.sort(matchedBIndex);
        Collections.sort(matchedYIndex);
        if (matchedBIndex.size() >= 2) {
            if (matchedBIndex.get(matchedBIndex.size()-1) + (1-matchedBions.get(matchedBIndex.get(matchedBIndex.size()-1))) > matchedBIndex.get(matchedBIndex.size()-2)+4) {
                matchedBions.remove(matchedBIndex.get(matchedBIndex.size()-1));
            }
        }
        if (matchedYIndex.size() >= 2) {
            if (matchedYIndex.get(0) - (1-matchedYions.get(matchedYIndex.get(0))) < matchedYIndex.get(1)-4) {
                matchedYions.remove(matchedYIndex.get(0));
            }
        }

        if (matchedBions.size() > 1) {  // should has at least two peaks to be trusted
            lb = Collections.max(matchedBions.keySet()) + 1;
        }
        if (matchedYions.size() > 1) {
            rb = Collections.min(matchedYions.keySet()) - 1;
        }

        if (rb - lb <= 0) {
            double bSumIntens = 0;
            for (double intes : matchedBions.values()) bSumIntens += intes;
            double ySumIntens = 0;
            for (double intes : matchedYions.values()) ySumIntens += intes;
            if (bSumIntens > ySumIntens) {
                rb = freeSeq.length() - 1;
            } else {
                lb = 0;
            }
        }
        if (rb < lb) {
            rb = freeSeq.length() - 1;
            lb = 0;
        }
        return IntStream.rangeClosed(lb, rb).boxed().collect(Collectors.toSet());
    }
    private Set<Integer> getInitModZone(String freeSeq, boolean isDecoy, MassTool massTool, int localMaxMS2Charge,
                                        double tagVecScore, int globalRank, SparseVector expProcessedPL, TreeMap<Double, Double> plMap, double cutMass, double ncPart){

        Peptide cleanPep = new Peptide(freeSeq, isDecoy, massTool);
        int lb = 0;  //lb included
        int rb = freeSeq.length() - 1;//rb included
        Map<Integer, Double> matchedBions = new HashMap<>();
        Map<Integer, Double> matchedYions = new HashMap<>();
        double[][] ionMatrix = cleanPep.getIonMatrixNow();

        if (ncPart == N_PART) { // is n part seq
            for (int i = 0; i < ionMatrix[1].length; i++) {
                ionMatrix[1][i] += cutMass;
            }
        } else {            //is c part seq
            for (int i = 0; i < ionMatrix[0].length; i++) {
                ionMatrix[0][i] += cutMass;
            }
        }

        Set<Integer> jRange = IntStream.rangeClosed(0, freeSeq.length()-1).boxed().collect(Collectors.toSet());

        double cleanScore = massTool.buildVectorAndCalXCorr(ionMatrix, 1, expProcessedPL, matchedBions, matchedYions, jRange) ;
        cleanPep.setScore(cleanScore);
        cleanPep.setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, cleanPep.getIonMatrix(), ms2Tolerance));
        cleanPep.matchedBions.putAll(matchedBions);
        cleanPep.matchedYions.putAll(matchedYions);

        List<Integer> matchedBIndex = new ArrayList<>(matchedBions.keySet());
        List<Integer> matchedYIndex = new ArrayList<>(matchedYions.keySet());
        Collections.sort(matchedBIndex);
        Collections.sort(matchedYIndex);
        if (matchedBIndex.size() >= 2) {
            if (matchedBIndex.get(matchedBIndex.size()-1) + (1-matchedBions.get(matchedBIndex.get(matchedBIndex.size()-1))) > matchedBIndex.get(matchedBIndex.size()-2)+4) {
                matchedBions.remove(matchedBIndex.get(matchedBIndex.size()-1));
            }
        }
        if (matchedYIndex.size() >= 2) {
            if (matchedYIndex.get(0) - (1-matchedYions.get(matchedYIndex.get(0))) < matchedYIndex.get(1)-4) {
                matchedYions.remove(matchedYIndex.get(0));
            }
        }

        if (matchedBions.size() > 1) {  // should has at least two peaks to be trusted
            lb = Collections.max(matchedBions.keySet()) + 1;
        }
        if (matchedYions.size() > 1) {
            rb = Collections.min(matchedYions.keySet()) - 1;
        }

        if (rb - lb <= 0) {
            double bSumIntens = 0;
            for (double intes : matchedBions.values()) bSumIntens += intes;
            double ySumIntens = 0;
            for (double intes : matchedYions.values()) ySumIntens += intes;
            if (bSumIntens > ySumIntens) {
                rb = freeSeq.length() - 1;
            } else {
                lb = 0;
            }
        }
        if (rb < lb) {
            rb = freeSeq.length() - 1;
            lb = 0;
        }
        if (ncPart == N_PART) { // if this is nPart and there is n-term enriched ptm, then lb of modzone on npart must be 0
            lb = 0;
        }
        return IntStream.rangeClosed(lb, rb).boxed().collect(Collectors.toSet());
    }

    private void updateIonMatrix(double [][] ionMatrix, double cutMass, byte ncPart){
        if (ncPart == N_PART) { // is n part seq
            for (int i = 0; i < ionMatrix[1].length; i++) {
                ionMatrix[1][i] += cutMass;
            }
        } else {            //is c part seq
            for (int i = 0; i < ionMatrix[0].length; i++) {
                ionMatrix[0][i] += cutMass;
            }
        }
    }
//    public ModPepPool settlePtmOnSide(int scanNum, SparseVector expProcessedPL, TreeMap<Double, Double> plMap, double precursorMass, String partSeq, boolean isDecoy,
//                                       Map<Integer, VarPtm[]> posVarPtmArraySrcMap, double cutMass, double deltaMass, int precursorCharge, byte ncPart, double ms1TolAbs) throws CloneNotSupportedException {
//        int localMaxMS2Charge = 1;
//        double tagVecScore = -0.99;
//        int globalRank = -1;
//
//        TreeMap<Double, Double> unUsedPlMap = new TreeMap<>(plMap);
//        Peptide lastPeptide = new Peptide(partSeq, isDecoy, massTool);
//        lastPeptide.setVarPTM(new PosMassMap(partSeq.length()));
////        lastPeptide.shouldPTM = true;
//
//        Set<Integer> modifiedZone = getInitModZone(partSeq, isDecoy, massTool, localMaxMS2Charge, tagVecScore, globalRank, expProcessedPL, plMap, cutMass, ncPart);// e.g. nABCDEFc, modifiedZone={1,2,3,4,5,6}
////        Set<Integer> modifiedZone = IntStream.rangeClosed(0, partSeq.length()-1).boxed().collect(Collectors.toSet());// e.g. nABCDEFc, modifiedZone={1,2,3,4,5,6}
//        ModPepPool modPepPoolGood = new ModPepPool(partSeq, 10);
//        ModPepPool modPepPoolBad = new ModPepPool(partSeq, 10);
//
//        double massToSettle = deltaMass;
//        Set<Integer> toModZone = new HashSet<>(modifiedZone);
//        double finalUnsettledMass = 0;
//        for (int loop = 1; loop <= 2; loop++) {
//            DividedZone dividedZone = dividePepNew(scanNum, toModZone, modPepPoolGood, modPepPoolBad, lastPeptide, posVarPtmArraySrcMap, massToSettle, cutMass, ncPart, partSeq, isDecoy,
//                    tagVecScore, globalRank, expProcessedPL, plMap, precursorCharge, 1, unUsedPlMap, 0);
//
//            finalUnsettledMass = massToSettle - dividedZone.ptmMass;
//            if (modPepPoolBad.peptideTreeSet.isEmpty() ) { //fixme, why the first sentence
//                break;
//            }
//            toModZone = dividedZone.toModZone;
//            double massSettled = dividedZone.ptmMass;
//            massToSettle -= massSettled;
//            lastPeptide = modPepPoolBad.getTopPepPtn();
//            modPepPoolBad = new ModPepPool(partSeq, 1);
//            if (dividedZone.toModZone.isEmpty() || loop == 2) {
//                break;
//            }
//        }
//
//        // two make up strategy if there is no good pool found
////        if (modPepPoolGood.peptideTreeSet.isEmpty() && ncPart == N_PART && posVarPtmArraySrcMap.containsKey(0)) {
//        if (ncPart == N_PART && posVarPtmArraySrcMap.containsKey(0)) { //what if I dont limit only good for dimethyl label dataset
//            // Assign pepN or protN varPtm to the nterm aa and try whether it can settle all the mass
//
//            for (VarPtm nVarPtm : posVarPtmArraySrcMap.get(0)) {
//                if (nVarPtm.priority != 1) continue;
//                ModPepPool tmpModPepPoolBad = new ModPepPool(partSeq, 1);
//                double tmpMassToSettle = deltaMass-nVarPtm.mass;
//                Set<Integer> tmpToModZone = new HashSet<>(modifiedZone);
//                tmpToModZone.remove(0);
//                Peptide tmpPeptide = new Peptide(partSeq, isDecoy, massTool);
//                PosMassMap tmpNewPosMassMap = new PosMassMap(partSeq.length());
//                tmpNewPosMassMap.put(0, nVarPtm.mass);
//                tmpPeptide.posVarPtmResMap.put(0, nVarPtm);
//                tmpPeptide.setVarPTM(tmpNewPosMassMap);
//                double[][] ionMatrix = tmpPeptide.getIonMatrixNow();
//                updateIonMatrix(ionMatrix, cutMass, ncPart);
//
//                Set<Integer> jRange = IntStream.rangeClosed(0, partSeq.length()-1).boxed().collect(Collectors.toSet());
//                double tmpScore = massTool.buildVectorAndCalXCorr(ionMatrix, 1, expProcessedPL, tmpPeptide.matchedBions, tmpPeptide.matchedYions, jRange);
//                tmpPeptide.setScore(tmpScore); // init score needs to be bring in the 2-loop
//
//                for (int loop = 1; loop <= 1; loop++) {
//                    DividedZone dividedZone = dividePepNew(scanNum, tmpToModZone, modPepPoolGood, tmpModPepPoolBad, tmpPeptide, posVarPtmArraySrcMap, tmpMassToSettle, cutMass, ncPart, partSeq, isDecoy,
//                            tagVecScore, globalRank, expProcessedPL, plMap, precursorCharge, 1, unUsedPlMap, 0);
//                    if (tmpModPepPoolBad.peptideTreeSet.isEmpty() || dividedZone.toModZone.isEmpty() || loop == 2) { //fixme, why the first sentence
//                        break;
//                    }
//                    tmpToModZone = dividedZone.toModZone;
//                    double massSettled = dividedZone.ptmMass;
////                    dividedZone.ptmMass;
//                    tmpMassToSettle -= massSettled;
//                    tmpPeptide = tmpModPepPoolBad.getTopPepPtn();
//                    tmpModPepPoolBad = new ModPepPool(partSeq, 1);
//                }
//            }
//        }
//        if (modPepPoolGood.peptideTreeSet.isEmpty() && (massToSettle > minPtmMass && massToSettle < maxPtmMass)) {
//            //  just assign the remaining massToSettle to any aa in toModZone, and record the best one
//            double testMass = finalUnsettledMass;
//            if (Math.abs(testMass-1*MassTool.PROTON) < ms1TolAbs
//                    || Math.abs(testMass+1*MassTool.PROTON) < ms1TolAbs ) {
//                modPepPoolGood.push(lastPeptide);
//                return modPepPoolGood;
//            }
//
//            for (int pos : toModZone) { // here must use last toModZone and last massToSettle
//                VarPtm fakeVarPtm = new VarPtm(massToSettle, partSeq.charAt(pos), 4, String.format("PIPI_%s", massToSettle), "PIPI_unsettled", -1);
//                Peptide fakePeptide = lastPeptide.clone();
////                        Peptide fakePeptide = new Peptide(partSeq, isDecoy, massTool, 1, tagVecScore, globalRank);
//                PosMassMap fakeNewPosMassMap = new PosMassMap(partSeq.length());
////                        if (!modPepPoolBad.peptideTreeSet.isEmpty()){
//                for (Integer coor : lastPeptide.getVarPTMs().keySet()) {
//                    fakeNewPosMassMap.put(coor, lastPeptide.getVarPTMs().get(coor)); // copy the ptms from modPepPoolBad
//                }
////                        }
//                fakeNewPosMassMap.put(pos, massToSettle);
//                fakePeptide.posVarPtmResMap.put(pos, fakeVarPtm);
//                fakePeptide.setVarPTM(fakeNewPosMassMap);
//
//                Map<Integer, Double> matchedBions = new HashMap<>();
//                Map<Integer, Double> matchedYions = new HashMap<>();
//                double[][] ionMatrix = fakePeptide.getIonMatrixNow();
//                updateIonMatrix(ionMatrix, cutMass, ncPart);
//
//                int lb = Math.max(0, Collections.min(toModZone));
//                int rb = Math.min(partSeq.length() - 1, Collections.max(toModZone));
//                Set<Integer> jRange = IntStream.rangeClosed(lb, rb).boxed().collect(Collectors.toSet());
//                double fakeScore = massTool.buildVectorAndCalXCorr(ionMatrix, 1, expProcessedPL, matchedBions, matchedYions, jRange);
//                fakePeptide.setScore(fakeScore);
//                fakePeptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, fakePeptide.getIonMatrix(), ms2Tolerance));
//                fakePeptide.matchedBions.putAll(matchedBions);
//                fakePeptide.matchedYions.putAll(matchedYions);
//                modPepPoolGood.push(fakePeptide);
//            }
//        }
//
////        GRBEnv env = new GRBEnv(true);
////        env.set(GRB.IntParam.OutputFlag,0);
////        env.start();
////        List<Map<Integer, Integer>> allSolsList = new ArrayList<>();
////        int numOfSols = checkCorrectMass(scanNum, env, posVarPtmArraySrcMap.keySet(), posVarPtmArraySrcMap, deltaMass, partSeq , allSolsList);
//
//        return modPepPoolGood;
//    }

    public ModPepPool settlePtmOnSide(int scanNum, SparseVector expProcessedPL, TreeMap<Double, Double> plMap, double precursorMass, String partSeq, boolean isDecoy,
                                      Map<Double, Integer> allMassMaxTimeMap, double cutMass, double deltaMass, int precursorCharge, byte ncPart, double ms1TolAbs) throws CloneNotSupportedException, GRBException {
        int localMaxMS2Charge = 1;
        double tagVecScore = -0.99;
        int globalRank = -1;

        TreeMap<Double, Double> unUsedPlMap = new TreeMap<>(plMap);

        GRBEnv env = new GRBEnv();

        findBestPtmMIP(scanNum, env, allMassMaxTimeMap, deltaMass, partSeq, cutMass, ncPart, expProcessedPL, ms1TolAbs);

        Peptide lastPeptide = new Peptide(partSeq, isDecoy, massTool);
        lastPeptide.setVarPTM(new PosMassMap(partSeq.length()));

        Set<Integer> modifiedZone = getInitModZone(partSeq, isDecoy, massTool, localMaxMS2Charge, tagVecScore, globalRank, expProcessedPL, plMap, cutMass, ncPart);// e.g. nABCDEFc, modifiedZone={1,2,3,4,5,6}
        ModPepPool modPepPoolGood = new ModPepPool(partSeq, 10);
        ModPepPool modPepPoolBad = new ModPepPool(partSeq, 10);
//
//        double massToSettle = deltaMass;
//        Set<Integer> toModZone = new HashSet<>(modifiedZone);
//        double finalUnsettledMass = 0;
//        for (int loop = 1; loop <= 2; loop++) {
//            DividedZone dividedZone = dividePepNew(scanNum, toModZone, modPepPoolGood, modPepPoolBad, lastPeptide, allMassMaxTimeMap, massToSettle, cutMass, ncPart, partSeq, isDecoy,
//                    tagVecScore, globalRank, expProcessedPL, plMap, precursorCharge, 1, unUsedPlMap, 0);
//
//            finalUnsettledMass = massToSettle - dividedZone.ptmMass;
//            if (modPepPoolBad.peptideTreeSet.isEmpty() ) { //fixme, why the first sentence
//                break;
//            }
//            toModZone = dividedZone.toModZone;
//            double massSettled = dividedZone.ptmMass;
//            massToSettle -= massSettled;
//            lastPeptide = modPepPoolBad.getTopPepPtn();
//            modPepPoolBad = new ModPepPool(partSeq, 1);
//            if (dividedZone.toModZone.isEmpty() || loop == 2) {
//                break;
//            }
//        }
//
//        // two make up strategy if there is no good pool found
//        if (ncPart == N_PART && allMassMaxTimeMap.containsKey(0)) { //what if I dont limit only good for dimethyl label dataset
//            // Assign pepN or protN varPtm to the nterm aa and try whether it can settle all the mass
//
//            for (VarPtm nVarPtm : allMassMaxTimeMap.get(0)) {
//                if (nVarPtm.priority != 1) continue;
//                ModPepPool tmpModPepPoolBad = new ModPepPool(partSeq, 1);
//                double tmpMassToSettle = deltaMass-nVarPtm.mass;
//                Set<Integer> tmpToModZone = new HashSet<>(modifiedZone);
//                tmpToModZone.remove(0);
//                Peptide tmpPeptide = new Peptide(partSeq, isDecoy, massTool);
//                PosMassMap tmpNewPosMassMap = new PosMassMap(partSeq.length());
//                tmpNewPosMassMap.put(0, nVarPtm.mass);
//                tmpPeptide.posVarPtmResMap.put(0, nVarPtm);
//                tmpPeptide.setVarPTM(tmpNewPosMassMap);
//                double[][] ionMatrix = tmpPeptide.getIonMatrixNow();
//                updateIonMatrix(ionMatrix, cutMass, ncPart);
//
//                Set<Integer> jRange = IntStream.rangeClosed(0, partSeq.length()-1).boxed().collect(Collectors.toSet());
//                double tmpScore = massTool.buildVectorAndCalXCorr(ionMatrix, 1, expProcessedPL, tmpPeptide.matchedBions, tmpPeptide.matchedYions, jRange);
//                tmpPeptide.setScore(tmpScore); // init score needs to be bring in the 2-loop
//
//                for (int loop = 1; loop <= 1; loop++) {
//                    DividedZone dividedZone = dividePepNew(scanNum, tmpToModZone, modPepPoolGood, tmpModPepPoolBad, tmpPeptide, allMassMaxTimeMap, tmpMassToSettle, cutMass, ncPart, partSeq, isDecoy,
//                            tagVecScore, globalRank, expProcessedPL, plMap, precursorCharge, 1, unUsedPlMap, 0);
//                    if (tmpModPepPoolBad.peptideTreeSet.isEmpty() || dividedZone.toModZone.isEmpty() || loop == 2) { //fixme, why the first sentence
//                        break;
//                    }
//                    tmpToModZone = dividedZone.toModZone;
//                    double massSettled = dividedZone.ptmMass;
//                    tmpMassToSettle -= massSettled;
//                    tmpPeptide = tmpModPepPoolBad.getTopPepPtn();
//                    tmpModPepPoolBad = new ModPepPool(partSeq, 1);
//                }
//            }
//        }
//        if (modPepPoolGood.peptideTreeSet.isEmpty() && (massToSettle > minPtmMass && massToSettle < maxPtmMass)) {
//            double testMass = finalUnsettledMass;
//            if (Math.abs(testMass-1*MassTool.PROTON) < ms1TolAbs
//                    || Math.abs(testMass+1*MassTool.PROTON) < ms1TolAbs ) {
//                modPepPoolGood.push(lastPeptide);
//                return modPepPoolGood;
//            }
//
//            for (int pos : toModZone) { // here must use last toModZone and last massToSettle
//                VarPtm fakeVarPtm = new VarPtm(massToSettle, partSeq.charAt(pos), 4, String.format("PIPI_%s", massToSettle), "PIPI_unsettled", -1);
//                Peptide fakePeptide = lastPeptide.clone();
//                PosMassMap fakeNewPosMassMap = new PosMassMap(partSeq.length());
//                for (Integer coor : lastPeptide.getVarPTMs().keySet()) {
//                    fakeNewPosMassMap.put(coor, lastPeptide.getVarPTMs().get(coor)); // copy the ptms from modPepPoolBad
//                }
//                fakeNewPosMassMap.put(pos, massToSettle);
//                fakePeptide.posVarPtmResMap.put(pos, fakeVarPtm);
//                fakePeptide.setVarPTM(fakeNewPosMassMap);
//
//                Map<Integer, Double> matchedBions = new HashMap<>();
//                Map<Integer, Double> matchedYions = new HashMap<>();
//                double[][] ionMatrix = fakePeptide.getIonMatrixNow();
//                updateIonMatrix(ionMatrix, cutMass, ncPart);
//
//                int lb = Math.max(0, Collections.min(toModZone));
//                int rb = Math.min(partSeq.length() - 1, Collections.max(toModZone));
//                Set<Integer> jRange = IntStream.rangeClosed(lb, rb).boxed().collect(Collectors.toSet());
//                double fakeScore = massTool.buildVectorAndCalXCorr(ionMatrix, 1, expProcessedPL, matchedBions, matchedYions, jRange);
//                fakePeptide.setScore(fakeScore);
//                fakePeptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, fakePeptide.getIonMatrix(), ms2Tolerance));
//                fakePeptide.matchedBions.putAll(matchedBions);
//                fakePeptide.matchedYions.putAll(matchedYions);
//                modPepPoolGood.push(fakePeptide);
//            }
//        }
        return modPepPoolGood;
    }

    private DividedZone dividePepNew(int scanNum, Set<Integer> modifiedZone, ModPepPool modPepPoolGood, ModPepPool modPepPoolBad, Peptide lastPeptide, Map<Integer, VarPtm[]> posVarPtmArraySrcMap,
                                     double totalDeltaMass, double cutMass, byte ncPart, String ptmFreePeptide, boolean isDecoy, double normalizedCrossCorr, int globalRank,
                                     SparseVector expProcessedPL, TreeMap<Double, Double> plMap, int precursorCharge, int localMaxMS2Charge, SortedMap<Double, Double> unUsedPlMap, int startRelPos) throws CloneNotSupportedException {
        // Sometimes, the precursor mass error may affects the digitized spectrum.
        double ptmMass = 0;
        for (int pos : modifiedZone) {
            if (!posVarPtmArraySrcMap.containsKey(pos)) continue;

            for (int ptmId = 0; ptmId < posVarPtmArraySrcMap.get(pos).length; ptmId++){
                Peptide peptide = lastPeptide.clone();
                PosMassMap newPosMassMap = new PosMassMap(ptmFreePeptide.length());
                for (Integer coor : lastPeptide.getVarPTMs().keySet()) {
                    newPosMassMap.put(coor, lastPeptide.getVarPTMs().get(coor)); // copy the ptms from lastPeptide
                }
                newPosMassMap.put(pos, posVarPtmArraySrcMap.get(pos)[ptmId].mass);
                peptide.setVarPTM(newPosMassMap);
                int thisPriority = posVarPtmArraySrcMap.get(pos)[ptmId].priority;


                double[][] ionMatrix = peptide.getIonMatrixNow();
                updateIonMatrix(ionMatrix, cutMass, ncPart);

//                int lb = Math.max(0, Collections.min(modifiedZone));
//                int rb = Math.min(ptmFreePeptide.length()-1, Collections.max(modifiedZone)); //
//                Set<Integer> jRange = IntStream.rangeClosed(lb, rb).boxed().collect(Collectors.toSet());
//                double score = massTool.buildVectorAndCalXCorr(ionMatrix, 1, expProcessedPL, matchedBions, matchedYions, jRange) ;
//                double scoreLb = 0;

                Set<Integer> jRange = IntStream.rangeClosed(0, ptmFreePeptide.length()-1).boxed().collect(Collectors.toSet());
                double score = massTool.buildVectorAndCalXCorr(ionMatrix, 1, expProcessedPL, peptide.matchedBions, peptide.matchedYions, jRange) ;
                double scoreLb = lastPeptide.getScore();

                if (Math.abs(totalDeltaMass - posVarPtmArraySrcMap.get(pos)[ptmId].mass) <= 0.01){
                    scoreLb = lastPeptide.getScore()-1;
                }
                if (thisPriority == 1) {
                    scoreLb = -1;
                    score *= 2; //todo change it to that the priority one is only to consider not replacing the original one
                }
                if (score > scoreLb) {
                    peptide.setScore(score);
                    peptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, peptide.getIonMatrix(), ms2Tolerance));
                    peptide.posVarPtmResMap.put(pos, posVarPtmArraySrcMap.get(pos)[ptmId]);
                    if (Math.abs(totalDeltaMass - posVarPtmArraySrcMap.get(pos)[ptmId].mass) <= 0.02){
                        modPepPoolGood.push(peptide);
                    } else {
                        modPepPoolBad.push(peptide);
                    }
                }
            }
        }

        //find keptZone n1Mass and freeZone
        if (modPepPoolBad.peptideTreeSet.isEmpty()) {
//            System.out.println(scanNum);
            int a= 1;
            return new DividedZone(new HashSet<>(), new HashSet<>(), ptmMass);
        }
        Peptide inFeasiTopPep = modPepPoolBad.getTopPepPtn();
        Set<Integer> keptZone;
        Set<Integer> freeZone;
        int n1Pos = -1;
        for (Integer coor : inFeasiTopPep.getVarPTMs().keySet()) {
            if (modifiedZone.contains(coor)) {
                n1Pos = coor;
                ptmMass = inFeasiTopPep.getVarPTMs().get(coor);
            }
        }
        boolean yBetterB;
        if (inFeasiTopPep.matchedYions.size() > inFeasiTopPep.matchedBions.size()) {
            yBetterB = true;
        } else if (inFeasiTopPep.matchedYions.size() < inFeasiTopPep.matchedBions.size()) { /// B >>> Y
            yBetterB = false;
        } else {
            double bIntes = 0;
            double yIntes = 0;
            for (double mz : inFeasiTopPep.matchedBions.values()) {
                bIntes += mz;
            }
            for (double mz : inFeasiTopPep.matchedYions.values()) {
                yIntes += mz;
            }
            yBetterB = yIntes > bIntes;
        }
        if (yBetterB) { // Y >>> B
            int farestPos = Collections.max(modifiedZone);
            int n1LB = Collections.min(modifiedZone);
            int n1RB = Collections.max(modifiedZone);
            for (int pos : inFeasiTopPep.matchedYions.keySet()) {
                if (pos <= n1Pos) {
                    if (pos > n1LB) {
                        n1LB = pos;
                    }
                    if (pos < farestPos){
                        farestPos = pos;
                    }
                } else if (pos > n1Pos) {
                    if (pos < n1RB) {
                        n1RB = pos;
                    }
                } else {
                    System.out.println("lsz wrong 1");
                }
            }
            keptZone = IntStream.rangeClosed(n1LB, n1RB).boxed().collect(Collectors.toSet());
            freeZone = IntStream.range(Collections.min(modifiedZone), farestPos).boxed().collect(Collectors.toSet());
        } else { /// B >>> Y
            int closest = Collections.min(modifiedZone);
            int n1LB = Collections.min(modifiedZone);
            int n1RB = Collections.max(modifiedZone);
            for (int pos : inFeasiTopPep.matchedBions.keySet()) {
                if (pos < n1Pos) {
                    if (pos > n1LB) {
                        n1LB = pos+1;
                    }
                } else if (pos >= n1Pos) {
                    if (pos < n1RB) {
                        n1RB = pos;
                    }
                    if (pos > closest){
                        closest = pos;
                    }
                } else {
                    System.out.println("lsz wrong 1");
                }
            }
            keptZone = IntStream.rangeClosed(n1LB, n1RB).boxed().collect(Collectors.toSet());
            freeZone = IntStream.rangeClosed(closest+1, Collections.max(modifiedZone)).boxed().collect(Collectors.toSet());
        }
        return new DividedZone(keptZone, freeZone, ptmMass);
    }

    private DividedZone dividePepNewComple(int scanNum, Set<Integer> modifiedZone, ModPepPool modPepPoolGood, ModPepPool modPepPoolBad, Peptide lastPeptide, Map<Integer, VarPtm[]> posVarPtmArraySrcMap,
                                     double totalDeltaMass, String ptmFreePeptide, SparseVector expProcessedPL, TreeMap<Double, Double> plMap,Map<Integer, VarPtm> refVarPtmMap) throws CloneNotSupportedException {
        // Sometimes, the precursor mass error may affects the digitized spectrum.
        double ptmMass = 0;

        Set<VarPtm> refVarPtmSet = new HashSet<>(refVarPtmMap.values());
        for (int pos : modifiedZone) {
            if (!posVarPtmArraySrcMap.containsKey(pos)) continue;

            for (int ptmId = 0; ptmId < posVarPtmArraySrcMap.get(pos).length; ptmId++){
                Peptide peptide = lastPeptide.clone();
                PosMassMap newPosMassMap = new PosMassMap(ptmFreePeptide.length());
                for (Integer coor : lastPeptide.getVarPTMs().keySet()) {
                    newPosMassMap.put(coor, lastPeptide.getVarPTMs().get(coor)); // copy the ptms from lastPeptide
                }
                VarPtm thisPtm = posVarPtmArraySrcMap.get(pos)[ptmId];
                newPosMassMap.put(pos, thisPtm.mass);
                peptide.setVarPTM(newPosMassMap);
                int thisPriority = thisPtm.priority;


                double[][] ionMatrix = peptide.getIonMatrixNow();

                Set<Integer> jRange = IntStream.rangeClosed(0, ptmFreePeptide.length()-1).boxed().collect(Collectors.toSet());
                double score = massTool.buildVectorAndCalXCorr(ionMatrix, 1, expProcessedPL, peptide.matchedBions, peptide.matchedYions, jRange) ;
                double scoreLb = lastPeptide.getScore();

                if (Math.abs(totalDeltaMass - thisPtm.mass) <= 0.01){
                    scoreLb = lastPeptide.getScore()-1;
                }
                if (thisPriority == 1 ) { // if this ptm comes from the complementary peptide it should be priorized
                    scoreLb = -1;
                    score *= 2; //todo change it to that the priority one is only to consider not replacing the original one
                }
                if (refVarPtmSet.contains(thisPtm)) { // this is very important for synthetic dataset
                    score *= 1;
                }
                if (score > scoreLb) {
                    peptide.setScore(score);
                    peptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, peptide.getIonMatrix(), ms2Tolerance));
                    peptide.posVarPtmResMap.put(pos, thisPtm);
                    if (Math.abs(totalDeltaMass - thisPtm.mass) <= 0.02){
                        modPepPoolGood.push(peptide);
                    } else {
                        modPepPoolBad.push(peptide);
                    }
                }
            }
        }

        //find keptZone n1Mass and freeZone
        if (modPepPoolBad.peptideTreeSet.isEmpty()) {
//            System.out.println(scanNum);
            int a= 1;
            return new DividedZone(new HashSet<>(), new HashSet<>(), ptmMass);
        }
        Peptide inFeasiTopPep = modPepPoolBad.getTopPepPtn();
        Set<Integer> keptZone;
        Set<Integer> freeZone;
        int n1Pos = -1;
        for (Integer coor : inFeasiTopPep.getVarPTMs().keySet()) {
            if (modifiedZone.contains(coor)) {
                n1Pos = coor;
                ptmMass = inFeasiTopPep.getVarPTMs().get(coor);
            }
        }
        boolean yBetterB;
        if (inFeasiTopPep.matchedYions.size() > inFeasiTopPep.matchedBions.size()) {
            yBetterB = true;
        } else if (inFeasiTopPep.matchedYions.size() < inFeasiTopPep.matchedBions.size()) { /// B >>> Y
            yBetterB = false;
        } else {
            double bIntes = 0;
            double yIntes = 0;
            for (double mz : inFeasiTopPep.matchedBions.values()) {
                bIntes += mz;
            }
            for (double mz : inFeasiTopPep.matchedYions.values()) {
                yIntes += mz;
            }
            yBetterB = yIntes > bIntes;
        }
        if (yBetterB) { // Y >>> B
            int farestPos = Collections.max(modifiedZone);
            int n1LB = Collections.min(modifiedZone);
            int n1RB = Collections.max(modifiedZone);
            for (int pos : inFeasiTopPep.matchedYions.keySet()) {
                if (pos <= n1Pos) {
                    if (pos > n1LB) {
                        n1LB = pos;
                    }
                    if (pos < farestPos){
                        farestPos = pos;
                    }
                } else if (pos > n1Pos) {
                    if (pos < n1RB) {
                        n1RB = pos;
                    }
                } else {
                    System.out.println("lsz wrong 1");
                }
            }
            keptZone = IntStream.rangeClosed(n1LB, n1RB).boxed().collect(Collectors.toSet());
            freeZone = IntStream.range(Collections.min(modifiedZone), farestPos).boxed().collect(Collectors.toSet());
        } else { /// B >>> Y
            int closest = Collections.min(modifiedZone);
            int n1LB = Collections.min(modifiedZone);
            int n1RB = Collections.max(modifiedZone);
            for (int pos : inFeasiTopPep.matchedBions.keySet()) {
                if (pos < n1Pos) {
                    if (pos > n1LB) {
                        n1LB = pos+1;
                    }
                } else if (pos >= n1Pos) {
                    if (pos < n1RB) {
                        n1RB = pos;
                    }
                    if (pos > closest){
                        closest = pos;
                    }
                } else {
                    System.out.println("lsz wrong 1");
                }
            }
            keptZone = IntStream.rangeClosed(n1LB, n1RB).boxed().collect(Collectors.toSet());
            freeZone = IntStream.rangeClosed(closest+1, Collections.max(modifiedZone)).boxed().collect(Collectors.toSet());
        }
        return new DividedZone(keptZone, freeZone, ptmMass);
    }

    private DividedZone dividePep(int scanNum, Set<Integer> modifiedZone, ModPepPool ptmGoodRes, ModPepPool ptmBadRes, ModPepPool ptmTemplate, Map<Integer, VarPtm[]> idxVarModMap, double totalDeltaMass
            , String ptmFreePeptide, boolean isDecoy, SparseVector expProcessedPL, TreeMap<Double, Double> plMap, int precursorCharge
            , int localMaxMS2Charge, SortedMap<Double, Double> unUsedPlMap) { // Sometimes, the precursor mass error may affects the digitized spectrum.

        double ptmMass = 0;
        for (int pos : modifiedZone) {
            if (!idxVarModMap.containsKey(pos)) continue;

            for (int ptmId = 0; ptmId < idxVarModMap.get(pos).length; ptmId++){
                Peptide peptide = new Peptide(ptmFreePeptide, isDecoy, massTool);
                PosMassMap newPtmPtn = new PosMassMap(ptmFreePeptide.length());
                for (Integer coor : ptmTemplate.getTopPepPtn().getVarPTMs().keySet()) {
                    newPtmPtn.put(coor, ptmTemplate.getTopPepPtn().getVarPTMs().get(coor));
                }
                newPtmPtn.put(pos, idxVarModMap.get(pos)[ptmId].mass);
                peptide.setVarPTM(newPtmPtn);

                Map<Integer, Double> matchedBions = new HashMap<>();
                Map<Integer, Double> matchedYions = new HashMap<>();
                double[][] temp = peptide.getIonMatrixNow();
                if (ptmId == 0 && pos == 22){
                    int a = 1;
                }
                int lb = Math.max(0, Collections.min(modifiedZone));
                int rb = Math.min(ptmFreePeptide.length()-1, Collections.max(modifiedZone));
                Set<Integer> jRange = IntStream.rangeClosed(lb, rb).boxed().collect(Collectors.toSet());
                double score = massTool.buildVectorAndCalXCorr(peptide.getIonMatrixNow(), 1, expProcessedPL, matchedBions, matchedYions, jRange) ;
                double scoreLb = 0;
                if (Math.abs(totalDeltaMass - idxVarModMap.get(pos)[ptmId].mass) <= 0.01){
                    scoreLb = -1;
                }
                if (score > scoreLb) {
//                    System.out.println("numPtmsOnPep " + numPtmsOnPep  + " XCorr " + score);
                    peptide.setScore(score);
                    peptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, peptide.getIonMatrix(), ms2Tolerance));
                    peptide.matchedBions.putAll(matchedBions);
                    peptide.matchedYions.putAll(matchedYions);
                    if (Math.abs(totalDeltaMass - idxVarModMap.get(pos)[ptmId].mass) <= 0.01){
                        ptmGoodRes.push(peptide);
                    } else {
                        ptmBadRes.push(peptide);
                    }
                }
            }
        }

        //find keptZone n1Mass and freeZone
        if (ptmBadRes.peptideTreeSet.isEmpty()) {
//            System.out.println(scanNum);
            int a= 1;
            return new DividedZone(new HashSet<>(), new HashSet<>(), ptmMass);
        }
        Peptide inFeasiTopPep = ptmBadRes.getTopPepPtn();
        Set<Integer> keptZone = new HashSet<>();
        Set<Integer> freeZone = new HashSet<>();
        int n1Pos = -1;
        for (Integer coor : inFeasiTopPep.getVarPTMs().keySet()) {
            if (modifiedZone.contains(coor)) {
                n1Pos = coor;
                ptmMass = inFeasiTopPep.getVarPTMs().get(coor);
            }
        }
        boolean yBetterB = true;
        if (inFeasiTopPep.matchedYions.size() > inFeasiTopPep.matchedBions.size()) {
            yBetterB = true;
        } else if (inFeasiTopPep.matchedYions.size() < inFeasiTopPep.matchedBions.size()) { /// B >>> Y
            yBetterB = false;
        } else {
            double bIntes = 0;
            double yIntes = 0;
            for (double mz : inFeasiTopPep.matchedBions.values()) {
                bIntes += mz;
            }
            for (double mz : inFeasiTopPep.matchedYions.values()) {
                yIntes += mz;
            }
            yBetterB = yIntes > bIntes;
        }
        if (yBetterB) { // Y >>> B
            int farestPos = Collections.max(modifiedZone);
            int n1LB = Collections.min(modifiedZone);
            int n1RB = Collections.max(modifiedZone);
            for (int pos : inFeasiTopPep.matchedYions.keySet()) {
                if (pos <= n1Pos) {
                    if (pos > n1LB) {
                        n1LB = pos;
                    }
                    if (pos < farestPos){
                        farestPos = pos;
                    }
                } else if (pos > n1Pos) {
                    if (pos < n1RB) {
                        n1RB = pos;
                    }
                } else {
                    System.out.println("lsz wrong 1");
                }
            }
            keptZone = IntStream.rangeClosed(n1LB, n1RB).boxed().collect(Collectors.toSet());
            freeZone = IntStream.range(Collections.min(modifiedZone), farestPos).boxed().collect(Collectors.toSet());
        } else { /// B >>> Y
            int closest = Collections.min(modifiedZone);
            int n1LB = Collections.min(modifiedZone);
            int n1RB = Collections.max(modifiedZone);
            for (int pos : inFeasiTopPep.matchedBions.keySet()) {
                if (pos < n1Pos) {
                    if (pos > n1LB) {
                        n1LB = pos+1;
                    }
                } else if (pos >= n1Pos) {
                    if (pos < n1RB) {
                        n1RB = pos;
                    }
                    if (pos > closest){
                        closest = pos;
                    }
                } else {
                    System.out.println("lsz wrong 1");
                }
            }
            keptZone = IntStream.rangeClosed(n1LB, n1RB).boxed().collect(Collectors.toSet());
            freeZone = IntStream.rangeClosed(closest+1, Collections.max(modifiedZone)).boxed().collect(Collectors.toSet());
        }
        return new DividedZone(keptZone, freeZone, ptmMass);
    }

    public double getMinPtmMass() {
        return minPtmMass;
    }

    public double getMaxPtmMass() {
        return maxPtmMass;
    }

    private void readModFromUnimod() throws Exception {
        SAXReader reader = new SAXReader();
        InputStream inputStream = getClass().getClassLoader().getResourceAsStream("unimod.xml"); // PTMs from Unimod except for AA substitutions, isotopic labellings
//        BufferedReader reader1 = new BufferedReader(new InputStreamReader(inputStream));
////        inputStream.getClass().
        Document document = reader.read(inputStream);
        Element rootElement = document.getRootElement();
        Iterator<Element> rootIter = rootElement.elementIterator();

        while (rootIter.hasNext()) {
            Element rootElem = rootIter.next();
            if (!rootElem.getName().contentEquals("modifications")) continue;

            Iterator<Element> modIter = rootElem.elementIterator();

            while (modIter.hasNext()) {
                Element modElem = modIter.next();

                String name = modElem.attributeValue("title");
                if (name.contentEquals("Diethylphosphothione")) {
                    int a =1;
                }
                double mass = Double.valueOf(modElem.element("delta").attributeValue("mono_mass"));
                if (mass < minPtmMass || mass > maxPtmMass) continue;
                if (Math.abs(mass - 56.0626) <0.01) {
                    int a = 1;
                }
                for (Element spec : modElem.elements("specificity")) {
                    String classification = spec.attributeValue("classification");
                    if ( classification.contains("glycos") || classification.contains("Other")) {
                        continue;
                    }
//                    if (!recordIdWithPsiName.contains(recordId) && !classification.contentEquals("AA substitution")) continue;

                    String siteStr = spec.attributeValue("site");
                    String positionStr = spec.attributeValue("position");
                    if (classification.contentEquals("Isotopic label") && !(name.contentEquals("Propionyl") && siteStr.contentEquals("K")) && !(name.contentEquals("Succinyl") && siteStr.contentEquals("K"))) { // only for synthetic ptm data, because the authors uses them
                        continue;
                    }
//                    if (classification.contentEquals("Isotopic label") && !(name.contentEquals("Succinyl") && siteStr.contentEquals("K")) ) {
//                        continue;
//                    }
                    int position = 0;
                    switch (positionStr) {
                        case "Protein N-term":
                            position = 0;
                            break;
                        case "Protein C-term":
                            position = 1;
                            break;
                        case "Any N-term":
                            position = 2;
                            break;
                        case "Any C-term":
                            position = 3;
                            break;
                        case "Anywhere":
                            position = 4;
                            break;
                    }
                    if (siteStr.contentEquals("N-term") || siteStr.contentEquals("C-term")) {
                        for (char site : aaCharSet) {
                            if (aaWithFixModSet.contains(site) && siteStr.contentEquals("C-term")) {
                                continue; // if aa is C, just ignore the mod that are not at N term
                            }

                            VarPtm temp = new VarPtm(mass, site, position, name, classification, 0);
                            if (finalPtmMap.containsKey(site)) {
                                finalPtmMap.get(site).add(temp);
                            } else {
                                List<VarPtm> varPtmSet = new LinkedList<>();
                                varPtmSet.add(temp);
                                finalPtmMap.put(site, varPtmSet);
                            }
                        }
                    } else {
                        char site = siteStr.charAt(0);
                        if (aaWithFixModSet.contains(site) && (position == 1 || position == 3 || position == 4)) {
                            continue;  // if aa is C, just ignore the mod that are not at N term
                        }
                        VarPtm temp = new VarPtm(mass, site, position, name, classification, 0);
                        if (finalPtmMap.containsKey(site)) {
                            finalPtmMap.get(site).add(temp);
                        } else {
                            List<VarPtm> varPtmSet = new LinkedList<>();
                            varPtmSet.add(temp);
                            finalPtmMap.put(site, varPtmSet);
                        }
                    }
                }
            }
        }
    }

    private double calculateMassFromComposition(String composition) throws Exception {
        String[] parts = composition.split(" ");
        double mass = 0;
        for (String part : parts) {
            Matcher matcher = pattern.matcher(part.trim());
            if (matcher.matches()) {
                String element = matcher.group(1);
                int num = 1;
                if (matcher.group(2) != null) {
                    num = Integer.valueOf(matcher.group(3));
                }
                mass += num * elementTable.get(element);
            } else {
                throw new Exception(String.format(Locale.US, "The composition %s cannot be recognized.", part));
            }
        }
        return mass;
    }


    public Set<Integer> getFixModIdxes(String ptmFreePeptide) {
        Set<Integer> outputSet = new HashSet<>(ptmFreePeptide.length(), 1);
        char[] tempArray = ptmFreePeptide.toCharArray();
        for (int i = 0; i < tempArray.length; ++i) {
            if (Math.abs(fixModMap.get(tempArray[i])) > 0.1) {
                outputSet.add(i);
            }
        }
        return outputSet;
    }

    private Map<Integer, Set<VarPtm>> getIdxVarModMap(String ptmFreePeptide, Set<Integer> fixModIdxes, char leftFlank, char rightFlank) {
        Map<Integer, Set<VarPtm>> idxVarModMap = new HashMap<>(ptmFreePeptide.length() + 1, 1);
        for (int i = 0; i < ptmFreePeptide.length(); ++i) {
            if (!fixModIdxes.contains(i)) {
                char aa = ptmFreePeptide.charAt(i);
                if (aa == 'n') {
                    if (finalPtmMap.containsKey('n')) {
                        Set<VarPtm> tempSet = new HashSet<>();
                        for (VarPtm modEntry : finalPtmMap.get('n')) {
                            if (!modEntry.onlyProteinTerminalIfnc || leftFlank == '-') {
                                if ((leftFlank != 'K' || Math.abs(massTable.get('K') - modEntry.mass) > ms2Tolerance) && (leftFlank != 'R' || Math.abs(massTable.get('R') - modEntry.mass) > ms2Tolerance) && (massTable.get(ptmFreePeptide.charAt(1)) + modEntry.mass > ms2Tolerance)) {  // Fixing missed cleavages caused issue in N-term and the mass of a modified amino acid cannot be 0 or negative.
                                    tempSet.add(modEntry);
                                }
                            }
                        }
                        if (!tempSet.isEmpty()) {
                            idxVarModMap.put(0, tempSet);
                        }
                    }
                } else if (aa == 'c') {
                    if (finalPtmMap.containsKey('c')) {
                        Set<VarPtm> tempSet = new HashSet<>();
                        for (VarPtm modEntry : finalPtmMap.get('c')) {
                            if (!modEntry.onlyProteinTerminalIfnc || rightFlank == '-') {
                                if ((rightFlank != 'K' || Math.abs(massTable.get('K') - modEntry.mass) > ms2Tolerance) && (rightFlank != 'R' || Math.abs(massTable.get('R') - modEntry.mass) > ms2Tolerance) && (massTable.get(ptmFreePeptide.charAt(ptmFreePeptide.length() - 2)) + modEntry.mass > ms2Tolerance)) {  // Fixing missed cleavages caused issue in C-term and the mass of a modified amino acid cannot be 0 or negative
                                    tempSet.add(modEntry);
                                }
                            }
                        }
                        if (!tempSet.isEmpty()) {
                            idxVarModMap.put(ptmFreePeptide.length() - 1, tempSet);
                        }
                    }
                } else {
                    if (finalPtmMap.containsKey(aa)) {
                        idxVarModMap.put(i, new HashSet<>(finalPtmMap.get(aa)));
                    }
                }
            }
        }
        return idxVarModMap;
    }
    public Map<Double, Integer> getIdxVarModMapNew(String partSeq, Set<Integer> fixModIdxes, int isNorC_Part, int isProtNorC_Term) {
//        partSeq, fixModIdxes, isNorC_Side, isProtNorC_Term
        Map<Double, Integer> allMassMaxTimeMap = new HashMap<>(partSeq.length() * 10, 1);
        boolean hasProt_N_TermPtm = false;
        boolean hasProt_C_TermPtm = false;
        if (isNorC_Part == N_PART && isProtNorC_Term == N_TERM_PROT) {
            hasProt_N_TermPtm = true;
        }
        if (isNorC_Part == C_PART && isProtNorC_Term == C_TERM_PROT) {
            hasProt_C_TermPtm = true;
        }
        for (int i = 0; i < partSeq.length(); ++i) {
            if (fixModIdxes.contains(i) && i != 0) continue;  // if that pos has fix mod but is not N term, dont let it
            char aa = partSeq.charAt(i);

            if (finalPtmMap.containsKey(aa)) {
                Set<Double> massesOnAa = new HashSet<>();
                List<VarPtm> srcSet = finalPtmMap.get(aa);
                if (i == 0) { //aa at seq n term
                    for (VarPtm varPtm : srcSet) {
                        if (massTable.get(partSeq.charAt(i)) + varPtm.mass < ms2Tolerance) continue;//the mass of a modified amino acid cannot be 0 or negative

                        if (varPtm.position == 4 || (varPtm.position == 2 && isNorC_Part == N_PART) || (varPtm.position == 0 && hasProt_N_TermPtm)) { // anywhere or pepN or (protN and pepPos at protN)
//                            massesOnAa.add(Double.parseDouble(df2.format(varPtm.mass)));
                            massesOnAa.add(varPtm.mass);
                        }
                    }
                } else if (i == partSeq.length()-1) { //aa at seq c term
                    for (VarPtm varPtm : srcSet) {
                        if (massTable.get(partSeq.charAt(i)) + varPtm.mass < ms2Tolerance) continue;//the mass of a modified amino acid cannot be 0 or negative

                        if (varPtm.position == 4 || (varPtm.position == 3 && isNorC_Part == C_PART) || (varPtm.position == 1 && hasProt_C_TermPtm)) { // anywhere or pepC or (protC and pepPos at protC)
//                            massesOnAa.add(Double.parseDouble(df2.format(varPtm.mass)));
                            massesOnAa.add(varPtm.mass);
                        }
                    }
                } else {//aa at middle
                    for (VarPtm varPtm : srcSet) {
                        if (massTable.get(partSeq.charAt(i)) + varPtm.mass < ms2Tolerance) continue;//the mass of a modified amino acid cannot be 0 or negative

                        if (varPtm.position == 4) { // anywhere
//                            massesOnAa.add(Double.parseDouble(df2.format(varPtm.mass)));
                            massesOnAa.add(varPtm.mass);
                        }
                    }
                }

                for (double mass : massesOnAa) {
                    if (allMassMaxTimeMap.containsKey(mass)) {
                        allMassMaxTimeMap.put(mass, allMassMaxTimeMap.get(mass)+1);
                    } else {
                        allMassMaxTimeMap.put(mass, 1);
                    }
                }
            }
        }
        return allMassMaxTimeMap;
    }

    public Map<Integer, VarPtm[]> getIdxVarModMapNewComple(String partSeq, Set<Integer> fixModIdxes, int isProtNorC_Term) {
//        partSeq, fixModIdxes, isNorC_Side, isProtNorC_Term
        Map<Integer, VarPtm[]> idxVarModMap = new HashMap<>(partSeq.length() + 1, 1);
        boolean hasProt_N_TermPtm = false;
        boolean hasProt_C_TermPtm = false;
        if (isProtNorC_Term == N_TERM_PROT) {
            hasProt_N_TermPtm = true;
        }
        if (isProtNorC_Term == C_TERM_PROT) {
            hasProt_C_TermPtm = true;
        }
        for (int i = 0; i < partSeq.length(); ++i) {
            if (fixModIdxes.contains(i) && i != 0) continue;  // if that pos has fix mod but is not N term, dont let it
            char aa = partSeq.charAt(i);

            if (finalPtmMap.containsKey(aa)) {
                Map<String, VarPtm> dstMap = new HashMap<>();
                List<VarPtm> srcSet = finalPtmMap.get(aa);
                if (i == 0) { //aa at seq n term
                    for (VarPtm varPtm : srcSet) {
                        if (massTable.get(partSeq.charAt(i)) + varPtm.mass < ms2Tolerance) continue;//the mass of a modified amino acid cannot be 0 or negative

                        if (varPtm.position == 4 || (varPtm.position == 2) || (varPtm.position == 0 && hasProt_N_TermPtm)) { // anywhere or pepN or (protN and pepPos at protN)
                            String varPtmStr = varPtm.getStr();
                            if (dstMap.containsKey(varPtmStr)) {
                                VarPtm oldVarPtm = dstMap.get(varPtmStr);
                                if (varPtm.priority > oldVarPtm.priority){
                                    dstMap.put(varPtmStr, varPtm);
                                }
                            } else {
                                dstMap.put(varPtmStr, varPtm);
                            }
                        }
                    }
                } else if (i == partSeq.length()-1) { //aa at seq c term
                    for (VarPtm varPtm : srcSet) {
                        if (massTable.get(partSeq.charAt(i)) + varPtm.mass < ms2Tolerance) continue;//the mass of a modified amino acid cannot be 0 or negative

                        if (varPtm.position == 4 || (varPtm.position == 3) || (varPtm.position == 1 && hasProt_C_TermPtm)) { // anywhere or pepC or (protC and pepPos at protC)
                            String varPtmStr = varPtm.getStr();
                            if (dstMap.containsKey(varPtmStr)) {
                                VarPtm oldVarPtm = dstMap.get(varPtmStr);
                                if (varPtm.priority > oldVarPtm.priority){
                                    dstMap.put(varPtmStr, varPtm);
                                }
                            } else {
                                dstMap.put(varPtmStr, varPtm);
                            }
                        }
                    }
                } else {//aa at middle
                    for (VarPtm varPtm : srcSet) {
                        if (massTable.get(partSeq.charAt(i)) + varPtm.mass < ms2Tolerance) continue;//the mass of a modified amino acid cannot be 0 or negative

                        if (varPtm.position == 4) { // anywhere
                            String varPtmStr = varPtm.getStr();
                            if (dstMap.containsKey(varPtmStr)) {
                                VarPtm oldVarPtm = dstMap.get(varPtmStr);
                                if (varPtm.priority > oldVarPtm.priority){
                                    dstMap.put(varPtmStr, varPtm);
                                }
                            } else {
                                dstMap.put(varPtmStr, varPtm);
                            }
                        }
                    }
                }
                if (!dstMap.isEmpty()) {
                    VarPtm[] modArray = new VarPtm[dstMap.size()];
                    dstMap.values().toArray(modArray);
                    Arrays.sort(modArray, Comparator.comparingDouble(VarPtm::getMass));
                    idxVarModMap.put(i, modArray);
                }
            }
        }
        return idxVarModMap;
    }

    private Map<Integer, Set<VarPtm>> getIdxVarModMapNewOld(String ptmFreePeptide, Set<Integer> fixModIdxes, int pepPos, int tagPosInPep, int tagLen) {
        Map<Integer, Set<VarPtm>> idxVarModMap = new HashMap<>(ptmFreePeptide.length() + 1, 1);
        boolean nHasProtPtm = false;
        boolean cHasProtPtm = false;
        if (pepPos == -1) {
            nHasProtPtm = true;
        } else if (pepPos == 1) {
            cHasProtPtm = true;
        }
        if (pepPos != 0) {
            int a = 1;
        }
        for (int i = 0; i < ptmFreePeptide.length(); ++i) {
            if (i >= tagPosInPep && i < tagPosInPep+tagLen) continue; // dont consider var mod in tag regions
            if (fixModIdxes.contains(i) && i != 0) continue;  // if that pos has fix mod but is not N term, dont let it
            char aa = ptmFreePeptide.charAt(i);

            if (finalPtmMap.containsKey(aa)) {
                Map<String, VarPtm> dstMap = new HashMap<>();
                List<VarPtm> srcSet = finalPtmMap.get(aa);
                if (i == 0) { //aa at seq n term
                    for (VarPtm varPtm : srcSet) {
                        if (massTable.get(ptmFreePeptide.charAt(i)) + varPtm.mass < ms2Tolerance) continue;//the mass of a modified amino acid cannot be 0 or negative

                        if (varPtm.position == 4 || varPtm.position == 2 || (varPtm.position == 0 && nHasProtPtm)) { // anywhere or pepN or (protN and pepPos at protN)
                            String varPtmStr = varPtm.getStr();
                            if (dstMap.containsKey(varPtmStr)) {
                                VarPtm oldVarPtm = dstMap.get(varPtmStr);
                                if (varPtm.priority > oldVarPtm.priority){
                                    dstMap.put(varPtmStr, varPtm);
                                }
                            } else {
                                dstMap.put(varPtmStr, varPtm);
                            }
                        }
                    }
                } else if (i == ptmFreePeptide.length()-1) { //aa at seq c term
                    for (VarPtm varPtm : srcSet) {
                        if (massTable.get(ptmFreePeptide.charAt(i)) + varPtm.mass < ms2Tolerance) continue;//the mass of a modified amino acid cannot be 0 or negative

                        if (varPtm.position == 4 || varPtm.position == 3 || (varPtm.position == 1 && cHasProtPtm)) { // anywhere or pepC or (protC and pepPos at protC)
                            String varPtmStr = varPtm.getStr();
                            if (dstMap.containsKey(varPtmStr)) {
                                VarPtm oldVarPtm = dstMap.get(varPtmStr);
                                if (varPtm.priority > oldVarPtm.priority){
                                    dstMap.put(varPtmStr, varPtm);
                                }
                            } else {
                                dstMap.put(varPtmStr, varPtm);
                            }
                        }
                    }
                } else {//aa at middle
                    for (VarPtm varPtm : srcSet) {
                        if (massTable.get(ptmFreePeptide.charAt(i)) + varPtm.mass < ms2Tolerance) continue;//the mass of a modified amino acid cannot be 0 or negative

                        if (varPtm.position == 4) { // anywhere
                            String varPtmStr = varPtm.getStr();
                            if (dstMap.containsKey(varPtmStr)) {
                                VarPtm oldVarPtm = dstMap.get(varPtmStr);
                                if (varPtm.priority > oldVarPtm.priority){
                                    dstMap.put(varPtmStr, varPtm);
                                }
                            } else {
                                dstMap.put(varPtmStr, varPtm);
                            }
                        }
                    }
                }
                if (!dstMap.isEmpty()) {
                    idxVarModMap.put(i, new HashSet<>(dstMap.values()));
                }
            }
        }
        return idxVarModMap;
    }


    private class DividedZone {
        public Set<Integer> settledZone;
        public Set<Integer> toModZone;
        //        public double bMass = 0;
//        public double yMass = 0;
        public double ptmMass = 0;

        public DividedZone(Set<Integer> n1Zone, Set<Integer> toModZone, double  ptmMass) {
            this.settledZone = n1Zone;
            this.toModZone = toModZone;
//            this.bMass = bMass;
//            this.yMass = yMass;
            this.ptmMass = ptmMass;
        }
    }
}
