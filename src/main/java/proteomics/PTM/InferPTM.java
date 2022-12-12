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
import proteomics.Types.*;

import java.io.*;
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

    private final MassTool massTool;
    private final Map<String, Double> elementTable;
    private final Map<Character, Double> massTable;
    private final Map<Character, Double> fixModMap;
    private Set<VarPtm> varPtmSet = new HashSet<>();
    private final double minPtmMass;
    private final double maxPtmMass;
    private final double ms2Tolerance;
    private Map<Character, Set<VarPtm>> finalPtmMap = new HashMap<>();
    private final Set<Character> aaCharSet = new HashSet<>(Arrays.asList('A','C','D','E','F','G','H','I','K','L','M','N','O','P','Q','R','S','T','U','V','W','Y'));
    public InferPTM(MassTool massTool, Map<Character, Double> fixModMap, Map<String, String> parameterMap) throws Exception{
        this.massTool = massTool;
        elementTable = massTool.getElementTable();
        massTable = massTool.getMassTable();
        this.fixModMap = fixModMap;
        this.minPtmMass = Double.valueOf(parameterMap.get("min_ptm_mass"));
        this.maxPtmMass = Double.valueOf(parameterMap.get("max_ptm_mass"));
        this.ms2Tolerance = Double.valueOf(parameterMap.get("ms2_tolerance"));
        // Generate a varModParamSet from the parameterMap
        for (String k : parameterMap.keySet()) {
            if (k.startsWith("mod")) {
                String v = parameterMap.get(k);
                if (!v.startsWith("0.0")) {
                    String[] temp = v.split(",");
//                    System.out.println(temp[1].charAt(0));
                    if (Math.abs(fixModMap.get(temp[1].charAt(0))) < 0.1) {
                        // fix modification and var modification cannot be coexist
//                            public VarModParam(double mass, char site, int position,  String name, String classification, int priority) {
                        varPtmSet.add(new VarPtm(Double.valueOf(temp[0]), temp[1].charAt(0), Integer.valueOf(temp[2]), temp[3], "ByUser", 1)); // var mods from the parameter file have the highest priority, those PTM can exist in peptide terminal.
                    }  // todo varPtmSet should be transformed to 5 mod set
                }
            }
        }

        // Multimap<Character, ModEntry>
        // Reading Mod table...
        readModFromUnimod();

        // update ptm table with the high priority mods in the parameter file.
        for (VarPtm varPtm : varPtmSet) {
            if (finalPtmMap.containsKey(varPtm.site)) {
                if (finalPtmMap.get(varPtm.site).contains(varPtm)) {
                    finalPtmMap.get(varPtm.site).remove(varPtm); // remove old one add new one with priority
                }
                finalPtmMap.get(varPtm.site).add(varPtm);
            } else {
                Set<VarPtm> tempSet = new HashSet<>();
                tempSet.add(varPtm);
                finalPtmMap.put(varPtm.site, tempSet);
            }
        }

        // remove abundant records that are on the same site with lower level positions
        for (Character site : finalPtmMap.keySet()) {
            Set<VarPtm> toRemove = new HashSet<>();
            Iterator<VarPtm> iterI = finalPtmMap.get(site).iterator();
            while (iterI.hasNext()) {
                VarPtm ptmI = iterI.next();
                Iterator<VarPtm> iterJ = finalPtmMap.get(site).iterator();
                String nameI = ptmI.toString().substring(0, ptmI.toString().length()-1);
                int posI = Integer.valueOf(ptmI.toString().charAt(ptmI.toString().length()-1));
                while (iterJ.hasNext()) {
                    VarPtm ptmJ = iterJ.next();
                    String nameJ = ptmJ.toString().substring(0, ptmJ.toString().length()-1);
                    int posJ = Integer.valueOf(ptmJ.toString().charAt(ptmJ.toString().length()-1));
                    if (!ptmI.toString().contentEquals(ptmJ.toString()) && nameI.contentEquals(nameJ)) {
                        if (posI < posJ) {
                            toRemove.add(ptmI);
                        } else {
                            toRemove.add(ptmJ);
                        }
                    }
                }
            }
            finalPtmMap.get(site).removeAll(toRemove);
        }
        int a = 1;
    }

    public PeptidePTMPattern findPtmNew1(int scanNum, SparseVector expProcessedPL, TreeMap<Double, Double> plMap, double precursorMass, Peptide candiPep, PeptideInfo peptideInfo
                                        , int precursorCharge, int localMaxMS2Charge, double localMS1ToleranceL, double localMS1ToleranceR) {
//        double ptmFreeMass = massTool.calResidueMass(ptmFreePeptide) + massTool.H2O;
        // use the old way of findPtm but blocking tag region in the way as I block region with fixed mod
        String freeSeq = candiPep.getFreeSeq();
        boolean isDecoy = candiPep.isDecoy;
        double tagVecScore = candiPep.getTagVecScore();
        int globalRank = candiPep.getGlobalRank();
        char leftFlank = peptideInfo.leftFlank;
        char rightFlank = peptideInfo.rightFlank;

        double nDeltaMass = candiPep.nDeltaMass;
        double cDeltaMass = candiPep.cDeltaMass;
        int tagPosInPep = candiPep.tagPosInPep;
        int tagLen = 0;

        double score = 0;
        if (candiPep.finderTag != null) {
            score = candiPep.finderTag.getTotalIntensity();
            tagLen = candiPep.finderTag.size();
        }

        PeptidePTMPattern peptidePTMPattern = new PeptidePTMPattern(freeSeq);
        Peptide peptide = new Peptide(freeSeq, isDecoy, massTool, localMaxMS2Charge, tagVecScore, globalRank);

        double totalDeltaMass = precursorMass - peptide.getTheoMass();
        peptide.absDeltaMass = totalDeltaMass;
        if (Math.abs(totalDeltaMass) < 0.01 ) {
            return peptidePTMPattern;
        }
        Set<Integer> fixModIdxes = getFixModIdxes(freeSeq, fixModMap);  // positions that has fixed mod on it. Those postions should not bear var mod then.

//        String ptmFreePeptideOrdinary = freeSeq.replaceAll("I","L");

        int pepPosInProt = 0; // none of the terms is protein term
        if (leftFlank == '-') {
            pepPosInProt = -1; //nterm is protein N term
        } else if (rightFlank == '-') {
            pepPosInProt = 1;  //cterm is protein C term
        }

        Map<Integer, Set<VarPtm>> idxVarModMap = getIdxVarModMapNew(freeSeq, fixModIdxes, pepPosInProt, tagPosInPep, tagLen); //todo no need to generate var mod list for aa again and again, make it stored.
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


//        String truthPep = "n"+"DNVFENNRLAFEVAEK"+"c";
//        Peptide p1p = new Peptide(truthPep, isDecoy, massTool, localMaxMS2Charge, tagVecScore, globalRank);
//        PositionDeltaMassMap fakePatten = new PositionDeltaMassMap(freeSeq.length());
////        fakePatten.put(new Coordinate(0, 1), (58.005));
////        fakePatten.put(new Coordinate(22, 23), (33.988));
////        p1p.setVarPTM(fakePatten);
//        Map<Integer, Double> matchedBions1 = new HashMap<>();
//        Map<Integer, Double> matchedYions1 = new HashMap<>();
//        double[][] temp11 = p1p.getIonMatrixNow();
//        Set<Integer> jRange1 = IntStream.rangeClosed(0, truthPep.length()-3).boxed().collect(Collectors.toSet());
//        double ss = massTool.buildVectorAndCalXCorr(p1p.getIonMatrixNow(), 1, expProcessedPL, matchedBions1, matchedYions1, jRange1) ;
//        //stub end

        if (lszDebugScanNum.contains(scanNum)) {
            int a = 1;
        }
        if (modifiedZone.size() == 0) {
//            System.out.println(scanNum + " is empty modifiedZone after tag 2");
            return peptidePTMPattern; //Some scans are not valid Scans. Will be deleted soon.
        }

        TreeMap<Double, Double> unUsedPlMap = new TreeMap<>(plMap);

        PeptidePTMPattern allPtmPattern = new PeptidePTMPattern(freeSeq,1);
        PeptidePTMPattern allPtmPatternBad = new PeptidePTMPattern(freeSeq,1);
        modifiedZone = IntStream.rangeClosed(0, freeSeq.length()-1).boxed().collect(Collectors.toSet());// e.g. nABCDEFc, modifiedZone={1,2,3,4,5,6}

        Peptide cleanPep = new Peptide(freeSeq, isDecoy, massTool, localMaxMS2Charge, tagVecScore, globalRank);

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
//        modifiedZone = IntStream.range(1, ptmFreePeptide.length()-1).boxed().collect(Collectors.toSet());// e.g. nABCDEFc, modifiedZone={1,2,3,4,5,6}

        allPtmPatternBad.push(cleanPep);
        allPtmPattern.push(cleanPep);
        allPtmPatternBad.bestPep = allPtmPatternBad.getTopPepPtn();
        if (modifiedZone.isEmpty()) {
//            allPtmPatternBad.bestPep.
            return allPtmPatternBad;
        }


        PeptidePTMPattern ptmInitialTemp = new PeptidePTMPattern(freeSeq, 1);
        Peptide cleanPeptide = new Peptide(freeSeq, isDecoy, massTool, localMaxMS2Charge, tagVecScore, globalRank);
        cleanPeptide.setVarPTM(new PositionDeltaMassMap(freeSeq.length()));
        cleanPeptide.shouldPTM = true;
        ptmInitialTemp.push(cleanPeptide);
//        long t1 = System.currentTimeMillis();

        // 1st
        PeptidePTMPattern ptmN1Good = new PeptidePTMPattern(freeSeq, 1);
        PeptidePTMPattern ptmN1Bad = new PeptidePTMPattern(freeSeq, 1);

        DividedZone z1Z2Res = dividePep(scanNum, modifiedZone, ptmN1Good, ptmN1Bad, ptmInitialTemp, idxVarModArrayMap, totalDeltaMass, freeSeq, isDecoy, tagVecScore, globalRank, expProcessedPL, plMap, precursorCharge, localMaxMS2Charge, unUsedPlMap);
        if (!ptmN1Good.getPeptideTreeSet().isEmpty()) {
            ptmN1Good.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN1Good.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 1*0.1);
            ptmN1Good.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, ptmN1Good.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
            ptmN1Good.getTopPepPtn().shouldPTM = true;

            allPtmPattern.push(ptmN1Good.getTopPepPtn());
        }
        if (!ptmN1Bad.getPeptideTreeSet().isEmpty()) {
            ptmN1Bad.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN1Bad.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 1*0.1);
            ptmN1Bad.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, ptmN1Bad.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
            ptmN1Bad.getTopPepPtn().shouldPTM = true;
            allPtmPatternBad.push(ptmN1Bad.getTopPepPtn());
        }
        Set<Integer> zone1 = z1Z2Res.keptZone;
        Set<Integer> zone2 = z1Z2Res.freeZone;
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
        PeptidePTMPattern ptmN2Good = new PeptidePTMPattern(freeSeq, 1);
        PeptidePTMPattern ptmN2Bad = new PeptidePTMPattern(freeSeq, 1);

        DividedZone z3Z4Res = dividePep(scanNum, zone2, ptmN2Good, ptmN2Bad ,ptmN1Bad, idxVarModArrayMap, zone2Mass, freeSeq, isDecoy, tagVecScore, globalRank, expProcessedPL, plMap, precursorCharge, localMaxMS2Charge, unUsedPlMap);
        if (!ptmN2Good.getPeptideTreeSet().isEmpty()) {
            ptmN2Good.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN2Good.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 2*0.1);
            ptmN2Good.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, ptmN2Good.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
            ptmN2Good.getTopPepPtn().shouldPTM = true;

            allPtmPattern.push(ptmN2Good.getTopPepPtn());
        }
        if (!ptmN2Bad.getPeptideTreeSet().isEmpty()) {
            ptmN2Bad.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN2Bad.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 2*0.1);
            ptmN2Bad.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, ptmN2Bad.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
            ptmN2Bad.getTopPepPtn().shouldPTM = true;
            allPtmPatternBad.push(ptmN2Bad.getTopPepPtn());
        }
        Set<Integer> zone3 = z3Z4Res.keptZone;
        Set<Integer> zone4 = z3Z4Res.freeZone;
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
        PeptidePTMPattern ptmN3Good = new PeptidePTMPattern(freeSeq, 1);
        PeptidePTMPattern ptmN3Bad = new PeptidePTMPattern(freeSeq, 1);

        DividedZone z5Z6Res = dividePep(scanNum, zone4, ptmN3Good, ptmN3Bad, ptmN2Bad, idxVarModArrayMap, zone4Mass, freeSeq, isDecoy, tagVecScore, globalRank, expProcessedPL, plMap, precursorCharge, localMaxMS2Charge, unUsedPlMap);
        if (!ptmN3Good.peptideTreeSet.isEmpty()) {
            ptmN3Good.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN3Good.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 3*0.1);
            ptmN3Good.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, ptmN3Good.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
            ptmN3Good.getTopPepPtn().shouldPTM = true;

            allPtmPattern.push(ptmN3Good.getTopPepPtn());
        }
        if (!ptmN3Bad.getPeptideTreeSet().isEmpty()) {
            ptmN3Bad.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN3Bad.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 3*0.1);
            ptmN3Bad.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, ptmN3Bad.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
            ptmN3Bad.getTopPepPtn().shouldPTM = true;
            allPtmPatternBad.push(ptmN3Bad.getTopPepPtn());
        }
        Set<Integer> zone5 = z5Z6Res.keptZone;
        Set<Integer> zone6 = z5Z6Res.freeZone;
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

        PeptidePTMPattern ptmN4Good = new PeptidePTMPattern(freeSeq, 1);
        PeptidePTMPattern ptmN4Bad = new PeptidePTMPattern(freeSeq, 1);
        DividedZone z7Z8Res = dividePep(scanNum, zone6, ptmN4Good, ptmN4Bad ,ptmN3Bad, idxVarModArrayMap, zone6Mass, freeSeq, isDecoy, tagVecScore, globalRank, expProcessedPL, plMap, precursorCharge, localMaxMS2Charge, unUsedPlMap);
        if (!ptmN4Good.getPeptideTreeSet().isEmpty()) {
            ptmN4Good.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN4Good.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 4*0.1);
            ptmN4Good.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, ptmN4Good.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
            ptmN4Good.getTopPepPtn().shouldPTM = true;
            allPtmPattern.push(ptmN4Good.getTopPepPtn());
        }
        if (!ptmN4Bad.getPeptideTreeSet().isEmpty()) {
            ptmN4Bad.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN4Bad.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 4*0.1);
            ptmN4Bad.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, ptmN4Bad.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
            ptmN4Bad.getTopPepPtn().shouldPTM = true;

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

    public PeptidePTMPattern findPtmNew2(int scanNum, SparseVector expProcessedPL, TreeMap<Double, Double> plMap, double precursorMass, Peptide candiPep, PeptideInfo peptideInfo
            , int precursorCharge, int localMaxMS2Charge, double localMS1ToleranceL, double localMS1ToleranceR) {
//        double ptmFreeMass = massTool.calResidueMass(ptmFreePeptide) + massTool.H2O;
        String freeSeq = candiPep.getFreeSeq();
        boolean isDecoy = candiPep.isDecoy;
        double tagVecScore = candiPep.getTagVecScore();
        int globalRank = candiPep.getGlobalRank();
        char leftFlank = peptideInfo.leftFlank;
        char rightFlank = peptideInfo.rightFlank;

        double nDeltaMass = candiPep.nDeltaMass;
        double cDeltaMass = candiPep.cDeltaMass;
        int tagPosInPep = candiPep.tagPosInPep;

        double score = candiPep.finderTag.getTotalIntensity();

        PeptidePTMPattern peptidePTMPattern = new PeptidePTMPattern(freeSeq);
        Peptide peptide = new Peptide(freeSeq, isDecoy, massTool, localMaxMS2Charge, tagVecScore, globalRank);

        double totalDeltaMass = precursorMass - peptide.getTheoMass();
        peptide.absDeltaMass = totalDeltaMass;
        if (Math.abs(totalDeltaMass) < 0.01 ) {
            return peptidePTMPattern;
        }
        Set<Integer> fixModIdxes = getFixModIdxes(freeSeq, fixModMap);  // positions that has fixed mod on it. Those postions should not bear var mod then.

//        String ptmFreePeptideOrdinary = freeSeq.replaceAll("I","L");

        int pepPosInProt = 0; // none of the terms is protein term
        if (leftFlank == '-') {
            pepPosInProt = -1; //nterm is protein N term
        } else if (rightFlank == '-') {
            pepPosInProt = 1;  //cterm is protein C term
        }

        Map<Integer, Set<VarPtm>> idxVarModMap = getIdxVarModMapNew(freeSeq, fixModIdxes, pepPosInProt, tagPosInPep, candiPep.finderTag.size()); //todo no need to generate var mod list for aa again and again, make it stored.
        Map<Integer, VarPtm[]> idxVarModArrayMap = new HashMap<>();
        for (int id : idxVarModMap.keySet()){
            VarPtm[] modArray = new VarPtm[idxVarModMap.get(id).size()];
            idxVarModMap.get(id).toArray(modArray);
            Arrays.sort(modArray, Comparator.comparingDouble(VarPtm::getMass));
            idxVarModArrayMap.put(id, modArray);
        }

        if (nDeltaMass >= 0.2) { // nterm has unsettled mass and must has untagged amino acid for var mod

        }

        if (cDeltaMass >= 0.2) { // nterm has unsettled mass and must has untagged amino acid for var mod

        }
        Set<Integer> modifiedZone = new HashSet<>(idxVarModArrayMap.keySet());// e.g. nABCDEFc, modifiedZone={1,2,3,4,5,6}

        if (pepPosInProt != 0){
            int a = 1;
        }


//        String truthPep = "n"+"DNVFENNRLAFEVAEK"+"c";
//        Peptide p1p = new Peptide(truthPep, isDecoy, massTool, localMaxMS2Charge, tagVecScore, globalRank);
//        PositionDeltaMassMap fakePatten = new PositionDeltaMassMap(freeSeq.length());
////        fakePatten.put(new Coordinate(0, 1), (58.005));
////        fakePatten.put(new Coordinate(22, 23), (33.988));
////        p1p.setVarPTM(fakePatten);
//        Map<Integer, Double> matchedBions1 = new HashMap<>();
//        Map<Integer, Double> matchedYions1 = new HashMap<>();
//        double[][] temp11 = p1p.getIonMatrixNow();
//        Set<Integer> jRange1 = IntStream.rangeClosed(0, truthPep.length()-3).boxed().collect(Collectors.toSet());
//        double ss = massTool.buildVectorAndCalXCorr(p1p.getIonMatrixNow(), 1, expProcessedPL, matchedBions1, matchedYions1, jRange1) ;
//        //stub end

        if (lszDebugScanNum.contains(scanNum)) {
            int a = 1;
        }
        if (modifiedZone.size() == 0) {
//            System.out.println(scanNum + " is empty modifiedZone after tag 2");
            return peptidePTMPattern; //Some scans are not valid Scans. Will be deleted soon.
        }

        TreeMap<Double, Double> unUsedPlMap = new TreeMap<>(plMap);

        PeptidePTMPattern allPtmPattern = new PeptidePTMPattern(freeSeq,1);
        PeptidePTMPattern allPtmPatternBad = new PeptidePTMPattern(freeSeq,1);
        modifiedZone = IntStream.rangeClosed(0, freeSeq.length()-1).boxed().collect(Collectors.toSet());// e.g. nABCDEFc, modifiedZone={1,2,3,4,5,6}

        Peptide cleanPep = new Peptide(freeSeq, isDecoy, massTool, localMaxMS2Charge, tagVecScore, globalRank);

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
//        modifiedZone = IntStream.range(1, ptmFreePeptide.length()-1).boxed().collect(Collectors.toSet());// e.g. nABCDEFc, modifiedZone={1,2,3,4,5,6}

        allPtmPatternBad.push(cleanPep);
        allPtmPattern.push(cleanPep);
        allPtmPatternBad.bestPep = allPtmPatternBad.getTopPepPtn();
        if (modifiedZone.isEmpty()) {
            return allPtmPatternBad;
        }


        PeptidePTMPattern ptmInitialTemp = new PeptidePTMPattern(freeSeq, 1);
        Peptide cleanPeptide = new Peptide(freeSeq, isDecoy, massTool, localMaxMS2Charge, tagVecScore, globalRank);
        cleanPeptide.setVarPTM(new PositionDeltaMassMap(freeSeq.length()));
        cleanPeptide.shouldPTM = true;
        ptmInitialTemp.push(cleanPeptide);
//        long t1 = System.currentTimeMillis();

        // 1st
        PeptidePTMPattern ptmN1Good = new PeptidePTMPattern(freeSeq, 1);
        PeptidePTMPattern ptmN1Bad = new PeptidePTMPattern(freeSeq, 1);

        DividedZone z1Z2Res = dividePep(scanNum, modifiedZone, ptmN1Good, ptmN1Bad, ptmInitialTemp, idxVarModArrayMap, totalDeltaMass, freeSeq, isDecoy, tagVecScore, globalRank, expProcessedPL, plMap, precursorCharge, localMaxMS2Charge, unUsedPlMap);
        if (!ptmN1Good.getPeptideTreeSet().isEmpty()) {
            ptmN1Good.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN1Good.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 1*0.1);
            ptmN1Good.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, ptmN1Good.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
            ptmN1Good.getTopPepPtn().shouldPTM = true;

            allPtmPattern.push(ptmN1Good.getTopPepPtn());
        }
        if (!ptmN1Bad.getPeptideTreeSet().isEmpty()) {
            ptmN1Bad.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN1Bad.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 1*0.1);
            ptmN1Bad.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, ptmN1Bad.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
            ptmN1Bad.getTopPepPtn().shouldPTM = true;
            allPtmPatternBad.push(ptmN1Bad.getTopPepPtn());
        }
        Set<Integer> zone1 = z1Z2Res.keptZone;
        Set<Integer> zone2 = z1Z2Res.freeZone;
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
        PeptidePTMPattern ptmN2Good = new PeptidePTMPattern(freeSeq, 1);
        PeptidePTMPattern ptmN2Bad = new PeptidePTMPattern(freeSeq, 1);

        DividedZone z3Z4Res = dividePep(scanNum, zone2, ptmN2Good, ptmN2Bad ,ptmN1Bad, idxVarModArrayMap, zone2Mass, freeSeq, isDecoy, tagVecScore, globalRank, expProcessedPL, plMap, precursorCharge, localMaxMS2Charge, unUsedPlMap);
        if (!ptmN2Good.getPeptideTreeSet().isEmpty()) {
            ptmN2Good.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN2Good.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 2*0.1);
            ptmN2Good.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, ptmN2Good.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
            ptmN2Good.getTopPepPtn().shouldPTM = true;

            allPtmPattern.push(ptmN2Good.getTopPepPtn());
        }
        if (!ptmN2Bad.getPeptideTreeSet().isEmpty()) {
            ptmN2Bad.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN2Bad.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 2*0.1);
            ptmN2Bad.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, ptmN2Bad.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
            ptmN2Bad.getTopPepPtn().shouldPTM = true;
            allPtmPatternBad.push(ptmN2Bad.getTopPepPtn());
        }
        Set<Integer> zone3 = z3Z4Res.keptZone;
        Set<Integer> zone4 = z3Z4Res.freeZone;
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
        PeptidePTMPattern ptmN3Good = new PeptidePTMPattern(freeSeq, 1);
        PeptidePTMPattern ptmN3Bad = new PeptidePTMPattern(freeSeq, 1);

        DividedZone z5Z6Res = dividePep(scanNum, zone4, ptmN3Good, ptmN3Bad, ptmN2Bad, idxVarModArrayMap, zone4Mass, freeSeq, isDecoy, tagVecScore, globalRank, expProcessedPL, plMap, precursorCharge, localMaxMS2Charge, unUsedPlMap);
        if (!ptmN3Good.peptideTreeSet.isEmpty()) {
            ptmN3Good.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN3Good.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 3*0.1);
            ptmN3Good.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, ptmN3Good.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
            ptmN3Good.getTopPepPtn().shouldPTM = true;

            allPtmPattern.push(ptmN3Good.getTopPepPtn());
        }
        if (!ptmN3Bad.getPeptideTreeSet().isEmpty()) {
            ptmN3Bad.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN3Bad.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 3*0.1);
            ptmN3Bad.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, ptmN3Bad.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
            ptmN3Bad.getTopPepPtn().shouldPTM = true;
            allPtmPatternBad.push(ptmN3Bad.getTopPepPtn());
        }
        Set<Integer> zone5 = z5Z6Res.keptZone;
        Set<Integer> zone6 = z5Z6Res.freeZone;
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

        PeptidePTMPattern ptmN4Good = new PeptidePTMPattern(freeSeq, 1);
        PeptidePTMPattern ptmN4Bad = new PeptidePTMPattern(freeSeq, 1);
        DividedZone z7Z8Res = dividePep(scanNum, zone6, ptmN4Good, ptmN4Bad ,ptmN3Bad, idxVarModArrayMap, zone6Mass, freeSeq, isDecoy, tagVecScore, globalRank, expProcessedPL, plMap, precursorCharge, localMaxMS2Charge, unUsedPlMap);
        if (!ptmN4Good.getPeptideTreeSet().isEmpty()) {
            ptmN4Good.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN4Good.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 4*0.1);
            ptmN4Good.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, ptmN4Good.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
            ptmN4Good.getTopPepPtn().shouldPTM = true;
            allPtmPattern.push(ptmN4Good.getTopPepPtn());
        }
        if (!ptmN4Bad.getPeptideTreeSet().isEmpty()) {
            ptmN4Bad.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN4Bad.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 4*0.1);
            ptmN4Bad.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, ptmN4Bad.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
            ptmN4Bad.getTopPepPtn().shouldPTM = true;

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

//    public PeptidePTMPattern findPTM(int scanNum, SparseVector expProcessedPL, TreeMap<Double, Double> plMap, double precursorMass, String ptmFreePeptide, boolean isDecoy,
//                                     double normalizedCrossCorr, char leftFlank, char rightFlank, int globalRank, int precursorCharge, int localMaxMS2Charge,
//                                     double localMS1ToleranceL, double localMS1ToleranceR) {
////        double ptmFreeMass = massTool.calResidueMass(ptmFreePeptide) + massTool.H2O;
//        PeptidePTMPattern peptidePTMPattern = new PeptidePTMPattern(ptmFreePeptide);
//        Peptide peptide = new Peptide(ptmFreePeptide, isDecoy, massTool, localMaxMS2Charge, normalizedCrossCorr, globalRank);
//        double totalDeltaMass = precursorMass - peptide.getTheoMass();
//        peptide.absDeltaMass = totalDeltaMass;
//        if (Math.abs(totalDeltaMass) < 0.01 ) {
//            return peptidePTMPattern;
//        }
//        Set<Integer> fixModIdxes = getFixModIdxes(ptmFreePeptide, fixModMap);  // record the indexes of the amino acids with fix mod, e.g. when there is no C, this is empty set.
//
//        String ptmFreePeptideOrdinary = ptmFreePeptide.replaceAll("I","L");
//
//        int pepPos = 0;
//        if (leftFlank == '-') {
//            pepPos = -1;
//        } else if (rightFlank == '-') {
//            pepPos = 1;
//        }
//
//        Map<Integer, Set<VarPtm>> idxVarModMap = getIdxVarModMapNew(ptmFreePeptide, fixModIdxes, pepPos);
//        Map<Integer, VarPtm[]> idxVarModArrayMap = new HashMap<>();
//        for (int id : idxVarModMap.keySet()){
//            VarPtm[] modArray = new VarPtm[idxVarModMap.get(id).size()];
//            idxVarModMap.get(id).toArray(modArray);
//            Arrays.sort(modArray, Comparator.comparingDouble(VarPtm::getMass));
//            idxVarModArrayMap.put(id, modArray);
//        }
//        Set<Integer> modifiedZone = new HashSet<>(idxVarModArrayMap.keySet());// e.g. nABCDEFc, modifiedZone={1,2,3,4,5,6}
//
//        if (pepPos != 0){
//            int a = 1;
//        }
//        String truthPep = "n"+"DNVFENNRLAFEVAEK"+"c";
//        Peptide p1p = new Peptide(truthPep, isDecoy, massTool, localMaxMS2Charge, normalizedCrossCorr, globalRank);
//        PositionDeltaMassMap fakePatten = new PositionDeltaMassMap(ptmFreePeptide.length());
////        fakePatten.put(new Coordinate(0, 1), (58.005));
////        fakePatten.put(new Coordinate(22, 23), (33.988));
////        p1p.setVarPTM(fakePatten);
//        Map<Integer, Double> matchedBions1 = new HashMap<>();
//        Map<Integer, Double> matchedYions1 = new HashMap<>();
//        double[][] temp11 = p1p.getIonMatrixNow();
//        Set<Integer> jRange1 = IntStream.rangeClosed(0, truthPep.length()-3).boxed().collect(Collectors.toSet());
//        double ss = massTool.buildVectorAndCalXCorr(p1p.getIonMatrixNow(), 1, expProcessedPL, matchedBions1, matchedYions1, jRange1) ;
//        //stub end
//
//        if (lszDebugScanNum.contains(scanNum)) {
//            int a = 1;
//        }
//        if (modifiedZone.size() == 0) {
////            System.out.println(scanNum + " is empty modifiedZone after tag 2");
//            return peptidePTMPattern; //Some scans are not valid Scans. Will be deleted soon.
//        }
//
//        TreeMap<Double, Double> unUsedPlMap = new TreeMap<>(plMap);
//
//        PeptidePTMPattern allPtmPattern = new PeptidePTMPattern(ptmFreePeptide,1);
//        PeptidePTMPattern allPtmPatternBad = new PeptidePTMPattern(ptmFreePeptide,1);
//        modifiedZone = IntStream.rangeClosed(0, ptmFreePeptide.length()-1).boxed().collect(Collectors.toSet());// e.g. nABCDEFc, modifiedZone={1,2,3,4,5,6}
//
//        Peptide cleanPep = new Peptide(ptmFreePeptide, isDecoy, massTool, localMaxMS2Charge, normalizedCrossCorr, globalRank);
//
//        int lb = 0;  //lb included
//        int rb = ptmFreePeptide.length() - 1;//rb included
//        Map<Integer, Double> matchedBions = new HashMap<>();
//        Map<Integer, Double> matchedYions = new HashMap<>();
//        double[][] temp1 = cleanPep.getIonMatrixNow();
//        Set<Integer> jRange = IntStream.rangeClosed(0, ptmFreePeptide.length()-1).boxed().collect(Collectors.toSet());
//        double[][] cleanPepIonMatrix = cleanPep.getIonMatrixNow();
//        double cleanScore = massTool.buildVectorAndCalXCorr(cleanPepIonMatrix, 1, expProcessedPL, matchedBions, matchedYions, jRange) ;
//        cleanPep.setScore(cleanScore);
//        cleanPep.setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, cleanPep.getIonMatrix(), ms2Tolerance));
//        cleanPep.matchedBions.putAll(matchedBions);
//        cleanPep.matchedYions.putAll(matchedYions);
//
//        List<Integer> matchedBIndex = new ArrayList<>(matchedBions.keySet());
//        List<Integer> matchedYIndex = new ArrayList<>(matchedYions.keySet());
//        Collections.sort(matchedBIndex);
//        Collections.sort(matchedYIndex);
//        if (matchedBIndex.size() >= 2) {
//            if (matchedBIndex.get(matchedBIndex.size()-1) + (1-matchedBions.get(matchedBIndex.get(matchedBIndex.size()-1))) > matchedBIndex.get(matchedBIndex.size()-2)+4) {
//                matchedBions.remove(matchedBIndex.get(matchedBIndex.size()-1));
//            }
//        }
//        if (matchedYIndex.size() >= 2) {
//            if (matchedYIndex.get(0) - (1-matchedYions.get(matchedYIndex.get(0))) < matchedYIndex.get(1)-4) {
//                matchedYions.remove(matchedYIndex.get(0));
//            }
//        }
//
//        if (!matchedBions.isEmpty()) {
//            lb = Collections.max(matchedBions.keySet()) + 1;
//        }
//        if (!matchedYions.isEmpty()) {
//            rb = Collections.min(matchedYions.keySet()) - 1;
//        }
//
//        if (rb - lb <= 0) {
//            double bSumIntens = 0;
//            for (double intes : matchedBions.values()) bSumIntens += intes;
//            double ySumIntens = 0;
//            for (double intes : matchedYions.values()) ySumIntens += intes;
//            if (bSumIntens > ySumIntens) {
//                rb = ptmFreePeptide.length() - 1;
//            } else {
//                lb = 0;
//            }
//        }
//        modifiedZone = IntStream.rangeClosed(lb, rb).boxed().collect(Collectors.toSet());
////        modifiedZone = IntStream.range(1, ptmFreePeptide.length()-1).boxed().collect(Collectors.toSet());// e.g. nABCDEFc, modifiedZone={1,2,3,4,5,6}
//
//        allPtmPatternBad.push(cleanPep);
//        allPtmPattern.push(cleanPep);
//        allPtmPatternBad.bestPep = allPtmPatternBad.getTopPepPtn();
//        if (modifiedZone.isEmpty()) {
//            return allPtmPatternBad;
//        }
//
//
//        PeptidePTMPattern ptmInitialTemp = new PeptidePTMPattern(ptmFreePeptide, 1);
//        Peptide cleanPeptide = new Peptide(ptmFreePeptide, isDecoy, massTool, localMaxMS2Charge, normalizedCrossCorr, globalRank);
//        cleanPeptide.setVarPTM(new PositionDeltaMassMap(ptmFreePeptide.length()));
//        cleanPeptide.shouldPTM = true;
//        ptmInitialTemp.push(cleanPeptide);
////        long t1 = System.currentTimeMillis();
//
//        // 1st
//        PeptidePTMPattern ptmN1Good = new PeptidePTMPattern(ptmFreePeptide, 1);
//        PeptidePTMPattern ptmN1Bad = new PeptidePTMPattern(ptmFreePeptide, 1);
//
//        DividedZone z1Z2Res = dividePep(scanNum, modifiedZone, ptmN1Good, ptmN1Bad, ptmInitialTemp, idxVarModArrayMap, totalDeltaMass, ptmFreePeptide, isDecoy, normalizedCrossCorr, globalRank, expProcessedPL, plMap, precursorCharge, localMaxMS2Charge, unUsedPlMap);
//        if (!ptmN1Good.getPeptideTreeSet().isEmpty()) {
//            ptmN1Good.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN1Good.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 1*0.1);
//            ptmN1Good.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, ptmN1Good.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
//            ptmN1Good.getTopPepPtn().shouldPTM = true;
//
//            allPtmPattern.push(ptmN1Good.getTopPepPtn());
//        }
//        if (!ptmN1Bad.getPeptideTreeSet().isEmpty()) {
//            ptmN1Bad.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN1Bad.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 1*0.1);
//            ptmN1Bad.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, ptmN1Bad.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
//            ptmN1Bad.getTopPepPtn().shouldPTM = true;
//            allPtmPatternBad.push(ptmN1Bad.getTopPepPtn());
//        }
//        Set<Integer> zone1 = z1Z2Res.keptZone;
//        Set<Integer> zone2 = z1Z2Res.freeZone;
//        if (ptmN1Bad.peptideTreeSet.isEmpty() || zone2.isEmpty()) {
//            if (allPtmPattern.peptideTreeSet.isEmpty()) {
//                return allPtmPatternBad;
//            }
//            allPtmPattern.bestPep = allPtmPatternBad.getTopPepPtn();
//            return allPtmPattern;
//        }
//        Map<Integer, Double> matchedBionsBad1 = new HashMap<>();
//        Map<Integer, Double> matchedYionsBad1 = new HashMap<>();
////        ptmN1Bad.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN1Bad.getTopPepPtn().getIonMatrixNow(), 1,  expProcessedPL, matchedBionsBad1, matchedYionsBad1, IntStream.range(0, ptmFreePeptide.length()-2).boxed().collect(Collectors.toSet())));
//
////        ptmN1Bad.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN1Bad.getTopPepPtn().getIonMatrixNow(), 1,  expProcessedPL) - 1*0.1);
////        ptmN1Bad.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, localMaxMS2Charge, ptmN1Bad.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
////        allPtmPatternBad.push(ptmN1Bad.getTopPepPtn());
//
//        double zone1Mass = ptmN1Bad.getTopPepPtn().getVarPTMs().values().iterator().next();
//        double zone2Mass = totalDeltaMass - zone1Mass;
//
//        // 2nd
//        PeptidePTMPattern ptmN2Good = new PeptidePTMPattern(ptmFreePeptide, 1);
//        PeptidePTMPattern ptmN2Bad = new PeptidePTMPattern(ptmFreePeptide, 1);
//
//        DividedZone z3Z4Res = dividePep(scanNum, zone2, ptmN2Good, ptmN2Bad ,ptmN1Bad, idxVarModArrayMap, zone2Mass, ptmFreePeptide, isDecoy, normalizedCrossCorr, globalRank, expProcessedPL, plMap, precursorCharge, localMaxMS2Charge, unUsedPlMap);
//        if (!ptmN2Good.getPeptideTreeSet().isEmpty()) {
//            ptmN2Good.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN2Good.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 2*0.1);
//            ptmN2Good.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, ptmN2Good.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
//            ptmN2Good.getTopPepPtn().shouldPTM = true;
//
//            allPtmPattern.push(ptmN2Good.getTopPepPtn());
//        }
//        if (!ptmN2Bad.getPeptideTreeSet().isEmpty()) {
//            ptmN2Bad.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN2Bad.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 2*0.1);
//            ptmN2Bad.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, ptmN2Bad.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
//            ptmN2Bad.getTopPepPtn().shouldPTM = true;
//            allPtmPatternBad.push(ptmN2Bad.getTopPepPtn());
//        }
//        Set<Integer> zone3 = z3Z4Res.keptZone;
//        Set<Integer> zone4 = z3Z4Res.freeZone;
//        if (ptmN2Bad.peptideTreeSet.isEmpty() || zone4.isEmpty()) {
//            if (allPtmPattern.peptideTreeSet.isEmpty()) {
//                return allPtmPatternBad;
//            }
//            allPtmPattern.bestPep = allPtmPatternBad.getTopPepPtn();
//            return allPtmPattern;
//        }
////        ptmN2Bad.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN2Bad.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 2*0.1);
////        ptmN2Bad.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, localMaxMS2Charge, ptmN2Bad.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
////        allPtmPatternBad.push(ptmN2Bad.getTopPepPtn());
//
//        double zone3Mass = z3Z4Res.ptmMass;
//        double zone4Mass = zone2Mass - zone3Mass;
//
//        // 3rd
//        PeptidePTMPattern ptmN3Good = new PeptidePTMPattern(ptmFreePeptide, 1);
//        PeptidePTMPattern ptmN3Bad = new PeptidePTMPattern(ptmFreePeptide, 1);
//
//        DividedZone z5Z6Res = dividePep(scanNum, zone4, ptmN3Good, ptmN3Bad, ptmN2Bad, idxVarModArrayMap, zone4Mass, ptmFreePeptide, isDecoy, normalizedCrossCorr, globalRank, expProcessedPL, plMap, precursorCharge, localMaxMS2Charge, unUsedPlMap);
//        if (!ptmN3Good.peptideTreeSet.isEmpty()) {
//            ptmN3Good.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN3Good.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 3*0.1);
//            ptmN3Good.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, ptmN3Good.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
//            ptmN3Good.getTopPepPtn().shouldPTM = true;
//
//            allPtmPattern.push(ptmN3Good.getTopPepPtn());
//        }
//        if (!ptmN3Bad.getPeptideTreeSet().isEmpty()) {
//            ptmN3Bad.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN3Bad.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 3*0.1);
//            ptmN3Bad.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, ptmN3Bad.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
//            ptmN3Bad.getTopPepPtn().shouldPTM = true;
//            allPtmPatternBad.push(ptmN3Bad.getTopPepPtn());
//        }
//        Set<Integer> zone5 = z5Z6Res.keptZone;
//        Set<Integer> zone6 = z5Z6Res.freeZone;
//        if (ptmN3Bad.peptideTreeSet.isEmpty() || zone6.isEmpty()) {
//            if (allPtmPattern.peptideTreeSet.isEmpty()) {
//                return allPtmPatternBad;
//            }
//            allPtmPattern.bestPep = allPtmPatternBad.getTopPepPtn();
//            return allPtmPattern;
//        }
////        ptmN3Bad.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN3Bad.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 3*0.1);
////        ptmN3Bad.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, localMaxMS2Charge, ptmN3Bad.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
////        allPtmPatternBad.push(ptmN3Bad.getTopPepPtn());
//
//        double zone5Mass = z5Z6Res.ptmMass;
//        double zone6Mass = zone4Mass - zone5Mass;
//
//        PeptidePTMPattern ptmN4Good = new PeptidePTMPattern(ptmFreePeptide, 1);
//        PeptidePTMPattern ptmN4Bad = new PeptidePTMPattern(ptmFreePeptide, 1);
//        DividedZone z7Z8Res = dividePep(scanNum, zone6, ptmN4Good, ptmN4Bad ,ptmN3Bad, idxVarModArrayMap, zone6Mass, ptmFreePeptide, isDecoy, normalizedCrossCorr, globalRank, expProcessedPL, plMap, precursorCharge, localMaxMS2Charge, unUsedPlMap);
//        if (!ptmN4Good.getPeptideTreeSet().isEmpty()) {
//            ptmN4Good.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN4Good.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 4*0.1);
//            ptmN4Good.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, ptmN4Good.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
//            ptmN4Good.getTopPepPtn().shouldPTM = true;
//            allPtmPattern.push(ptmN4Good.getTopPepPtn());
//        }
//        if (!ptmN4Bad.getPeptideTreeSet().isEmpty()) {
//            ptmN4Bad.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN4Bad.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 4*0.1);
//            ptmN4Bad.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, ptmN4Bad.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
//            ptmN4Bad.getTopPepPtn().shouldPTM = true;
//
//            allPtmPatternBad.push(ptmN4Bad.getTopPepPtn());
//        }
//
//
//
//        long end = System.currentTimeMillis();
//
////        System.out.println("N1 time, "+ (end - start));
//        if (allPtmPattern.peptideTreeSet.isEmpty()) {
//            return allPtmPatternBad;
//        }
//        allPtmPattern.bestPep = allPtmPatternBad.getTopPepPtn();
//        return allPtmPattern;
//    }

    private DividedZone dividePep(int scanNum, Set<Integer> modifiedZone, PeptidePTMPattern ptmGoodRes, PeptidePTMPattern ptmBadRes, PeptidePTMPattern ptmTemplate, Map<Integer, VarPtm[]> idxVarModMap, double totalDeltaMass, String ptmFreePeptide, boolean isDecoy, double normalizedCrossCorr, int globalRank, SparseVector expProcessedPL, TreeMap<Double, Double> plMap, int precursorCharge, int localMaxMS2Charge, SortedMap<Double, Double> unUsedPlMap) { // Sometimes, the precursor mass error may affects the digitized spectrum.
        double ptmMass = 0;
        for (int pos : modifiedZone) {
            if (!idxVarModMap.containsKey(pos)) continue;

            for (int ptmId = 0; ptmId < idxVarModMap.get(pos).length; ptmId++){
                Peptide peptide = new Peptide(ptmFreePeptide, isDecoy, massTool, localMaxMS2Charge, normalizedCrossCorr, globalRank);
                PositionDeltaMassMap newPtmPtn = new PositionDeltaMassMap(ptmFreePeptide.length());
                for (Coordinate coor : ptmTemplate.getTopPepPtn().getVarPTMs().keySet()) {
                    newPtmPtn.put(new Coordinate(coor.x, coor.y), ptmTemplate.getTopPepPtn().getVarPTMs().get(coor));
                }
                newPtmPtn.put(new Coordinate(pos, pos + 1), idxVarModMap.get(pos)[ptmId].mass);
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
//                Set<Integer> jRange = new HashSet<>();
//                for (int i : modifiedZone) {
//                    jRange.add(i-1);
//                }
                double score = massTool.buildVectorAndCalXCorr(peptide.getIonMatrixNow(), 1, expProcessedPL, matchedBions, matchedYions, jRange) ;
//                double score = massTool.buildVectorAndCalXCorr(peptide.getIonMatrixNow(), 1, expProcessedPL, matchedBions, matchedYions,
//                        IntStream.range(Collections.min(modifiedZone)-1, ptmFreePeptide.length()-2).boxed().collect(Collectors.toSet()),
//                        IntStream.range(0, Collections.max(modifiedZone)).boxed().collect(Collectors.toSet())) ;
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
        for (Coordinate coor : inFeasiTopPep.getVarPTMs().keySet()) {
            if (modifiedZone.contains(coor.x)) {
                n1Pos = coor.x;
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
        Document document = reader.read(new File("/home/slaiad/Code/PIPI/src/main/resources/unimod.xml"));
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
                Set<Integer> recordIdWithPsiName = new HashSet<> (Arrays.asList(1,2,3,4,5,6,7,8,9,10,11,12,13,17,20,21,23,24,25,26,27,28,29,30,31,34,35,36,37,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,58,59,60,61,62,63,64,65,66,89,90,91,92,93,94,95,97,105,106,107,108,118,119,121,122,123,124,126,127,128,129,130,131,134,135,136,137,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,170,171,172,176,178,184,185,186,187,188,193,194,195,196,197,198,199,200,205,206,207,208,209,211,212,213,214,243,253,254,255,256,258,259,260,261,262,264,267,268,269,270,271,272,275,276,278,280,281,284,285,286,288,289,290,291,292,293,294,295,298,299,301,302,303,305,307,308,309,310,311,312,313,314,316,318,319,320,323,324,325,327,329,330,332,333,335,337,340,342,344,345,348,349,350,351,352,354,359,360,362,363,364,365,366,368,369,371,372,374,375,376,377,378,379,380,381,382,385,387,388,389,390,391,392,393,394,395,396,397,398,400,401,402,403,405,407,408,409,410,411,412,413,414,415,416,417,419,420,421,422,423,424,425,426,428,429,431,432,433,434,435,436,437,438,439,440,442,443,444,445,447,448,449,450,451,452,453,454,455,457,464,472,476,477,478,481,488,490,493,494,495,498,499,500,501,510,512,513,514,515,518,519,520,522,523,526,528,529,530,531,532,533,534,535,687,695,696,772));
                int recordId = Integer.valueOf(modElem.attributeValue("record_id"));
                double mass = Double.valueOf(modElem.element("delta").attributeValue("mono_mass"));
                if (mass < minPtmMass || mass > maxPtmMass) continue;
                if (Math.abs(mass - 56.0626) <0.01) {
                    int a = 1;
                }
                for (Element spec : modElem.elements("specificity")) {
                    String classification = spec.attributeValue("classification");
                    if (classification.contentEquals("Isotopic label") || classification.contains("glycos") || classification.contains("Other")) continue;

//                    if (!recordIdWithPsiName.contains(recordId) && !classification.contentEquals("AA substitution")) continue;

                    String siteStr = spec.attributeValue("site");
                    String positionStr = spec.attributeValue("position");
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
                            VarPtm temp = new VarPtm(mass, site, position, name, classification, 0);
                            if (finalPtmMap.containsKey(site)) {
                                finalPtmMap.get(site).add(temp);
                            } else {
                                Set<VarPtm> varPtmSet = new HashSet<>();
                                varPtmSet.add(temp);
                                finalPtmMap.put(site, varPtmSet);
                            }
                        }
                    } else {
                        char site = siteStr.charAt(0);
                        VarPtm temp = new VarPtm(mass, site, position, name, classification, 0);
                        if (finalPtmMap.containsKey(site)) {
                            finalPtmMap.get(site).add(temp);
                        } else {
                            Set<VarPtm> varPtmSet = new HashSet<>();
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


    private Set<Integer> getFixModIdxes(String ptmFreePeptide, Map<Character, Double> fixModMap) {
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

    private Map<Integer, Set<VarPtm>> getIdxVarModMapNew(String ptmFreePeptide, Set<Integer> fixModIdxes, int pepPos, int tagPosInPep, int tagLen) {
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
            if (fixModIdxes.contains(i)) continue;
            char aa = ptmFreePeptide.charAt(i);

            if (finalPtmMap.containsKey(aa)) {
                Set<VarPtm> dstSet = new HashSet<>();
                Set<VarPtm> srcSet = finalPtmMap.get(aa);
                if (i == 0) { //aa at seq n term
                    for (VarPtm varPtm : srcSet) {
                        if (massTable.get(ptmFreePeptide.charAt(i)) + varPtm.mass < ms2Tolerance) continue;//the mass of a modified amino acid cannot be 0 or negative

                        if (varPtm.position == 4 || varPtm.position == 2 || (varPtm.position == 0 && nHasProtPtm)) { // anywhere or pepN or (protN and pepPos at protN)
                            dstSet.add(varPtm);
                        }
                    }
                } else if (i == ptmFreePeptide.length()-1) { //aa at seq c term
                    for (VarPtm varPtm : srcSet) {
                        if (massTable.get(ptmFreePeptide.charAt(i)) + varPtm.mass < ms2Tolerance) continue;//the mass of a modified amino acid cannot be 0 or negative

                        if (varPtm.position == 4 || varPtm.position == 3 || (varPtm.position == 1 && cHasProtPtm)) { // anywhere or pepC or (protC and pepPos at protC)
                            dstSet.add(varPtm);
                        }
                    }
                } else {//aa at middle
                    for (VarPtm varPtm : srcSet) {
                        if (massTable.get(ptmFreePeptide.charAt(i)) + varPtm.mass < ms2Tolerance) continue;//the mass of a modified amino acid cannot be 0 or negative

                        if (varPtm.position == 4) { // anywhere
                            dstSet.add(varPtm);
                        }
                    }
                }
                if (!dstSet.isEmpty()) {
                    idxVarModMap.put(i, dstSet);
                }
            }
        }
        return idxVarModMap;
    }

    private class DividedZone {
        public Set<Integer> keptZone;
        public Set<Integer> freeZone;
        //        public double bMass = 0;
//        public double yMass = 0;
        public double ptmMass = 0;

        public DividedZone(Set<Integer> n1Zone, Set<Integer> freeZone, double  ptmMass) {
            this.keptZone = n1Zone;
            this.freeZone = freeZone;
//            this.bMass = bMass;
//            this.yMass = yMass;
            this.ptmMass = ptmMass;
        }
    }
}
