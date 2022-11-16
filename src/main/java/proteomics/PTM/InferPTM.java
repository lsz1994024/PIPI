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
import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
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

import org.dom4j.Attribute;
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
    private Set<VarModParam> varModParamSet = new HashSet<>();
    private final double minPtmMass;
    private final double maxPtmMass;
    private final double ms2Tolerance;
    private final Multimap<Character, ModEntry> siteModMap;

    private Map<Character, Set<VarModParam>> finalPtmMap = new HashMap<>();

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
                    String[] temp = v.split("@");
                    if (Math.abs(fixModMap.get(temp[1].charAt(0))) < 0.1) {
                        // fix modification and var modification cannot be coexist
                        varModParamSet.add(new VarModParam(Double.valueOf(temp[0]), temp[1].charAt(0), 1, false)); // var mods from the parameter file have the highest priority, those PTM can exist in peptide terminal.
                    }
                }
            } else if (k.contentEquals("Nterm")) {
                if (Math.abs(fixModMap.get('n')) < 0.1) {
                    // fix modification and var modification cannot be coexist
                    if (!parameterMap.get(k).startsWith("0.0")) {
                        String[] tempArray = parameterMap.get(k).split(",");
                        for (String tempString : tempArray) {
                            varModParamSet.add(new VarModParam(Double.valueOf(tempString.trim()), 'n', 1, false)); // var mods from the parameter file have the highest priority
                        }
                    }
                }
            } else if (k.contentEquals("Cterm")) {
                if (Math.abs(fixModMap.get('c')) < 0.1) {
                    // fix modification and var modification cannot be coexist
                    if (!parameterMap.get(k).startsWith("0.0")) {
                        String[] tempArray = parameterMap.get(k).split(",");
                        for (String tempString : tempArray) {
                            varModParamSet.add(new VarModParam(Double.valueOf(tempString.trim()), 'c', 1, false)); // var mods from the parameter file have the highest priority
                        }
                    }
                }
            }
        }

        // Generate a site mod map including all PTMs in Unimod and amino acid substitutions.
//        siteModMap = readUnimodAndGenerateAAS(minPtmMass, maxPtmMass);
        siteModMap = HashMultimap.create(); //to prove this siteModMap is useless and the unimod.xml.tsv

        // Building an amino acid substitution matrix.
        Map<Character, Set<VarModParam>> aasMap = buildAASMap(siteModMap);

        // Reading Mod table...
        Map<Character, Set<VarModParam>> ptmMap = readModXml();
        int numPtm = 0;
        for (Set<VarModParam> varSet : ptmMap.values()) {
            numPtm += varSet.size();
        }

        // update ptm table with the high priority mods in the parameter file.
        for (VarModParam varModParam : varModParamSet) {
            varModParam.position = "Anywhere";
            varModParam.classification = "User define";
            varModParam.name = "user";
            if (ptmMap.containsKey(varModParam.site)) {
//                ptmMap.get(varModParam.site).remove(varModParam);
                ptmMap.get(varModParam.site).add(varModParam);
            } else {
                Set<VarModParam> tempSet = new HashSet<>();
                tempSet.add(varModParam);
                ptmMap.put(varModParam.site, tempSet);
            }
        }

        // merge PTM map and AA substitution map
//        for (char aa : aasMap.keySet()) {
//            for (VarModParam modEntry : aasMap.get(aa)) {
//                if (finalPtmMap.containsKey(aa)) {
//                    if (finalPtmMap.get(aa).contains(modEntry)) {
//                        if (modEntry.priority == 1 || !modEntry.onlyProteinTerminalIfnc) {
//                            finalPtmMap.get(aa).remove(modEntry); // VarModParam only differs by mass and site.
//                            finalPtmMap.get(aa).add(modEntry); // temp is a high priority PTM, it is safe to overwrite the original PTM.
//                        }
//                    } else {
//                        finalPtmMap.get(aa).add(modEntry);
//                    }
//                } else {
//                    Set<VarModParam> tempSet = new HashSet<>();
//                    tempSet.add(modEntry);
//                    finalPtmMap.put(aa, tempSet);
//                }
//            }
//        }
        for (char aa : ptmMap.keySet()) {
            for (VarModParam modEntry : ptmMap.get(aa)) {
                if (finalPtmMap.containsKey(aa)) {
                    if (finalPtmMap.get(aa).contains(modEntry)) {
                        if (modEntry.priority == 1 || !modEntry.onlyProteinTerminalIfnc) {
                            finalPtmMap.get(aa).remove(modEntry); // VarModParam only differs by mass and site.
                            finalPtmMap.get(aa).add(modEntry); // temp is a high priority PTM, it is safe to overwrite the original PTM.
                        }
                    } else {
                        finalPtmMap.get(aa).add(modEntry);
                    }
                } else {
                    Set<VarModParam> tempSet = new HashSet<>();
                    tempSet.add(modEntry);
                    finalPtmMap.put(aa, tempSet);
                }
            }
        }
        int a = 1;
    }



    public PeptidePTMPattern findPTM(int scanNum, SparseVector expProcessedPL, TreeMap<Double, Double> plMap, double precursorMass, String ptmFreePeptide, boolean isDecoy,
                                     double normalizedCrossCorr, char leftFlank, char rightFlank, int globalRank, int precursorCharge, int localMaxMS2Charge,
                                     double localMS1ToleranceL, double localMS1ToleranceR) {
//        double ptmFreeMass = massTool.calResidueMass(ptmFreePeptide) + massTool.H2O;
        PeptidePTMPattern peptidePTMPattern = new PeptidePTMPattern(ptmFreePeptide);
        Peptide peptide = new Peptide(ptmFreePeptide, isDecoy, massTool, localMaxMS2Charge, normalizedCrossCorr, globalRank);
        double totalDeltaMass = precursorMass - peptide.getTheoMass();
        peptide.absDeltaMass = totalDeltaMass;
        if (Math.abs(totalDeltaMass) < 0.01 ) {
            return peptidePTMPattern;
        }
        Set<Integer> fixModIdxes = getFixModIdxes(ptmFreePeptide, fixModMap);  // record the indexes of the amino acids with fix mod, e.g. when there is no C, this is empty set.

        String ptmFreePeptideOrdinary = ptmFreePeptide.replaceAll("I","L");

        Map<Integer, Set<VarModParam>> idxVarModMap = getIdxVarModMapNew(ptmFreePeptide, fixModIdxes, leftFlank, rightFlank);
        Map<Integer, VarModParam[]> idxVarModArrayMap = new HashMap<>();
        for (int id : idxVarModMap.keySet()){
            VarModParam[] modArray = new VarModParam[idxVarModMap.get(id).size()];
            idxVarModMap.get(id).toArray(modArray);
            Arrays.sort(modArray, Comparator.comparingDouble(VarModParam::getMass));
            idxVarModArrayMap.put(id, modArray);
        }
        Set<Integer> modifiedZone = new HashSet<>(idxVarModArrayMap.keySet());// e.g. nABCDEFc, modifiedZone={1,2,3,4,5,6}


        String truthPep = "n"+"DNVFENNRLAFEVAEK"+"c";
        Peptide p1p = new Peptide(truthPep, isDecoy, massTool, localMaxMS2Charge, normalizedCrossCorr, globalRank);
        PositionDeltaMassMap fakePatten = new PositionDeltaMassMap(ptmFreePeptide.length());
//        fakePatten.put(new Coordinate(0, 1), (58.005));
//        fakePatten.put(new Coordinate(22, 23), (33.988));
//        p1p.setVarPTM(fakePatten);
        Map<Integer, Double> matchedBions1 = new HashMap<>();
        Map<Integer, Double> matchedYions1 = new HashMap<>();
        double[][] temp11 = p1p.getIonMatrixNow();
        Set<Integer> jRange1 = IntStream.rangeClosed(0, truthPep.length()-3).boxed().collect(Collectors.toSet());
        double ss = massTool.buildVectorAndCalXCorr(p1p.getIonMatrixNow(), 1, expProcessedPL, matchedBions1, matchedYions1, jRange1) ;
        //stub end

        if (lszDebugScanNum.contains(scanNum)) {
            int a = 1;
        }
        if (modifiedZone.size() == 0) {
//            System.out.println(scanNum + " is empty modifiedZone after tag 2");
            return peptidePTMPattern; //Some scans are not valid Scans. Will be deleted soon.
        }

        TreeMap<Double, Double> unUsedPlMap = new TreeMap<>(plMap);

        PeptidePTMPattern allPtmPattern = new PeptidePTMPattern(ptmFreePeptide,1);
        PeptidePTMPattern allPtmPatternBad = new PeptidePTMPattern(ptmFreePeptide,1);
        modifiedZone = IntStream.rangeClosed(0, ptmFreePeptide.length()-1).boxed().collect(Collectors.toSet());// e.g. nABCDEFc, modifiedZone={1,2,3,4,5,6}

        Peptide cleanPep = new Peptide(ptmFreePeptide, isDecoy, massTool, localMaxMS2Charge, normalizedCrossCorr, globalRank);

        int lb = 0;  //lb included
        int rb = ptmFreePeptide.length() - 1;//rb included
        Map<Integer, Double> matchedBions = new HashMap<>();
        Map<Integer, Double> matchedYions = new HashMap<>();
        double[][] temp1 = cleanPep.getIonMatrixNow();
        Set<Integer> jRange = IntStream.rangeClosed(0, ptmFreePeptide.length()-3).boxed().collect(Collectors.toSet());
        double[][] cleanPepIonMatrix = cleanPep.getIonMatrixNow();
        double cleanScore = massTool.buildVectorAndCalXCorr(cleanPepIonMatrix, 1, expProcessedPL, matchedBions, matchedYions, jRange) ;
        cleanPep.setScore(cleanScore);
        cleanPep.setMatchedPeakNum(Score.getMatchedIonNum(plMap, localMaxMS2Charge, cleanPep.getIonMatrix(), ms2Tolerance));
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
                rb = ptmFreePeptide.length() - 1;
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


        PeptidePTMPattern ptmInitialTemp = new PeptidePTMPattern(ptmFreePeptide, 1);
        Peptide cleanPeptide = new Peptide(ptmFreePeptide, isDecoy, massTool, localMaxMS2Charge, normalizedCrossCorr, globalRank);
        cleanPeptide.setVarPTM(new PositionDeltaMassMap(ptmFreePeptide.length()));
        cleanPeptide.shouldPTM = true;
        ptmInitialTemp.push(cleanPeptide);
//        long t1 = System.currentTimeMillis();

        // 1st
        PeptidePTMPattern ptmN1Good = new PeptidePTMPattern(ptmFreePeptide, 1);
        PeptidePTMPattern ptmN1Bad = new PeptidePTMPattern(ptmFreePeptide, 1);

        DividedZone z1Z2Res = dividePep(scanNum, modifiedZone, ptmN1Good, ptmN1Bad, ptmInitialTemp, idxVarModArrayMap, totalDeltaMass, ptmFreePeptide, isDecoy, normalizedCrossCorr, globalRank, expProcessedPL, plMap, precursorCharge, localMaxMS2Charge, unUsedPlMap);
        if (!ptmN1Good.getPeptideTreeSet().isEmpty()) {
            ptmN1Good.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN1Good.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 1*0.1);
            ptmN1Good.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, localMaxMS2Charge, ptmN1Good.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
            ptmN1Good.getTopPepPtn().shouldPTM = true;

            allPtmPattern.push(ptmN1Good.getTopPepPtn());
        }
        if (!ptmN1Bad.getPeptideTreeSet().isEmpty()) {
            ptmN1Bad.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN1Bad.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 1*0.1);
            ptmN1Bad.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, localMaxMS2Charge, ptmN1Bad.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
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
        PeptidePTMPattern ptmN2Good = new PeptidePTMPattern(ptmFreePeptide, 1);
        PeptidePTMPattern ptmN2Bad = new PeptidePTMPattern(ptmFreePeptide, 1);

        DividedZone z3Z4Res = dividePep(scanNum, zone2, ptmN2Good, ptmN2Bad ,ptmN1Bad, idxVarModArrayMap, zone2Mass, ptmFreePeptide, isDecoy, normalizedCrossCorr, globalRank, expProcessedPL, plMap, precursorCharge, localMaxMS2Charge, unUsedPlMap);
        if (!ptmN2Good.getPeptideTreeSet().isEmpty()) {
            ptmN2Good.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN2Good.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 2*0.1);
            ptmN2Good.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, localMaxMS2Charge, ptmN2Good.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
            ptmN2Good.getTopPepPtn().shouldPTM = true;

            allPtmPattern.push(ptmN2Good.getTopPepPtn());
        }
        if (!ptmN2Bad.getPeptideTreeSet().isEmpty()) {
            ptmN2Bad.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN2Bad.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 2*0.1);
            ptmN2Bad.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, localMaxMS2Charge, ptmN2Bad.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
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
        PeptidePTMPattern ptmN3Good = new PeptidePTMPattern(ptmFreePeptide, 1);
        PeptidePTMPattern ptmN3Bad = new PeptidePTMPattern(ptmFreePeptide, 1);

        DividedZone z5Z6Res = dividePep(scanNum, zone4, ptmN3Good, ptmN3Bad, ptmN2Bad, idxVarModArrayMap, zone4Mass, ptmFreePeptide, isDecoy, normalizedCrossCorr, globalRank, expProcessedPL, plMap, precursorCharge, localMaxMS2Charge, unUsedPlMap);
        if (!ptmN3Good.peptideTreeSet.isEmpty()) {
            ptmN3Good.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN3Good.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 3*0.1);
            ptmN3Good.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, localMaxMS2Charge, ptmN3Good.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
            ptmN3Good.getTopPepPtn().shouldPTM = true;

            allPtmPattern.push(ptmN3Good.getTopPepPtn());
        }
        if (!ptmN3Bad.getPeptideTreeSet().isEmpty()) {
            ptmN3Bad.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN3Bad.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 3*0.1);
            ptmN3Bad.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, localMaxMS2Charge, ptmN3Bad.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
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

        PeptidePTMPattern ptmN4Good = new PeptidePTMPattern(ptmFreePeptide, 1);
        PeptidePTMPattern ptmN4Bad = new PeptidePTMPattern(ptmFreePeptide, 1);
        DividedZone z7Z8Res = dividePep(scanNum, zone6, ptmN4Good, ptmN4Bad ,ptmN3Bad, idxVarModArrayMap, zone6Mass, ptmFreePeptide, isDecoy, normalizedCrossCorr, globalRank, expProcessedPL, plMap, precursorCharge, localMaxMS2Charge, unUsedPlMap);
        if (!ptmN4Good.getPeptideTreeSet().isEmpty()) {
            ptmN4Good.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN4Good.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 4*0.1);
            ptmN4Good.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, localMaxMS2Charge, ptmN4Good.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
            ptmN4Good.getTopPepPtn().shouldPTM = true;
            allPtmPattern.push(ptmN4Good.getTopPepPtn());
        }
        if (!ptmN4Bad.getPeptideTreeSet().isEmpty()) {
            ptmN4Bad.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN4Bad.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 4*0.1);
            ptmN4Bad.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, localMaxMS2Charge, ptmN4Bad.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
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

    private DividedZone dividePep(int scanNum, Set<Integer> modifiedZone, PeptidePTMPattern ptmGoodRes, PeptidePTMPattern ptmBadRes, PeptidePTMPattern ptmTemplate, Map<Integer, VarModParam[]> idxVarModMap, double totalDeltaMass, String ptmFreePeptide, boolean isDecoy, double normalizedCrossCorr, int globalRank, SparseVector expProcessedPL, TreeMap<Double, Double> plMap, int precursorCharge, int localMaxMS2Charge, SortedMap<Double, Double> unUsedPlMap) { // Sometimes, the precursor mass error may affects the digitized spectrum.
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
                int lb = Math.max(0, Collections.min(modifiedZone)-1);
                int rb = Math.min(ptmFreePeptide.length()-3, Collections.max(modifiedZone)-1);
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
                    peptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, localMaxMS2Charge, peptide.getIonMatrix(), ms2Tolerance));
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

//    public static Multimap<Character, ModEntry> readUnimodAndGenerateAAS(double minPtmMass, double maxPtmMass) throws IOException {
//        Multimap<Character, ModEntry> siteModMap = HashMultimap.create();
//        InputStream inputStream = OutputPeff.class.getClassLoader().getResourceAsStream("unimod.xml.tsv");
//        BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));
//        String line;
//        while ((line = reader.readLine()) != null) {
//            line = line.trim();
//            if (!line.isEmpty() && !line.startsWith("accession")) {
//                String[] parts = line.split("\t");
//                if (!parts[6].contentEquals("Isotopic label") && !parts[6].contentEquals("AA substitution")) { // We don't consider isotopic label and amino acid substitution..
//                    char site = parts[1].charAt(0);
//                    if (parts[1].trim().contentEquals("N-term")) {
//                        site = 'n';
//                    } else if (parts[1].trim().contentEquals("C-term")) {
//                        site = 'c';
//                    }
//                    if (site != 'n' && site != 'c') { // We don't consider peptide terminal modifications
//                        double ptmDeltaMass =  Double.valueOf(parts[4].trim());
//                        if (ptmDeltaMass >= minPtmMass && ptmDeltaMass <= maxPtmMass) { // only record the PTM within the delta mass range
//                            siteModMap.put(site, new ModEntry("UNIMOD:" + parts[0].trim(), parts[2].trim().contentEquals("null") ? parts[3].trim() : parts[2].trim(), ptmDeltaMass)); // if there is no PSI-MS name, use the internal name in the Unimod
//                        }
//                    }
//                }
//            }
//        }
//        reader.close();
//
//        return siteModMap;
//    }

    public Multimap<Character, ModEntry> getSiteModMap() {
        return siteModMap;
    }

    public double getMinPtmMass() {
        return minPtmMass;
    }

    public double getMaxPtmMass() {
        return maxPtmMass;
    }

    private Map<Character, Set<VarModParam>> readModFile() throws Exception {
        Map<Character, Set<VarModParam>> siteModMap = new HashMap<>();
        InputStream inputStream = getClass().getClassLoader().getResourceAsStream("modTable.tsv"); // PTMs from Unimod except for AA substitutions, isotopic labellings
        BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));
        String line;
        while ((line = reader.readLine()) != null) {
            line = line.trim();
            if (!line.isEmpty() && !line.startsWith("site")) {
                String[] parts = line.split("\t");
                if (!parts[1].trim().contentEquals("null")) { // We don't consider those without PSI-MS names.
                    char site = parts[0].charAt(0);
                    if (parts[0].trim().contentEquals("N-term")) {
                        site = 'n';
                    } else if (parts[0].trim().contentEquals("C-term")) {
                        site = 'c';
                    }
                    double mass = calculateMassFromComposition(parts[3].trim());
                    if (mass >= minPtmMass && mass <= maxPtmMass) {
                        if (site == 'n' || site == 'c' || massTable.get(site) + mass > ms2Tolerance) { // The mass of a modified amino acid cannot be 0 or negative.
                            int priority = Integer.valueOf(parts[4]);
                            VarModParam temp = new VarModParam(mass, site, priority, true); // natural N/C-terminal PTMs only happens on protein terminal.
                            if (siteModMap.containsKey(site)) {
                                if (siteModMap.get(site).contains(temp)) {
                                    if (temp.priority == 1) { // temp is a high priority PTM, it is safe to overwrite the original PTM.
                                        siteModMap.get(site).add(temp);
                                    }
                                } else {
                                    siteModMap.get(site).add(temp);
                                }
                            } else {
                                Set<VarModParam> tempSet = new HashSet<>();
                                tempSet.add(temp);
                                siteModMap.put(site, tempSet);
                            }
                        }
                    }
                }
            }
        }
        reader.close();
        return siteModMap;
    }

    private Map<Character, Set<VarModParam>> readModXml() throws Exception {
        Map<Character, Set<VarModParam>> siteModMap = new HashMap<>();

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
                Set<Integer> recordIdWithPsiName = new HashSet<> (Arrays.asList(1,2,3,4,5,6,7,8,9,10,11,12,13,17,20,21,23,24,25,26,27,28,29,30,31,34,35,36,37,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,58,59,60,61,62,63,64,65,66,89,90,91,92,93,94,95,97,105,106,107,108,118,119,121,122,123,124,126,127,128,129,130,131,134,135,136,137,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,170,171,172,176,178,184,185,186,187,188,193,194,195,196,197,198,199,200,205,206,207,208,209,211,212,213,214,243,253,254,255,256,258,259,260,261,262,264,267,268,269,270,271,272,275,276,278,280,281,284,285,286,288,289,290,291,292,293,294,295,298,299,301,302,303,305,307,308,309,310,311,312,313,314,316,318,319,320,323,324,325,327,329,330,332,333,335,337,340,342,344,345,348,349,350,351,352,354,359,360,362,363,364,365,366,368,369,371,372,374,375,376,377,378,379,380,381,382,385,387,388,389,390,391,392,393,394,395,396,397,398,400,401,402,403,405,407,408,409,410,411,412,413,414,415,416,417,419,420,421,422,423,424,425,426,428,429,431,432,433,434,435,436,437,438,439,440,442,443,444,445,447,448,449,450,451,452,453,454,455,457,464,472,476,477,478,481,488,490,493,494,495,498,499,500,501,510,512,513,514,515,518,519,520,522,523,526,528,529,530,531,532,533,534,535,687,695,696,772));
                int recordId = Integer.valueOf(modElem.attributeValue("record_id"));
                double mass = Double.valueOf(modElem.element("delta").attributeValue("mono_mass"));
                if (mass < -250 || mass > 250) continue;

                for (Element spec : modElem.elements("specificity")) {
                    String classification = spec.attributeValue("classification");
//                    if (classification.contentEquals("Isotopic label") || classification.contentEquals("AA substitution")) continue;
                    if (!recordIdWithPsiName.contains(Integer.valueOf(modElem.attributeValue("record_id"))) && !classification.contentEquals("AA substitution")) continue;
                    if (classification.contentEquals("Isotopic label") || classification.contains("glycos") || classification.contains("Other")) continue;

                    char site = spec.attributeValue("site").charAt(0);
                    boolean onlyProteinTerminalIfnc = false;
                    String position = spec.attributeValue("position");
                    if (spec.attributeValue("site").contentEquals("N-term")) {
                        site = 'n';
                        if (position.contentEquals("Protein N-term")) onlyProteinTerminalIfnc = true;
                    } else if (spec.attributeValue("site").contentEquals("C-term")) {
                        site = 'c';
                        if (position.contentEquals("Protein C-term")) onlyProteinTerminalIfnc = true;
                    }
                    if (position == null) {
                        int a = 1;
                    }
                    VarModParam temp = new VarModParam(mass, site, position, onlyProteinTerminalIfnc, classification, name);

                    if (siteModMap.containsKey(site)) {
                        siteModMap.get(site).add(temp);
                    } else {
                        Set<VarModParam> tempSet = new HashSet<>();
                        tempSet.add(temp);
                        siteModMap.put(site, tempSet);
                    }
                }
            }
        }
        return siteModMap;
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

    private Map<Character, Set<VarModParam>> buildAASMap(Multimap<Character, ModEntry> siteModMap) {
        int[][] pam1Matrix = new int[][]{
                {9867 , 2    , 9    , 10   , 3    , 8    , 17   , 21   , 2    , 6    , 4    , 2    , 6    , 2    , 22   , 35   , 32   , 0    , 2    , 18},
                {1    , 9913 , 1    , 0    , 1    , 10   , 0    , 0    , 10   , 3    , 1    , 19   , 4    , 1    , 4    , 6    , 1    , 8    , 0    , 1},
                {4    , 1    , 9822 , 36   , 0    , 4    , 6    , 6    , 21   , 3    , 1    , 13   , 0    , 1    , 2    , 20   , 9    , 1    , 4    , 1},
                {6    , 0    , 42   , 9859 , 0    , 6    , 53   , 6    , 4    , 1    , 0    , 3    , 0    , 0    , 1    , 5    , 3    , 0    , 0    , 1},
                {1    , 1    , 0    , 0    , 9973 , 0    , 0    , 0    , 1    , 1    , 0    , 0    , 0    , 0    , 1    , 5    , 1    , 0    , 3    , 2},
                {3    , 9    , 4    , 5    , 0    , 9876 , 27   , 1    , 23   , 1    , 3    , 6    , 4    , 0    , 6    , 2    , 2    , 0    , 0    , 1},
                {10   , 0    , 7    , 56   , 0    , 35   , 9865 , 4    , 2    , 3    , 1    , 4    , 1    , 0    , 3    , 4    , 2    , 0    , 1    , 2},
                {21   , 1    , 12   , 11   , 1    , 3    , 7    , 9935 , 1    , 0    , 1    , 2    , 1    , 1    , 3    , 21   , 3    , 0    , 0    , 5},
                {1    , 8    , 18   , 3    , 1    , 20   , 1    , 0    , 9912 , 0    , 1    , 1    , 0    , 2    , 3    , 1    , 1    , 1    , 4    , 1},
                {2    , 2    , 3    , 1    , 2    , 1    , 2    , 0    , 0    , 9872 , 9    , 2    , 12   , 7    , 0    , 1    , 7    , 0    , 1    , 33},
                {3    , 1    , 3    , 0    , 0    , 6    , 1    , 1    , 4    , 22   , 9947 , 2    , 45   , 13   , 3    , 1    , 3    , 4    , 2    , 15},
                {2    , 37   , 25   , 6    , 0    , 12   , 7    , 2    , 2    , 4    , 1    , 9926 , 20   , 0    , 3    , 8    , 11   , 0    , 1    , 1},
                {1    , 1    , 0    , 0    , 0    , 2    , 0    , 0    , 0    , 5    , 8    , 4    , 9874 , 1    , 0    , 1    , 2    , 0    , 0    , 4},
                {1    , 1    , 1    , 0    , 0    , 0    , 0    , 1    , 2    , 8    , 6    , 0    , 4    , 9946 , 0    , 2    , 1    , 3    , 28   , 0},
                {13   , 5    , 2    , 1    , 1    , 8    , 3    , 2    , 5    , 1    , 2    , 2    , 1    , 1    , 9926 , 12   , 4    , 0    , 0    , 2},
                {28   , 11   , 34   , 7    , 11   , 4    , 6    , 16   , 2    , 2    , 1    , 7    , 4    , 3    , 17   , 9840 , 38   , 5    , 2    , 2},
                {22   , 2    , 13   , 4    , 1    , 3    , 2    , 2    , 1    , 11   , 2    , 8    , 6    , 1    , 5    , 32   , 9871 , 0    , 2    , 9},
                {0    , 2    , 0    , 0    , 0    , 0    , 0    , 0    , 0    , 0    , 0    , 0    , 0    , 1    , 0    , 1    , 0    , 9976 , 1    , 0},
                {1    , 0    , 3    , 0    , 3    , 0    , 1    , 0    , 4    , 1    , 1    , 0    , 0    , 21   , 0    , 1    , 1    , 2    , 9945 , 1},
                {13   , 2    , 1    , 1    , 3    , 2    , 2    , 3    , 3    , 57   , 11   , 1    , 17   , 1    , 3    , 2    , 10   , 0    , 2    , 9901},
        };
        char[] aaArray = new char[]{'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'};
        Map<Character, Set<VarModParam>> aasMap = new HashMap<>();
        for (int i = 0; i < aaArray.length; ++i) {
            for (int j = 0; j < aaArray.length; ++j) {
                if (i != j) {
                    if (!(aaArray[i] == 'I' && aaArray[j] == 'L') && !(aaArray[i] == 'L' && aaArray[j] == 'I')) { // "I" and "L" have the same mass. don't consider such a AA substitution.
                        double ptmDeltaMass = massTable.get(aaArray[j]) - massTable.get(aaArray[i]);
                        if (ptmDeltaMass >= minPtmMass && ptmDeltaMass <= maxPtmMass) {
                            VarModParam temp = new VarModParam(ptmDeltaMass, aaArray[i], 0, true);
                            if (aasMap.containsKey(aaArray[i])) {
                                aasMap.get(aaArray[i]).add(temp); // all AAs have the same priority
                            } else {
                                Set<VarModParam> tempSet = new HashSet<>();
                                tempSet.add(temp);
                                aasMap.put(aaArray[i], tempSet);
                            }

                            siteModMap.put(aaArray[i], new ModEntry("AAS", aaArray[i] + "->" + aaArray[j], ptmDeltaMass)); // add the amino acid substitution to the map.
                        }
                    }
                }
            }
        }
        return aasMap;
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

    private Map<Integer, Set<VarModParam>> getIdxVarModMap(String ptmFreePeptide, Set<Integer> fixModIdxes, char leftFlank, char rightFlank) {
        Map<Integer, Set<VarModParam>> idxVarModMap = new HashMap<>(ptmFreePeptide.length() + 1, 1);
        for (int i = 0; i < ptmFreePeptide.length(); ++i) {
            if (!fixModIdxes.contains(i)) {
                char aa = ptmFreePeptide.charAt(i);
                if (aa == 'n') {
                    if (finalPtmMap.containsKey('n')) {
                        Set<VarModParam> tempSet = new HashSet<>();
                        for (VarModParam modEntry : finalPtmMap.get('n')) {
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
                        Set<VarModParam> tempSet = new HashSet<>();
                        for (VarModParam modEntry : finalPtmMap.get('c')) {
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

    private Map<Integer, Set<VarModParam>> getIdxVarModMapNew(String ptmFreePeptide, Set<Integer> fixModIdxes, char leftFlank, char rightFlank) {
        Map<Integer, Set<VarModParam>> idxVarModMap = new HashMap<>(ptmFreePeptide.length() + 1, 1);
        for (int i = 0; i < ptmFreePeptide.length(); ++i) {
            if (!fixModIdxes.contains(i)) {
                char aa = ptmFreePeptide.charAt(i);
                if (aa == 'n') {
                    if (finalPtmMap.containsKey('n')) {
                        Set<VarModParam> tempSet = new HashSet<>();
                        for (VarModParam modEntry : finalPtmMap.get('n')) {
                            if ( (leftFlank != '-' && !modEntry.position.contentEquals("Protein N-term"))
                                    || (leftFlank == '-')) {
                                if ((leftFlank != 'K' || Math.abs(massTable.get('K') - modEntry.mass) > ms2Tolerance)
                                        && (leftFlank != 'R' || Math.abs(massTable.get('R') - modEntry.mass) > ms2Tolerance)
                                        && (massTable.get(ptmFreePeptide.charAt(1)) + modEntry.mass > ms2Tolerance)) {
                                    // Fixing missed cleavages caused issue in N-term and the mass of a modified amino acid cannot be 0 or negative. lsz ????
                                    int a = 1;
                                }
                                if (massTable.get(ptmFreePeptide.charAt(1)) + modEntry.mass > ms2Tolerance){
                                    // Fixing missed cleavages caused issue in N-term and the mass of a modified amino acid cannot be 0 or negative. lsz ????
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
                        Set<VarModParam> tempSet = new HashSet<>();
                        for (VarModParam modEntry : finalPtmMap.get('c')) {
                            if ( (rightFlank != '-' && !modEntry.position.contentEquals("Protein C-term"))
                                    || (rightFlank == '-')) {
                                if ((rightFlank != 'K' || Math.abs(massTable.get('K') - modEntry.mass) > ms2Tolerance)
                                        && (rightFlank != 'R' || Math.abs(massTable.get('R') - modEntry.mass) > ms2Tolerance)
                                        && (massTable.get(ptmFreePeptide.charAt(ptmFreePeptide.length() - 2)) + modEntry.mass > ms2Tolerance)) {
                                    // Fixing missed cleavages caused issue in C-term and the mass of a modified amino acid cannot be 0 or negative
                                }
                                if ((massTable.get(ptmFreePeptide.charAt(ptmFreePeptide.length() - 2)) + modEntry.mass > ms2Tolerance)) {
                                    // Fixing missed cleavages caused issue in C-term and the mass of a modified amino acid cannot be 0 or negative
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
                        Set<VarModParam> tempSet = new HashSet<>();
                        for (VarModParam modEntry : finalPtmMap.get(aa)) {
                            if (modEntry.position == null) {
                                int a = 1;
                            }
                            if ((massTable.get(ptmFreePeptide.charAt(i)) + modEntry.mass <= ms2Tolerance)) {
                                // Fixing missed cleavages caused issue in C-term and the mass of a modified amino acid cannot be 0 or negative
                                continue;
                            }
                            if (modEntry.position.contentEquals("Anywhere")) {
                                tempSet.add(modEntry);
                            }
                            if (modEntry.position.contentEquals("Any N-term") && (i == 1)) {
                                tempSet.add(modEntry);
                            }
                            if (modEntry.position.contentEquals("Any C-term") && (i == ptmFreePeptide.length()-2)) {
                                tempSet.add(modEntry);
                            }
                            if (modEntry.position.contentEquals("Protein N-term") && (i == 1) && (leftFlank == '_')) {
                                tempSet.add(modEntry);
                            }
                            if (modEntry.position.contentEquals("Any C-term") && (i == ptmFreePeptide.length()-2) && (rightFlank == '_')) {
                                tempSet.add(modEntry);
                            }
                        }
                        if (!tempSet.isEmpty()) {
                            idxVarModMap.put(i, tempSet);
                        }
                    }
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
