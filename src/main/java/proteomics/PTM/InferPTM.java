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
import proteomics.OutputPeff;
import proteomics.Types.*;

import java.io.*;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;


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
        siteModMap = readUnimodAndGenerateAAS(minPtmMass, maxPtmMass);

        // Building an amino acid substitution matrix.
        Map<Character, Set<VarModParam>> aasMap = buildAASMap(siteModMap);

        // Reading Mod table...
        Map<Character, Set<VarModParam>> ptmMap = readModFile();

        // update ptm table with the high priority mods in the parameter file.
        for (VarModParam varModParam : varModParamSet) {
            if (ptmMap.containsKey(varModParam.aa)) {
                ptmMap.get(varModParam.aa).remove(varModParam);
                ptmMap.get(varModParam.aa).add(varModParam);
            } else {
                Set<VarModParam> tempSet = new HashSet<>();
                tempSet.add(varModParam);
                ptmMap.put(varModParam.aa, tempSet);
            }
        }

        // merge PTM map and AA substitution map
        for (char aa : aasMap.keySet()) {
            for (VarModParam modEntry : aasMap.get(aa)) {
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
    }

    public PeptidePTMPattern tryPTM(SparseVector expProcessedPL, TreeMap<Double, Double> plMap, double precursorMass, String ptmFreePeptide, boolean isDecoy, double normalizedCrossCorr, char leftFlank, char rightFlank, int globalRank, int precursorCharge, int localMaxMS2Charge, double localMS1ToleranceL, double localMS1ToleranceR) {
        double ptmFreeMass = massTool.calResidueMass(ptmFreePeptide) + massTool.H2O;
        double deltaMass = precursorMass - ptmFreeMass;
        double leftMassBound = deltaMass + localMS1ToleranceL;
        double rightMassBound = deltaMass + localMS1ToleranceR;
        Set<Integer> fixModIdxes = getFixModIdxes(ptmFreePeptide, fixModMap);

        PeptidePTMPattern peptidePTMPattern = new PeptidePTMPattern(ptmFreePeptide);

        // try different PTM combinations
        Set<String> checkedPtmPattern1 = new HashSet<>();
        Set<String> checkedPtmPattern2 = new HashSet<>();
        Set<String> checkedPtmPattern3 = new HashSet<>();
        Set<String> checkedPtmPattern4 = new HashSet<>();
        Set<String> checkedPtmPattern5 = new HashSet<>();

        Map<Integer, Set<VarModParam>> idxVarModMap = getIdxVarModMap(ptmFreePeptide, fixModIdxes, leftFlank, rightFlank);

        try1PTMs(idxVarModMap, leftMassBound, rightMassBound, ptmFreePeptide, isDecoy, normalizedCrossCorr, globalRank, checkedPtmPattern1, peptidePTMPattern, expProcessedPL, plMap, precursorCharge, localMaxMS2Charge);

        if (idxVarModMap.size() > 1) {
            // Try 2 PTMs
            try2PTMs(idxVarModMap, leftMassBound, rightMassBound, ptmFreePeptide, isDecoy, normalizedCrossCorr, globalRank, checkedPtmPattern2, peptidePTMPattern, expProcessedPL, plMap, precursorCharge, localMaxMS2Charge);
        }

        if (idxVarModMap.size() > 2) {
            // Try 3 PTMs
            try3PTMs(idxVarModMap, leftMassBound, rightMassBound, ptmFreePeptide, isDecoy, normalizedCrossCorr, globalRank, checkedPtmPattern3, peptidePTMPattern, expProcessedPL, plMap, precursorCharge, localMaxMS2Charge);
        }

        if (idxVarModMap.size() > 3) {
            // Try 4 PTMs
            try4PTMs(idxVarModMap, leftMassBound, rightMassBound, ptmFreePeptide, isDecoy, normalizedCrossCorr, globalRank, checkedPtmPattern4, peptidePTMPattern, expProcessedPL, plMap, precursorCharge, localMaxMS2Charge);
        }

        if (idxVarModMap.size() > 4) {
            // Try 5 PTMs
            try5PTMs(idxVarModMap, leftMassBound, rightMassBound, ptmFreePeptide, isDecoy, normalizedCrossCorr, globalRank, checkedPtmPattern5, peptidePTMPattern, expProcessedPL, plMap, precursorCharge, localMaxMS2Charge);
        }

        return peptidePTMPattern;
    }

    public PeptidePTMPattern findPTM(int scanNum, SparseVector expProcessedPL, TreeMap<Double, Double> plMap, double precursorMass, String ptmFreePeptide, boolean isDecoy, double normalizedCrossCorr, char leftFlank, char rightFlank, int globalRank, int precursorCharge, int localMaxMS2Charge, double localMS1ToleranceL, double localMS1ToleranceR, List<ThreeExpAA> expAaLists) {
//        double ptmFreeMass = massTool.calResidueMass(ptmFreePeptide) + massTool.H2O;
        PeptidePTMPattern peptidePTMPattern = new PeptidePTMPattern(ptmFreePeptide);
        Peptide peptide = new Peptide(ptmFreePeptide, isDecoy, massTool, localMaxMS2Charge, normalizedCrossCorr, globalRank);
        double totalDeltaMass = precursorMass - peptide.getTheoMass();
        if (Math.abs(totalDeltaMass) < 0.01 ) {
            return peptidePTMPattern;
        }
        Set<Integer> fixModIdxes = getFixModIdxes(ptmFreePeptide, fixModMap);  // record the indexes of the amino acids with fix mod, e.g. when there is no C, this is empty set.

//        double a = peptide.getTheoMass();
        double[][] theoIonsMatrix = peptide.getTheoIonMatrix();// Note that if you set varMod on it, it returns b y ions with var mod
        String ptmFreePeptideOrdinary = ptmFreePeptide.replaceAll("I","#").replaceAll("L","#");

        Map<Integer, Set<VarModParam>> idxVarModMap = getIdxVarModMap(ptmFreePeptide, fixModIdxes, leftFlank, rightFlank);
        Map<Integer, VarModParam[]> idxVarModArrayMap = new HashMap<>();
        for (int id : idxVarModMap.keySet()){
            VarModParam[] modArray = new VarModParam[idxVarModMap.get(id).size()];
            idxVarModMap.get(id).toArray(modArray);
            Arrays.sort(modArray, Comparator.comparingDouble(VarModParam::getMass));
            idxVarModArrayMap.put(id, modArray);
        }
        Set<Integer> modifiedZone = new HashSet<>(idxVarModArrayMap.keySet());// e.g. nABCDEFc, modifiedZone={1,2,3,4,5,6}

        expAaLists.sort(Comparator.comparing(ThreeExpAA::getTotalIntensity).reversed()); // sort tags by total intensity descending order

        Set<Integer> tag2Ids = new HashSet<>();

        //find feasible zone, clean tags and M_modTags, only for b ions now, wait for Y ions
        Map<Double, Set<Integer>> extraMassTagIdMap = new HashMap<>(); // records m-mod tags and the tag ids.
        Set<String> seqOfGoodTags = new HashSet<>();
        Set<Integer> notGoodTag3IdSet = new HashSet<>();

        int maxTailPosOfCleanTags = 0;
        int minHeadPosOfMmodTags = Collections.max(modifiedZone) + 1;
        double cleanZoneScore = 0d;
        double tagZoneScore = 0d;
        for (int tagId = 0; tagId < expAaLists.size(); tagId++) {
            ThreeExpAA tagInfo = expAaLists.get(tagId);
            String tagSeq = tagInfo.getPtmFreeAAString();
            int alignBPos = ptmFreePeptideOrdinary.indexOf(tagSeq);
            if (alignBPos == -1) {
                if (ptmFreePeptideOrdinary.indexOf(tagSeq.substring(0, 2)) != -1) {
                    tag2Ids.add(tagId); // if no clean tag or M-mod tag is found, use tag 2 to narrow down modifiedZone
                }
                continue;
            }

            while (alignBPos != -1) {
                double deltaMass = tagInfo.getHeadLocation() - theoIonsMatrix[0][alignBPos - 1];
                if (Math.abs(deltaMass) <= 0.01 && alignBPos + 3 <= Collections.max(modifiedZone)) {
                    //if the tag is at the C term, it can not be right so skip.
                    if (alignBPos + 3 - 1 > maxTailPosOfCleanTags) {
                        maxTailPosOfCleanTags = alignBPos + 3 - 1;
                    }
                    cleanZoneScore += tagInfo.getTotalIntensity();
                    tagInfo.isGoodTag3 = true;
                    seqOfGoodTags.add(tagSeq);
                    break;
                } else if (Math.abs(deltaMass - totalDeltaMass) <= 0.01 && alignBPos > Collections.min(modifiedZone)) {
                    //skip if it has ptm but it is at the n ter
                    if (alignBPos < minHeadPosOfMmodTags) {
                        minHeadPosOfMmodTags = alignBPos;
                    }
                    tagZoneScore += tagInfo.getTotalIntensity();
                    tagInfo.isGoodTag3 = true;
                    seqOfGoodTags.add(tagSeq);
                    break;
                } else {
                    tagInfo.bAlignPosMassMap.put(alignBPos, deltaMass);
                    alignBPos = ptmFreePeptideOrdinary.indexOf(tagSeq, alignBPos + 1);
                }

            }
            if (!tagInfo.isGoodTag3) {
                notGoodTag3IdSet.add(tagId);
            }
        }

        for (int tagId : notGoodTag3IdSet) {
            ThreeExpAA tagInfo = expAaLists.get(tagId);
            if (seqOfGoodTags.contains(tagInfo.getPtmFreeAAString())) {
                continue;
            }
            for (int alignPos : tagInfo.bAlignPosMassMap.keySet()){
                double deltaMass = tagInfo.bAlignPosMassMap.get(alignPos);
                boolean isNewMass = true;
                for (double mass : extraMassTagIdMap.keySet()) {
                    if (Math.abs(deltaMass - mass) <= 0.01){

                        int newMassTimes = 1+extraMassTagIdMap.get(mass).size();
                        double newAveMass = (deltaMass+mass*(newMassTimes-1)) / newMassTimes;
                        Set<Integer> newMassIdSet = new HashSet<>(extraMassTagIdMap.get(mass));
                        newMassIdSet.add(tagId);
                        extraMassTagIdMap.remove(mass);
                        extraMassTagIdMap.put(newAveMass, newMassIdSet);
                        isNewMass = false;
                        break;
                    }
                }
                if (isNewMass) {
                    Set<Integer> massIdSet = new HashSet<>();
                    massIdSet.add(tagId);
                    extraMassTagIdMap.put(deltaMass, massIdSet);
                }
            }
        }

        Set<Integer> tempMiddleZone = IntStream.range(maxTailPosOfCleanTags+1, minHeadPosOfMmodTags).boxed().collect(Collectors.toSet());
        tempMiddleZone.retainAll(modifiedZone);
        if (!tempMiddleZone.isEmpty()){ // not interfering between tagZone and cleanZone
            for (int i = 0; i < maxTailPosOfCleanTags; i++){
                modifiedZone.remove(i);
            }
            for (int i = minHeadPosOfMmodTags; i < ptmFreePeptide.length()-1; i++){
                modifiedZone.remove(i);
            }
        } else {
            if (tagZoneScore > cleanZoneScore){
                for (int i = minHeadPosOfMmodTags; i < ptmFreePeptide.length()-1; i++){
                    modifiedZone.remove(i);
                }
            } else {
                for (int i = 0; i < maxTailPosOfCleanTags; i++){
                    modifiedZone.remove(i);
                }
            }
        }

        if (modifiedZone.isEmpty()) {
//            System.out.println(scanNum + " is empty modifiedZone before tag 2");
            modifiedZone.addAll(idxVarModArrayMap.keySet());
        } else if (modifiedZone.size() == idxVarModArrayMap.size() && extraMassTagIdMap.size() != 0) {
            // use tag2
//            System.out.println(scanNum + " has tag 3 but uses tag 2");
            for (int tagId : tag2Ids){
                ThreeExpAA tagInfo = expAaLists.get(tagId);
                String tagSeq = tagInfo.getPtmFreeAAString();
                StringBuilder sb = new StringBuilder(tagSeq);
                sb.reverse();
                String tagSeqReverse = sb.toString();
                String tagSeqPrefix = tagSeq.substring(0,2);
                int alignPos = ptmFreePeptideOrdinary.indexOf(tagSeqPrefix);
                if (alignPos == -1) {
                    continue; // pep does not contain the tag
                }
                int headPeakId = alignPos - 1;
                double deltaMass = tagInfo.getHeadLocation() - theoIonsMatrix[0][headPeakId];
                if (Math.abs(deltaMass) <= 0.01){
                    if (alignPos+2 > Collections.max(modifiedZone)){
                        continue; //if the tag is at the C term, it can not be right so skip.
                    }
                    for (int i = 1; i < alignPos+2; i++){
                        modifiedZone.remove(i);
                    }
                } else if(Math.abs(deltaMass - totalDeltaMass) <= 0.01){
                    if (alignPos <= Collections.min(modifiedZone)) {
                        continue;
                    }
                    for (int i = alignPos; i < ptmFreePeptide.length()-1; i++){
                        modifiedZone.remove(i);
                    }
                }
            }
        }

        String truthPep = "nGAEASAASEEEAGPQATEPSTPSGPESGPTPASAEQNEc";
        Peptide p1p = new Peptide(truthPep, isDecoy, massTool, localMaxMS2Charge, normalizedCrossCorr, globalRank);
        PositionDeltaMassMap fakePatten = new PositionDeltaMassMap(ptmFreePeptide.length());
        fakePatten.put(new Coordinate(18, 19), (-31.990));
        fakePatten.put(new Coordinate(22, 23), (33.988));
        p1p.setVarPTM(fakePatten);
        Map<Integer, Double> matchedBions1 = new HashMap<>();
        Map<Integer, Double> matchedYions1 = new HashMap<>();
        double[][] temp11 = p1p.getIonMatrixNow();
        Set<Integer> jRange1 = new HashSet<>();
        for (int i : modifiedZone) {
            jRange1.add(i-1);
        }
        double ss = massTool.buildVectorAndCalXCorr(p1p.getIonMatrixNow(), 1, expProcessedPL, matchedBions1, matchedYions1, IntStream.range(0, truthPep.length()-2).boxed().collect(Collectors.toSet())) ;
        //stub end

        if (modifiedZone.size() == 0) {
//            System.out.println(scanNum + " is empty modifiedZone after tag 2");
            return peptidePTMPattern; //Some scans are not valid Scans. Will be deleted soon.
        }

        TreeMap<Double, Double> unUsedPlMap = new TreeMap<>(plMap);

        PeptidePTMPattern allPtmPattern = new PeptidePTMPattern(ptmFreePeptide,1);
        PeptidePTMPattern allPtmPatternBad = new PeptidePTMPattern(ptmFreePeptide,1);
        modifiedZone = IntStream.range(1, ptmFreePeptide.length()-1).boxed().collect(Collectors.toSet());// e.g. nABCDEFc, modifiedZone={1,2,3,4,5,6}

        Peptide cleanPep = new Peptide(ptmFreePeptide, isDecoy, massTool, localMaxMS2Charge, normalizedCrossCorr, globalRank);

        int lb = 1;
        int rb = ptmFreePeptide.length() - 2;
        Map<Integer, Double> matchedBions = new HashMap<>();
        Map<Integer, Double> matchedYions = new HashMap<>();
        double[][] temp1 = cleanPep.getIonMatrixNow();
        Set<Integer> jRange = new HashSet<>();
        for (int i : modifiedZone) {
            jRange.add(i-1);
        }
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
            rb = Collections.min(matchedYions.keySet());
        }


//        matchedBIndex.addAll(matchedBions.keySet());
        if (rb - lb < 1) {
            double bSumIntens = 0;
            for (double intes : matchedBions.values()) bSumIntens += intes;
            double ySumIntens = 0;
            for (double intes : matchedYions.values()) ySumIntens += intes;
            if (bSumIntens > ySumIntens) {
                rb = ptmFreePeptide.length() - 2;
            } else {
                lb = 1;
            }
        }
        modifiedZone = IntStream.range(lb, rb).boxed().collect(Collectors.toSet());
//        modifiedZone = IntStream.range(1, ptmFreePeptide.length()-1).boxed().collect(Collectors.toSet());// e.g. nABCDEFc, modifiedZone={1,2,3,4,5,6}

        allPtmPatternBad.push(cleanPep);
        allPtmPattern.push(cleanPep);
        allPtmPatternBad.bestPep = allPtmPatternBad.getTopPepPtn();
        if (modifiedZone.isEmpty()) {
            return allPtmPatternBad;
        }
//
//        for (int j = rb - 2; j < cleanPepIonMatrix[0].length; j++) {
//            cleanPepIonMatrix[0][j] += totalDeltaMass;
//        }
//        for (int j = 0; j < lb; j++) {
//            if (j >= cleanPepIonMatrix[1].length){
//                int a = 1;
//            }
//            cleanPepIonMatrix[1][j] += totalDeltaMass;
//        }
//
//        Map<Integer, Double> matchedBionsNew = new HashMap<>();
//        Map<Integer, Double> matchedYionsNew = new HashMap<>();
//        double cleanScoreNew = massTool.buildVectorAndCalXCorr(cleanPepIonMatrix, 1, expProcessedPL, matchedBionsNew, matchedYionsNew, jRange) ;
//        if (cleanScoreNew > cleanScore) {
//            cleanPep.setScore(cleanScoreNew - cleanP);
//            cleanPep.setMatchedPeakNum(Score.getMatchedIonNum(plMap, localMaxMS2Charge, cleanPepIonMatrix, ms2Tolerance));
//        }



        PeptidePTMPattern ptmInitialTemp = new PeptidePTMPattern(ptmFreePeptide, 1);
        Peptide cleanPeptide = new Peptide(ptmFreePeptide, isDecoy, massTool, localMaxMS2Charge, normalizedCrossCorr, globalRank);
        cleanPeptide.setVarPTM(new PositionDeltaMassMap(ptmFreePeptide.length()));
        ptmInitialTemp.push(cleanPeptide);
//        long t1 = System.currentTimeMillis();

        // 1st
        PeptidePTMPattern ptmN1Good = new PeptidePTMPattern(ptmFreePeptide, 1);
        PeptidePTMPattern ptmN1Bad = new PeptidePTMPattern(ptmFreePeptide, 1);

        DividedZone z1Z2Res = dividePep(scanNum, modifiedZone, ptmN1Good, ptmN1Bad, ptmInitialTemp, idxVarModArrayMap, totalDeltaMass, ptmFreePeptide, isDecoy, normalizedCrossCorr, globalRank, expProcessedPL, plMap, precursorCharge, localMaxMS2Charge, unUsedPlMap);
        if (!ptmN1Good.getPeptideTreeSet().isEmpty()) {
            ptmN1Good.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN1Good.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 1*0.1);
            ptmN1Good.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, localMaxMS2Charge, ptmN1Good.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
            allPtmPattern.push(ptmN1Good.getTopPepPtn());
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
        ptmN1Bad.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN1Bad.getTopPepPtn().getIonMatrixNow(), 1,  expProcessedPL) - 1*0.1);
        ptmN1Bad.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, localMaxMS2Charge, ptmN1Bad.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
        allPtmPatternBad.push(ptmN1Bad.getTopPepPtn());

        double zone1Mass = ptmN1Bad.getTopPepPtn().getVarPTMs().values().iterator().next();
        double zone2Mass = totalDeltaMass - zone1Mass;

        // 2nd
        PeptidePTMPattern ptmN2Good = new PeptidePTMPattern(ptmFreePeptide, 1);
        PeptidePTMPattern ptmN2Bad = new PeptidePTMPattern(ptmFreePeptide, 1);

        DividedZone z3Z4Res = dividePep(scanNum, zone2, ptmN2Good, ptmN2Bad ,ptmN1Bad, idxVarModArrayMap, zone2Mass, ptmFreePeptide, isDecoy, normalizedCrossCorr, globalRank, expProcessedPL, plMap, precursorCharge, localMaxMS2Charge, unUsedPlMap);
        if (!ptmN2Good.getPeptideTreeSet().isEmpty()) {
            ptmN2Good.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN2Good.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 2*0.1);
            ptmN2Good.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, localMaxMS2Charge, ptmN2Good.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
            allPtmPattern.push(ptmN2Good.getTopPepPtn());
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
        ptmN2Bad.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN2Bad.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 2*0.1);
        ptmN2Bad.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, localMaxMS2Charge, ptmN2Bad.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
        allPtmPatternBad.push(ptmN2Bad.getTopPepPtn());

        double zone3Mass = z3Z4Res.ptmMass;
        double zone4Mass = zone2Mass - zone3Mass;

        // 3rd
        PeptidePTMPattern ptmN3Good = new PeptidePTMPattern(ptmFreePeptide, 1);
        PeptidePTMPattern ptmN3Bad = new PeptidePTMPattern(ptmFreePeptide, 1);

        DividedZone z5Z6Res = dividePep(scanNum, zone4, ptmN3Good, ptmN3Bad, ptmN2Bad, idxVarModArrayMap, zone4Mass, ptmFreePeptide, isDecoy, normalizedCrossCorr, globalRank, expProcessedPL, plMap, precursorCharge, localMaxMS2Charge, unUsedPlMap);
        if (!ptmN3Good.peptideTreeSet.isEmpty()) {
            ptmN3Good.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN3Good.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 3*0.1);
            ptmN3Good.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, localMaxMS2Charge, ptmN3Good.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
            allPtmPattern.push(ptmN3Good.getTopPepPtn());
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
        ptmN3Bad.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN3Bad.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 3*0.1);
        ptmN3Bad.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, localMaxMS2Charge, ptmN3Bad.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
        allPtmPatternBad.push(ptmN3Bad.getTopPepPtn());

        double zone5Mass = z5Z6Res.ptmMass;
        double zone6Mass = zone4Mass - zone5Mass;

        PeptidePTMPattern ptmN4Good = new PeptidePTMPattern(ptmFreePeptide, 1);
        PeptidePTMPattern ptmN4Bad = new PeptidePTMPattern(ptmFreePeptide, 1);
        DividedZone z7Z8Res = dividePep(scanNum, zone6, ptmN4Good, ptmN4Bad ,ptmN3Bad, idxVarModArrayMap, zone6Mass, ptmFreePeptide, isDecoy, normalizedCrossCorr, globalRank, expProcessedPL, plMap, precursorCharge, localMaxMS2Charge, unUsedPlMap);
        if (!ptmN4Good.getPeptideTreeSet().isEmpty()) {
            ptmN4Good.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN4Good.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 4*0.1);
            ptmN4Good.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, localMaxMS2Charge, ptmN4Good.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
            allPtmPattern.push(ptmN4Good.getTopPepPtn());
        }
        if (!ptmN4Bad.getPeptideTreeSet().isEmpty()) {
            ptmN4Bad.getTopPepPtn().setScore(massTool.buildVectorAndCalXCorr(ptmN4Bad.getTopPepPtn().getIonMatrixNow(), 1, expProcessedPL) - 4*0.1);
            ptmN4Bad.getTopPepPtn().setMatchedPeakNum(Score.getMatchedIonNum(plMap, localMaxMS2Charge, ptmN4Bad.getTopPepPtn().getIonMatrixNow(), ms2Tolerance));
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
                Set<Integer> jRange = new HashSet<>();
                for (int i : modifiedZone) {
                    jRange.add(i-1);
                }
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

    public static Multimap<Character, ModEntry> readUnimodAndGenerateAAS(double minPtmMass, double maxPtmMass) throws IOException {
        Multimap<Character, ModEntry> siteModMap = HashMultimap.create();
        InputStream inputStream = OutputPeff.class.getClassLoader().getResourceAsStream("unimod.xml.tsv");
        BufferedReader reader = new BufferedReader(new InputStreamReader(inputStream));
        String line;
        while ((line = reader.readLine()) != null) {
            line = line.trim();
            if (!line.isEmpty() && !line.startsWith("accession")) {
                String[] parts = line.split("\t");
                if (!parts[6].contentEquals("Isotopic label") && !parts[6].contentEquals("AA substitution")) { // We don't consider isotopic label and amino acid substitution..
                    char site = parts[1].charAt(0);
                    if (parts[1].trim().contentEquals("N-term")) {
                        site = 'n';
                    } else if (parts[1].trim().contentEquals("C-term")) {
                        site = 'c';
                    }
                    if (site != 'n' && site != 'c') { // We don't consider peptide terminal modifications
                        double ptmDeltaMass =  Double.valueOf(parts[4].trim());
                        if (ptmDeltaMass >= minPtmMass && ptmDeltaMass <= maxPtmMass) { // only record the PTM within the delta mass range
                            siteModMap.put(site, new ModEntry("UNIMOD:" + parts[0].trim(), parts[2].trim().contentEquals("null") ? parts[3].trim() : parts[2].trim(), ptmDeltaMass)); // if there is no PSI-MS name, use the internal name in the Unimod
                        }
                    }
                }
            }
        }
        reader.close();

        return siteModMap;
    }

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




    private void try1PTMs(Map<Integer, Set<VarModParam>> idxVarModMap, double leftMassBound, double rightMassBound, String ptmFreePeptide, boolean isDecoy, double normalizedCrossCorr, int globalRank, Set<String> checkedPtmPattern, PeptidePTMPattern peptidePTMPattern, SparseVector expProcessedPL, TreeMap<Double, Double> plMap, int precursorCharge, int localMaxMS2Charge) { // Sometimes, the precursor mass error may affects the digitized spectrum.
        Integer[] idxArray = idxVarModMap.keySet().toArray(new Integer[0]);
        Arrays.sort(idxArray);
        for (int i = 0; i < idxArray.length - 1; ++i) {
            for (VarModParam modEntry1 : idxVarModMap.get(idxArray[i])) {
                if (modEntry1.mass <= rightMassBound && modEntry1.mass >= leftMassBound) {
                    if (!checkedPtmPattern.contains(idxArray[i] + "-" + Math.round(modEntry1.mass * 1000))) {
                        checkedPtmPattern.add(idxArray[i] + "-" + Math.round(modEntry1.mass * 1000));//use checkedPtmPattern to avoid recording repeat pattern, and mass*1000 as an integer to be key.
                        PositionDeltaMassMap positionDeltaMassMap = new PositionDeltaMassMap(ptmFreePeptide.length());
                        positionDeltaMassMap.put(new Coordinate(idxArray[i], idxArray[i] + 1), modEntry1.mass);
                        Peptide peptide = new Peptide(ptmFreePeptide, isDecoy, massTool, localMaxMS2Charge, normalizedCrossCorr, globalRank);
                        peptide.setVarPTM(positionDeltaMassMap); // since this is try 1 PTM, it is easy to traverse all options and find the top 5
                        double score = massTool.buildVectorAndCalXCorr(peptide.getIonMatrix(), precursorCharge, expProcessedPL);
                        if (score > 0) {
                            peptide.setScore(score);
                            peptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, localMaxMS2Charge, peptide.getIonMatrix(), ms2Tolerance));
                            peptidePTMPattern.push(peptide);
                        }
                    }
                }
            }
        }
    }



    private void try2PTMs(Map<Integer, Set<VarModParam>> idxVarModMap, double leftMassBound, double rightMassBound, String ptmFreePeptide, boolean isDecoy, double normalizedCrossCorr, int globalRank, Set<String> checkedPtmPattern, PeptidePTMPattern peptidePTMPattern, SparseVector expProcessedPL, TreeMap<Double, Double> plMap, int precursorCharge, int localMaxMS2Charge) {
        Integer[] idxArray = idxVarModMap.keySet().toArray(new Integer[0]);
        Arrays.sort(idxArray);
        for (int i = 0; i < idxArray.length - 1; ++i) {
            for (VarModParam modEntry1 : idxVarModMap.get(idxArray[i])) {
                for (int j = i + 1; j < idxArray.length; ++j) {
                    for (VarModParam modEntry2 : idxVarModMap.get(idxArray[j])) {
                        if (modEntry1.mass + modEntry2.mass <= rightMassBound && modEntry1.mass + modEntry2.mass >= leftMassBound) {
                            if (!checkedPtmPattern.contains(idxArray[i] + "-" + Math.round(modEntry1.mass * 1000) + "-" + idxArray[j] + "-" + Math.round(modEntry2.mass * 1000))) {
                                checkedPtmPattern.add(idxArray[i] + "-" + Math.round(modEntry1.mass * 1000) + "-" + idxArray[j] + "-" + Math.round(modEntry2.mass * 1000));
                                PositionDeltaMassMap positionDeltaMassMap = new PositionDeltaMassMap(ptmFreePeptide.length());
                                positionDeltaMassMap.put(new Coordinate(idxArray[i], idxArray[i] + 1), modEntry1.mass);
                                positionDeltaMassMap.put(new Coordinate(idxArray[j], idxArray[j] + 1), modEntry2.mass);
                                Peptide peptide = new Peptide(ptmFreePeptide, isDecoy, massTool, localMaxMS2Charge, normalizedCrossCorr, globalRank);
                                peptide.setVarPTM(positionDeltaMassMap);
                                double score = massTool.buildVectorAndCalXCorr(peptide.getIonMatrix(), precursorCharge, expProcessedPL);
                                if (score > 0) {
                                    peptide.setScore(score);
                                    peptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, localMaxMS2Charge, peptide.getIonMatrix(), ms2Tolerance));
                                    peptidePTMPattern.push(peptide);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    private void try3PTMs(Map<Integer, Set<VarModParam>> idxVarModMap, double leftMassBound, double rightMassBound, String ptmFreePeptide, boolean isDecoy, double normalizedCrossCorr, int globalRank, Set<String> checkedPtmPattern, PeptidePTMPattern peptidePTMPattern, SparseVector expProcessedPL, TreeMap<Double, Double> plMap, int precursorCharge, int localMaxMS2Charge) {
        Integer[] idxArray = idxVarModMap.keySet().toArray(new Integer[0]);
        Arrays.sort(idxArray);
        for (int i = 0; i < idxArray.length - 2; ++i) {
            for (VarModParam modEntry1 : idxVarModMap.get(idxArray[i])) {
                for (int j = i + 1; j < idxArray.length - 1; ++j) {
                    for (VarModParam modEntry2 : idxVarModMap.get(idxArray[j])) {
                        if (modEntry1.priority + modEntry2.priority > 0) {
                            if (Math.abs(modEntry1.mass + modEntry2.mass) >= ptmMassTolerance) { // two self cancelled PTM masses are not allowed.
                                for (int k = j + 1; k < idxArray.length; ++k) {
                                    for (VarModParam modEntry3 : idxVarModMap.get(idxArray[k])) {
                                        if (modEntry1.priority + modEntry2.priority + modEntry3.priority > 1) {
                                            if (Math.abs(modEntry1.mass + modEntry3.mass) >= ptmMassTolerance && Math.abs(modEntry2.mass + modEntry3.mass) >= ptmMassTolerance) {
                                                if (modEntry1.mass + modEntry2.mass + modEntry3.mass <= rightMassBound && modEntry1.mass + modEntry2.mass + modEntry3.mass >= leftMassBound) {
                                                    if (!checkedPtmPattern.contains(idxArray[i] + "-" + Math.round(modEntry1.mass * 1000) + "-" + idxArray[j] + "-" + Math.round(modEntry2.mass * 1000) + "-" + idxArray[k] + "-" + Math.round(modEntry3.mass * 1000))) {
                                                        checkedPtmPattern.add(idxArray[i] + "-" + Math.round(modEntry1.mass * 1000) + "-" + idxArray[j] + "-" + Math.round(modEntry2.mass * 1000) + "-" + idxArray[k] + "-" + Math.round(modEntry3.mass * 1000));
                                                        PositionDeltaMassMap positionDeltaMassMap = new PositionDeltaMassMap(ptmFreePeptide.length());
                                                        positionDeltaMassMap.put(new Coordinate(idxArray[i], idxArray[i] + 1), modEntry1.mass);
                                                        positionDeltaMassMap.put(new Coordinate(idxArray[j], idxArray[j] + 1), modEntry2.mass);
                                                        positionDeltaMassMap.put(new Coordinate(idxArray[k], idxArray[k] + 1), modEntry3.mass);
                                                        Peptide peptide = new Peptide(ptmFreePeptide, isDecoy, massTool, localMaxMS2Charge, normalizedCrossCorr, globalRank);
                                                        peptide.setVarPTM(positionDeltaMassMap);
                                                        double score = massTool.buildVectorAndCalXCorr(peptide.getIonMatrix(), precursorCharge, expProcessedPL);
                                                        if (score > 0) {
                                                            peptide.setScore(score);
                                                            peptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, localMaxMS2Charge, peptide.getIonMatrix(), ms2Tolerance));
                                                            peptidePTMPattern.push(peptide);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    private void try4PTMs(Map<Integer, Set<VarModParam>> idxVarModMap, double leftMassBound, double rightMassBound, String ptmFreePeptide, boolean isDecoy, double normalizedCrossCorr, int globalRank, Set<String> checkedPtmPattern, PeptidePTMPattern peptidePTMPattern, SparseVector expProcessedPL, TreeMap<Double, Double> plMap, int precursorCharge, int localMaxMS2Charge) { // only allow one low priority PTM
        Integer[] idxArray = idxVarModMap.keySet().toArray(new Integer[0]);
        Arrays.sort(idxArray);
        for (int i = 0; i < idxArray.length - 3; ++i) {
            for (VarModParam modEntry1 : idxVarModMap.get(idxArray[i])) {
                for (int j = i + 1; j < idxArray.length - 2; ++j) {
                    for (VarModParam modEntry2 : idxVarModMap.get(idxArray[j])) {
                        if (modEntry1.priority + modEntry2.priority > 0) {
                            if (Math.abs(modEntry1.mass + modEntry2.mass) >= ptmMassTolerance) {
                                for (int k = j + 1; k < idxArray.length - 1; ++k) {
                                    for (VarModParam modEntry3 : idxVarModMap.get(idxArray[k])) {
                                        if (modEntry1.priority + modEntry2.priority + modEntry3.priority > 1) {
                                            if (Math.abs(modEntry1.mass + modEntry3.mass) >= ptmMassTolerance && Math.abs(modEntry2.mass + modEntry3.mass) >= ptmMassTolerance) {
                                                for (int l = k + 1; l < idxArray.length; ++l) {
                                                    for (VarModParam modEntry4 : idxVarModMap.get(idxArray[l])) {
                                                        if (modEntry1.priority + modEntry2.priority + modEntry3.priority + modEntry4.priority > 2) {
                                                            if (Math.abs(modEntry1.mass + modEntry4.mass) >= ptmMassTolerance && Math.abs(modEntry2.mass + modEntry4.mass) >= ptmMassTolerance && Math.abs(modEntry3.mass + modEntry4.mass) >= ptmMassTolerance) {
                                                                if (modEntry1.mass + modEntry2.mass + modEntry3.mass + modEntry4.mass <= rightMassBound && modEntry1.mass + modEntry2.mass + modEntry3.mass + modEntry4.mass >= leftMassBound) {
                                                                    if (!checkedPtmPattern.contains(idxArray[i] + "-" + Math.round(modEntry1.mass * 1000) + "-" + idxArray[j] + "-" + Math.round(modEntry2.mass * 1000) + "-" + idxArray[k] + "-" + Math.round(modEntry3.mass * 1000) + "-" + idxArray[l] + "-" + Math.round(modEntry4.mass * 1000))) {
                                                                        checkedPtmPattern.add(idxArray[i] + "-" + Math.round(modEntry1.mass * 1000) + "-" + idxArray[j] + "-" + Math.round(modEntry2.mass * 1000) + "-" + idxArray[k] + "-" + Math.round(modEntry3.mass * 1000) + "-" + idxArray[l] + "-" + Math.round(modEntry4.mass * 1000));
                                                                        PositionDeltaMassMap positionDeltaMassMap = new PositionDeltaMassMap(ptmFreePeptide.length());
                                                                        positionDeltaMassMap.put(new Coordinate(idxArray[i], idxArray[i] + 1), modEntry1.mass);
                                                                        positionDeltaMassMap.put(new Coordinate(idxArray[j], idxArray[j] + 1), modEntry2.mass);
                                                                        positionDeltaMassMap.put(new Coordinate(idxArray[k], idxArray[k] + 1), modEntry3.mass);
                                                                        positionDeltaMassMap.put(new Coordinate(idxArray[l], idxArray[l] + 1), modEntry4.mass);
                                                                        Peptide peptide = new Peptide(ptmFreePeptide, isDecoy, massTool, localMaxMS2Charge, normalizedCrossCorr, globalRank);
                                                                        peptide.setVarPTM(positionDeltaMassMap);
                                                                        double score = massTool.buildVectorAndCalXCorr(peptide.getIonMatrix(), precursorCharge, expProcessedPL);
                                                                        if (score > 0) {
                                                                            peptide.setScore(score);
                                                                            peptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, localMaxMS2Charge, peptide.getIonMatrix(), ms2Tolerance));
                                                                            peptidePTMPattern.push(peptide);
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    private void try5PTMs(Map<Integer, Set<VarModParam>> idxVarModMap, double leftMassBound, double rightMassBound, String ptmFreePeptide, boolean isDecoy, double normalizedCrossCorr, int globalRank, Set<String> checkedPtmPattern, PeptidePTMPattern peptidePTMPattern, SparseVector expProcessedPL, TreeMap<Double, Double> plMap, int precursorCharge, int localMaxMS2Charge) { // only allow one low priority PTM
        Integer[] idxArray = idxVarModMap.keySet().toArray(new Integer[0]);
        Arrays.sort(idxArray);
        for (int i = 0; i < idxArray.length - 4; ++i) {
            for (VarModParam modEntry1 : idxVarModMap.get(idxArray[i])) {
                for (int j = i + 1; j < idxArray.length - 3; ++j) {
                    for (VarModParam modEntry2 : idxVarModMap.get(idxArray[j])) {
                        if (modEntry1.priority + modEntry2.priority > 0) {
                            if (Math.abs(modEntry1.mass + modEntry2.mass) >= ptmMassTolerance) {
                                for (int k = j + 1; k < idxArray.length - 2; ++k) {
                                    for (VarModParam modEntry3 : idxVarModMap.get(idxArray[k])) {
                                        if (modEntry1.priority + modEntry2.priority + modEntry3.priority > 1) {
                                            if (Math.abs(modEntry1.mass + modEntry3.mass) >= ptmMassTolerance && Math.abs(modEntry2.mass + modEntry3.mass) >= ptmMassTolerance) {
                                                for (int l = k + 1; l < idxArray.length - 1; ++l) {
                                                    for (VarModParam modEntry4 : idxVarModMap.get(idxArray[l])) {
                                                        if (modEntry1.priority + modEntry2.priority + modEntry3.priority + modEntry4.priority > 2) {
                                                            if (Math.abs(modEntry1.mass + modEntry4.mass) >= ptmMassTolerance && Math.abs(modEntry2.mass + modEntry4.mass) >= ptmMassTolerance && Math.abs(modEntry3.mass + modEntry4.mass) >= ptmMassTolerance) {
                                                                for (int m = l + 1; m < idxArray.length; ++m) {
                                                                    for (VarModParam modEntry5 : idxVarModMap.get(idxArray[m])) {
                                                                        if (modEntry1.priority + modEntry2.priority + modEntry3.priority + modEntry4.priority + modEntry5.priority > 3) {
                                                                            if (Math.abs(modEntry1.mass + modEntry5.mass) >= ptmMassTolerance && Math.abs(modEntry2.mass + modEntry5.mass) >= ptmMassTolerance && Math.abs(modEntry3.mass + modEntry5.mass) >= ptmMassTolerance && Math.abs(modEntry4.mass + modEntry5.mass) >= ptmMassTolerance) {
                                                                                if (modEntry1.mass + modEntry2.mass + modEntry3.mass + modEntry4.mass + modEntry5.mass <= rightMassBound && modEntry1.mass + modEntry2.mass + modEntry3.mass + modEntry4.mass + modEntry5.mass >= leftMassBound) {
                                                                                    if (!checkedPtmPattern.contains(idxArray[i] + "-" + Math.round(modEntry1.mass * 1000) + "-" + idxArray[j] + "-" + Math.round(modEntry2.mass * 1000) + "-" + idxArray[k] + "-" + Math.round(modEntry3.mass * 1000) + "-" + idxArray[l] + "-" + Math.round(modEntry4.mass * 1000) + "-" + idxArray[m] + "-" + Math.round(modEntry5.mass * 1000))) {
                                                                                        checkedPtmPattern.add(idxArray[i] + "-" + Math.round(modEntry1.mass * 1000) + "-" + idxArray[j] + "-" + Math.round(modEntry2.mass * 1000) + "-" + idxArray[k] + "-" + Math.round(modEntry3.mass * 1000) + "-" + idxArray[l] + "-" + Math.round(modEntry4.mass * 1000) + "-" + idxArray[m] + "-" + Math.round(modEntry5.mass * 1000));
                                                                                        PositionDeltaMassMap positionDeltaMassMap = new PositionDeltaMassMap(ptmFreePeptide.length());
                                                                                        positionDeltaMassMap.put(new Coordinate(idxArray[i], idxArray[i] + 1), modEntry1.mass);
                                                                                        positionDeltaMassMap.put(new Coordinate(idxArray[j], idxArray[j] + 1), modEntry2.mass);
                                                                                        positionDeltaMassMap.put(new Coordinate(idxArray[k], idxArray[k] + 1), modEntry3.mass);
                                                                                        positionDeltaMassMap.put(new Coordinate(idxArray[l], idxArray[l] + 1), modEntry4.mass);
                                                                                        positionDeltaMassMap.put(new Coordinate(idxArray[m], idxArray[m] + 1), modEntry5.mass);
                                                                                        Peptide peptide = new Peptide(ptmFreePeptide, isDecoy, massTool, localMaxMS2Charge, normalizedCrossCorr, globalRank);
                                                                                        peptide.setVarPTM(positionDeltaMassMap);
                                                                                        double score = massTool.buildVectorAndCalXCorr(peptide.getIonMatrix(), precursorCharge, expProcessedPL);
                                                                                        if (score > 0) {
                                                                                            peptide.setScore(score);
                                                                                            peptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, localMaxMS2Charge, peptide.getIonMatrix(), ms2Tolerance));
                                                                                            peptidePTMPattern.push(peptide);
                                                                                        }
                                                                                    }
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
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
