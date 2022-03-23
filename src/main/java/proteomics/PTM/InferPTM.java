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
import gurobi.*;
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
//        Map<Character, Set<VarModParam>> ptmMap = readModFile();
        Map<Character, Set<VarModParam>> ptmMap = new HashMap<>();

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
            PTMs(idxVarModMap, leftMassBound, rightMassBound, ptmFreePeptide, isDecoy, normalizedCrossCorr, globalRank, checkedPtmPattern2, peptidePTMPattern, expProcessedPL, plMap, precursorCharge, localMaxMS2Charge);
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

    public PeptidePTMPattern findPTM(GRBEnv env, int scanNum, SparseVector expProcessedPL, TreeMap<Double, Double> plMap, double precursorMass, String ptmFreePeptide, boolean isDecoy, double normalizedCrossCorr, char leftFlank, char rightFlank, int globalRank, int precursorCharge, int localMaxMS2Charge, double localMS1ToleranceL, double localMS1ToleranceR, List<ThreeExpAA> expAaLists) {
        double ptmFreeMass = massTool.calResidueMass(ptmFreePeptide) + massTool.H2O;

        TreeMap<Double, Double> unUsedPlMap = new TreeMap<>(plMap);
        double maxLMz = 0;
        double minRMz = plMap.lastKey();
        PeptidePTMPattern allPtmPattern = new PeptidePTMPattern(ptmFreePeptide);

        double totalDeltaMass = precursorMass - ptmFreeMass;
        if (Math.abs(totalDeltaMass) < 0.01 ) {
            return allPtmPattern;
        }
        Set<Integer> fixModIdxes = getFixModIdxes(ptmFreePeptide, fixModMap);  // record the indexes of the amino acids with fix mod, e.g. when there is no C, this is empty set.

        Peptide peptide = new Peptide(ptmFreePeptide, isDecoy, massTool, localMaxMS2Charge, normalizedCrossCorr, globalRank);
        double[][] theoIonsMatrix = peptide.getTheoIonMatrix();// Note that if you set varMod on it, it returns b y ions with var mod
        String ptmFreePeptideOrdinary = ptmFreePeptide.replaceAll("I","#").replaceAll("L","#");

        Map<Integer, Set<VarModParam>> idxVarModMap = getIdxVarModMap(ptmFreePeptide, fixModIdxes, leftFlank, rightFlank);
        String aasWithNoPtm = "CIKLPR";
        // idxVarModMap records: in this peptide, each amino acid could bear which set of modifications
        Map<Integer, VarModParam[]> idxVarModArrayMap = new HashMap<>();
        for (int id : idxVarModMap.keySet()){
            if (aasWithNoPtm.contains(ptmFreePeptide.substring(id, id+1))){
                continue;
            }
            VarModParam[] modArray = new VarModParam[idxVarModMap.get(id).size()];
            idxVarModMap.get(id).toArray(modArray);
            Arrays.sort(modArray, Comparator.comparingDouble(VarModParam::getMass));
            idxVarModArrayMap.put(id, modArray);
        }
        Set<Integer> modifiedZone = IntStream.range(1, ptmFreePeptide.length()-1).boxed().collect(Collectors.toSet());// e.g. nABCDEFc, modifiedZone={1,2,3,4,5,6}

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
                    for (double mz : tagInfo.mzs) {
                        unUsedPlMap.remove(mz);
                    }
                    if (tagInfo.mzs[3] > maxLMz) maxLMz = tagInfo.mzs[3];
                    break;
                } else if (Math.abs(deltaMass - totalDeltaMass) <= 0.01 && alignBPos > Collections.min(modifiedZone)) {
                    //skip if it has ptm but it is at the n ter
                    if (alignBPos < minHeadPosOfMmodTags) {
                        minHeadPosOfMmodTags = alignBPos;
                    }
                    tagZoneScore += tagInfo.getTotalIntensity();
                    tagInfo.isGoodTag3 = true;
                    seqOfGoodTags.add(tagSeq);
                    for (double mz : tagInfo.mzs) {
                        unUsedPlMap.remove(mz);
                    }
                    if (tagInfo.mzs[0] < minRMz) minRMz = tagInfo.mzs[0];
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
            System.out.println(scanNum + " is empty modifiedZone before tag 2");
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
                    for (double mz : tagInfo.mzs) {
                        if (mz == tagInfo.mzs[3]) break;
                        unUsedPlMap.remove(mz);
                    }
                    if (tagInfo.mzs[2] > maxLMz) maxLMz = tagInfo.mzs[2];
                } else if(Math.abs(deltaMass - totalDeltaMass) <= 0.01){
                    if (alignPos <= Collections.min(modifiedZone)) {
                        continue;
                    }
                    for (int i = alignPos; i < ptmFreePeptide.length()-1; i++){
                        modifiedZone.remove(i);
                    }
                    for (double mz : tagInfo.mzs) {
                        if (mz == tagInfo.mzs[3]) break;
                        unUsedPlMap.remove(mz);
                    }
                    if (tagInfo.mzs[0] < minRMz) minRMz = tagInfo.mzs[0];

                }
            }
        }


        if (modifiedZone.size() == 0) {
            System.out.println(scanNum + " is empty modifiedZone after tag 2");
            return allPtmPattern; //Some scans are not valid Scans. Will be deleted soon.
        }

        // find the most frequently appearing extra mass as proof for multiple PTMs.
        double top1ExtraMass = -9999d;
        if (extraMassTagIdMap.size() != 0) {
            Map<Double, Integer> extraMassSizeMap = new HashMap<>(extraMassTagIdMap.entrySet().size());
            for (Map.Entry<Double, Set<Integer>> entry : extraMassTagIdMap.entrySet()) {
                extraMassSizeMap.put(entry.getKey(), entry.getValue().size());
            }
            List<Map.Entry<Double, Integer>> tempList = new ArrayList<>(extraMassSizeMap.entrySet());
            Collections.sort(tempList, Comparator.comparingInt(Map.Entry<Double, Integer>::getValue)); //get the key of the largest value in extraMassMap
            if (tempList.get(tempList.size()-1).getValue() > 1){
                top1ExtraMass = tempList.get(tempList.size()-1).getKey();
            }
        }

        Set<Integer> extraMassTagZone = new HashSet<>();
        Set<Integer> constraintZone = new HashSet<>(idxVarModArrayMap.keySet());// e.g. nABCDEFc, modifiedZone={1,2,3,4,5,6}

        List<Integer> extraTagIds = new ArrayList<>();
        mainLoop:
        if (top1ExtraMass != -9999d) {
            for (int tagId : extraMassTagIdMap.get(top1ExtraMass)){
                ThreeExpAA tagInfo = expAaLists.get(tagId);
                for (int alignPos : tagInfo.bAlignPosMassMap.keySet()) {
                    if( top1ExtraMass - Math.abs(tagInfo.bAlignPosMassMap.get(alignPos)) <= 0.01){
                        extraTagIds.add(tagId);
                        constraintZone.retainAll( IntStream.range(1, alignPos).boxed().collect(Collectors.toSet()) );
                        for (int i = alignPos; i < alignPos+3; i++) {
                            extraMassTagZone.add(i);
                        }
                        break;
                    }
                }
            }
        }

        constraintZone.retainAll(modifiedZone);  //only consider one extraMass now

        if (maxLMz > minRMz){
            int a = 1;
        }

        List<Double> mzToRemove = new ArrayList<>();
        for (double mz : unUsedPlMap.keySet()) {
            if (unUsedPlMap.get(mz) < 0.0 ) {
                mzToRemove.add(mz);
            }
        }

        long start = System.currentTimeMillis();
        for (double mz : mzToRemove){
            unUsedPlMap.remove(mz);
        }

//        positionOneMap.clear();
//        Map<Integer, Integer> positionOneMap = new HashMap<>();
//        positionOneMap.put(3,0);
//        positionOneMap.put(16,1);
//        positionOneMap.put(17,9);
//        PositionDeltaMassMap positionDeltaMassMap1 = new PositionDeltaMassMap(ptmFreePeptide.length());
//        Peptide p1p = new Peptide(ptmFreePeptide, isDecoy, massTool, localMaxMS2Charge, normalizedCrossCorr, globalRank);
//        for (Map.Entry<Integer, Integer> entry : positionOneMap.entrySet()){
//            positionDeltaMassMap1.put(new Coordinate(entry.getKey(), entry.getKey() + 1), idxVarModArrayMap.get(entry.getKey())[entry.getValue()].mass);
//        }
//        p1p.setVarPTM(positionDeltaMassMap1);
//        double[][] temp = p1p.getIonMatrix();
//        double ss = massTool.buildVectorAndCalXCorr(p1p.getIonMatrix(), precursorCharge, expProcessedPL) ;
        //stub end

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
        double cleanScore = massTool.buildVectorAndCalXCorr(cleanPep.getIonMatrixNow(), 1, expProcessedPL, matchedBions, matchedYions, jRange) ;
        cleanPep.setScore(cleanScore);
        cleanPep.setMatchedPeakNum(Score.getMatchedIonNum(plMap, localMaxMS2Charge, cleanPep.getIonMatrix(), ms2Tolerance));
        cleanPep.matchedBions.putAll(matchedBions);
        cleanPep.matchedYions.putAll(matchedYions);
        if (!matchedBions.isEmpty()) {
            lb = Collections.max(matchedBions.keySet()) + 1;
        }
        if (!matchedYions.isEmpty()) {
            rb = Collections.min(matchedYions.keySet());
        }
        modifiedZone = IntStream.range(lb, rb).boxed().collect(Collectors.toSet());

        PeptidePTMPattern ptmZ1Good = new PeptidePTMPattern(ptmFreePeptide, 1);
        PeptidePTMPattern ptmZ1Bad = new PeptidePTMPattern(ptmFreePeptide, 1);
        PeptidePTMPattern ptmInitialTemp = new PeptidePTMPattern(ptmFreePeptide, 1);
        Peptide cleanPeptide = new Peptide(ptmFreePeptide, isDecoy, massTool, localMaxMS2Charge, normalizedCrossCorr, globalRank);
        cleanPeptide.setVarPTM(new PositionDeltaMassMap(ptmFreePeptide.length()));
        ptmInitialTemp.push(cleanPeptide);
//        long t1 = System.currentTimeMillis();
        if (modifiedZone.isEmpty()) {
            return allPtmPattern;
        }
        DividedZone z1Z2Res = dividePep(scanNum, modifiedZone, ptmZ1Good, ptmZ1Bad, ptmInitialTemp, idxVarModArrayMap, totalDeltaMass, ptmFreePeptide, isDecoy, normalizedCrossCorr, globalRank, expProcessedPL, plMap, precursorCharge, localMaxMS2Charge, unUsedPlMap);
//        long t2 = System.currentTimeMillis();
//        System.out.println("End dividePep, "+ (t2-t1));

        if (!ptmZ1Good.getPeptideTreeSet().isEmpty()) {
            allPtmPattern.push(ptmZ1Good.getTopPepPtn());
        }
        Set<Integer> zone1 = z1Z2Res.keptZone;
        Set<Integer> zone2 = z1Z2Res.freeZone;
        double zone1Mass = ptmZ1Bad.getTopPepPtn().getVarPTMs().values().iterator().next();
        double zone2Mass = totalDeltaMass - zone1Mass;

        //for zone2 if numPtmZone1 == 2

        PeptidePTMPattern ptmZ3Good = new PeptidePTMPattern(ptmFreePeptide, 1);
        PeptidePTMPattern ptmZ3Bad = new PeptidePTMPattern(ptmFreePeptide, 1);

        if (zone2.isEmpty()) {
            return allPtmPattern;
        }
        DividedZone z3Z4Res = dividePep(scanNum, zone2, ptmZ3Good, ptmZ3Bad ,ptmZ1Bad, idxVarModArrayMap, zone2Mass, ptmFreePeptide, isDecoy, normalizedCrossCorr, globalRank, expProcessedPL, plMap, precursorCharge, localMaxMS2Charge, unUsedPlMap);
        if (!ptmZ3Good.getPeptideTreeSet().isEmpty()) {
            allPtmPattern.push(ptmZ3Good.getTopPepPtn());
        }
        Set<Integer> zone3 = z3Z4Res.keptZone;
        Set<Integer> zone4 = z3Z4Res.freeZone;
        double zone3Mass = z3Z4Res.ptmMass;
        double zone4Mass = zone2Mass - zone3Mass;

        //solve zone 4
        PeptidePTMPattern ptmZ5Good = new PeptidePTMPattern(ptmFreePeptide, 1);
        PeptidePTMPattern ptmZ5Bad = new PeptidePTMPattern(ptmFreePeptide, 1);

        if (!zone4.isEmpty()) {
            DividedZone z5Z6Res = dividePep(scanNum, zone4, ptmZ5Good, ptmZ5Bad, ptmZ3Bad, idxVarModArrayMap, zone4Mass, ptmFreePeptide, isDecoy, normalizedCrossCorr, globalRank, expProcessedPL, plMap, precursorCharge, localMaxMS2Charge, unUsedPlMap);
        }
        if (!ptmZ5Good.peptideTreeSet.isEmpty()) {
            allPtmPattern.push(ptmZ5Good.getTopPepPtn());
        }

        long end = System.currentTimeMillis();

//        System.out.println("N1 time, "+ (end - start));
        for (Peptide pep : allPtmPattern.peptideTreeSet) {
            double score = massTool.buildVectorAndCalXCorr(pep.getIonMatrixNow(), 1, expProcessedPL);
            pep.setScore(score);
            peptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, localMaxMS2Charge, peptide.getIonMatrixNow(), ms2Tolerance));

        }
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
                if (ptmId == 0 && pos == 7){
                    int a = 1;
                }
                Set<Integer> jRange = new HashSet<>();
                for (int i : modifiedZone) {
                    jRange.add(i-1);
                }
                double score = massTool.buildVectorAndCalXCorr(peptide.getIonMatrixNow(), 1, expProcessedPL, matchedBions, matchedYions, jRange) ;
                if (score > -1) {
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
            System.out.println(scanNum);
            int a= 1;
            return new DividedZone(new HashSet<>(), new HashSet<>(), ptmMass);
        }
        Peptide inFeasiTopPep = ptmBadRes.getTopPepPtn();
        Set<Integer> keptZone = new HashSet<>();
        Set<Integer> freeZone = new HashSet<>();
//        double bMassNew = 0;
//        double yMassNew = 0;
        int n1Pos = -1;
        for (Coordinate coor : inFeasiTopPep.getVarPTMs().keySet()) {
            if (modifiedZone.contains(coor.x)) {
                n1Pos = coor.x;
                ptmMass = inFeasiTopPep.getVarPTMs().get(coor);
            }
        }
//        int n1Pos = inFeasiTopPep.getVarPTMs().keySet().iterator().next().x-offset;
        if (inFeasiTopPep.matchedYions.size() > inFeasiTopPep.matchedBions.size()) { // Y >>> B
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

    private int collectPtmCombs(int scanNum, GRBEnv env, double[] bIons1,double[] yIons1, int numPtmsOnPep, Set<Integer> modifiedZone, Map<Integer, VarModParam[]> idxVarModMap, double totalDeltaMass, String ptmFreePeptide, boolean isDecoy, double normalizedCrossCorr, int globalRank, PeptidePTMPattern peptidePTMPattern, Peptide templatePep, SparseVector expProcessedPL, TreeMap<Double, Double> plMap, int precursorCharge, int localMaxMS2Charge, Map<Double, Double> unUsedPlMap) { // Sometimes, the precursor mass error may affects the digitized spectrum.
        List<Double> matchedPeaks = new ArrayList<>();
        double[][] theo = templatePep.getIonMatrixNow();
        double[] bIons = theo[0];
        double[] yIons = theo[1];
        try {
            GRBModel model = new GRBModel(env);
            double t = 0.01;

            //variables
            Map<Integer, GRBVar[]> posVarsMap = new HashMap<>();
            GRBLinExpr totalPtmsOnPepConstr = new GRBLinExpr();
            GRBLinExpr totalMassOnPepConstr = new GRBLinExpr();
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
                }
                totalMassOnPepConstr.addTerms(coeffMassAaArray, posVarsMap.get(pos));

            }

            //constraints
            model.addConstr(totalPtmsOnPepConstr, GRB.EQUAL, numPtmsOnPep, "CON_totalNumPtmsOnPep");

            model.addConstr(totalMassOnPepConstr, GRB.GREATER_EQUAL, totalDeltaMass - t, "CON_totalMassOnPep_GEQ");
            model.addConstr(totalMassOnPepConstr, GRB.LESS_EQUAL, totalDeltaMass + t, "CON_totalMassOnPep_LEQ"); //or put this to constraints as a model.addRange


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
                    model.addConstr(bVarMass, GRB.LESS_EQUAL, Drec*(bIons[pos-1]-peakList.get(pId).getKey())+tau, "CON_b_pos"+pos+ptmFreePeptide.charAt(pos)+"_peak"+pId+"_LEQ1");
                    model.addConstr(bVarMass2, GRB.LESS_EQUAL, Drec*(peakList.get(pId).getKey()-bIons[pos-1])+tau, "CON_b_pos"+pos+ptmFreePeptide.charAt(pos)+"_peak"+pId+"_LEQ2");
                }
                bPosPeakIdMap.put(pos, peakIdList);
                bPosVars.put(pos, varsMap);
                model.addConstr(bVarSumAtPos, GRB.LESS_EQUAL, 1, "CON_b_pos"+pos+ptmFreePeptide.charAt(pos)+"_max1Peak");
            }

            for (int pos1 : modifiedZone) {   //for No cross constraints
                if (pos1 == Collections.max(modifiedZone)) continue;
                int pos2 = pos1 + 1;
                for (int pId1 : bPosPeakIdMap.get(pos1)) {
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

            //obj
            model.setObjective(numMatchedPeaks, GRB.MAXIMIZE);


            //settings
//            model.set(GRB.IntParam.MIPFocus, 1);
            // solve the model
//            model.write("/home/slaiad/Data/Simulation_Data/test.lp");
//            model.write("/home/slaiad/Data/Simulation_Data/test.mps");
            int numOfSols = 0;
            int maxLastNPeaks = 0;

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
                    return 0;
            }

            Map<Integer, Integer> positionOneMap = new HashMap<>();
            for (int pos : modifiedZone) {
                if (!posVarsMap.containsKey(pos)) continue;
                GRBVar[] varsResArray = posVarsMap.get(pos);
                for (int varId = 0; varId < varsResArray.length; varId++ ) {  //should I be careful about the binary variable to be really 0/1 instead of decimal
                    if (1 == Math.round(varsResArray[varId].get(GRB.DoubleAttr.X))) {
                        positionOneMap.put(pos, varId);
                    }
                }
            }

            //stub test true other
//                positionOneMap.clear();
//                positionOneMap.put(1,5);
//                positionOneMap.put(10,13);
//                positionOneMap.put(12,1);
            //stub end

            PositionDeltaMassMap templatePtm = new PositionDeltaMassMap(ptmFreePeptide.length());
            Peptide peptide = new Peptide(ptmFreePeptide, isDecoy, massTool, localMaxMS2Charge, normalizedCrossCorr, globalRank);
            if (templatePep.hasVarPTM()) {
                for (Coordinate coor : templatePep.getVarPTMs().keySet()){
                    templatePtm.put(new Coordinate(coor.x, coor.y),  templatePep.getVarPTMs().get(coor));
                }
            }
            peptide.setVarPTM(templatePtm);
//            Peptide peptide = peptidePTMPattern.getTopPepPtn();
            PositionDeltaMassMap positionDeltaMassMap = new PositionDeltaMassMap(ptmFreePeptide.length());
            for (Map.Entry<Integer, Integer> entry : positionOneMap.entrySet()){
                positionDeltaMassMap.put(new Coordinate(entry.getKey(), entry.getKey() + 1), idxVarModMap.get(entry.getKey())[entry.getValue()].mass);
            }
            peptide.addVarPTM(positionDeltaMassMap);

            double[][] temp = peptide.getIonMatrixNow();
            double score = massTool.buildVectorAndCalXCorr(peptide.getIonMatrixNow(), 1, expProcessedPL) - 0.01*(numPtmsOnPep-1);
            double lastLpScore = 0;
            int currentNPeaks = 0;
            if (score > -11) {
                peptide.lpScore = model.get(GRB.DoubleAttr.ObjVal);
                System.out.println("numPtmsOnPep "+numPtmsOnPep+" lpScore "+peptide.lpScore + " XCorr " + score);
                peptide.setScore(score);
                peptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, localMaxMS2Charge, peptide.getIonMatrix(), ms2Tolerance));
                if (peptidePTMPattern.peptideTreeSet.isEmpty()){
                    peptidePTMPattern.push(peptide);
                } else if(peptide.getScore() > peptidePTMPattern.getTopPepPtn().getScore()) {
                    peptidePTMPattern.push(peptide);
                }
                for (int pos : bPosVars.keySet()) {
                    for (int peakId : bPosVars.get(pos).keySet()){
                        if (1 == Math.round(bPosVars.get(pos).get(peakId).get(GRB.DoubleAttr.X))){
                            matchedPeaks.add(peakList.get(peakId).getKey());
                            System.out.println("bMatch "+bPosVars.get(pos).get(peakId).get(GRB.StringAttr.VarName) + ", "+ peakList.get(peakId).getKey()+","+peakList.get(peakId).getValue());
                        }
                    }
                }
                for (int pos : yPosVars.keySet()) {
                    for (int peakId : yPosVars.get(pos).keySet()){
                        if (1 == Math.round(yPosVars.get(pos).get(peakId).get(GRB.DoubleAttr.X))){
                            matchedPeaks.add(peakList.get(peakId).getKey());

                            System.out.println("yMatch "+yPosVars.get(pos).get(peakId).get(GRB.StringAttr.VarName) + ", "+ peakList.get(peakId).getKey()+","+peakList.get(peakId).getValue());
                        }
                    }
                }
                System.out.println("this  "+ lastLpScore + "," + currentNPeaks + "," +peptide);

//                List<Double> massList = new ArrayList<>();
//                for (double mass : positionDeltaMassMap.values()){
//                    massList.add(mass);
//                }
//                Collections.sort(massList);
//                    System.out.println(numOfSols + ", " + score + ", "+peptide.getMatchedPeakNum()+", "+ massList);
            }

            model.dispose();
            if (peptidePTMPattern.peptideTreeSet.isEmpty()){
                int a = 1;
            }
            return 1;
        } catch (GRBException e) {
            System.out.println("Error code: " + e.getErrorCode() + ". " + e.getMessage());
        }
        return 0;
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



    private void PTMs(Map<Integer, Set<VarModParam>> idxVarModMap, double leftMassBound, double rightMassBound, String ptmFreePeptide, boolean isDecoy, double normalizedCrossCorr, int globalRank, Set<String> checkedPtmPattern, PeptidePTMPattern peptidePTMPattern, SparseVector expProcessedPL, TreeMap<Double, Double> plMap, int precursorCharge, int localMaxMS2Charge) {
        Integer[] idxArray = idxVarModMap.keySet().toArray(new Integer[0]);
        Arrays.sort(idxArray);
        for (int i = 0; i < idxArray.length - 1; ++i) {
            for (VarModParam modEntry1 : idxVarModMap.get(idxArray[i])) {
                for (int j = i + 1; j < idxArray.length; ++j) {
                    for (VarModParam modEn : idxVarModMap.get(idxArray[j])) {
                        if (modEntry1.mass + modEn.mass <= rightMassBound && modEntry1.mass + modEn.mass >= leftMassBound) {
                            if (!checkedPtmPattern.contains(idxArray[i] + "-" + Math.round(modEntry1.mass * 1000) + "-" + idxArray[j] + "-" + Math.round(modEn.mass * 1000))) {
                                checkedPtmPattern.add(idxArray[i] + "-" + Math.round(modEntry1.mass * 1000) + "-" + idxArray[j] + "-" + Math.round(modEn.mass * 1000));
                                PositionDeltaMassMap positionDeltaMassMap = new PositionDeltaMassMap(ptmFreePeptide.length());
                                positionDeltaMassMap.put(new Coordinate(idxArray[i], idxArray[i] + 1), modEntry1.mass);
                                positionDeltaMassMap.put(new Coordinate(idxArray[j], idxArray[j] + 1), modEn.mass);
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
                    for (VarModParam modEn : idxVarModMap.get(idxArray[j])) {
                        if (modEntry1.priority + modEn.priority > 0) {
                            if (Math.abs(modEntry1.mass + modEn.mass) >= ptmMassTolerance) { // two self cancelled PTM masses are not allowed.
                                for (int k = j + 1; k < idxArray.length; ++k) {
                                    for (VarModParam modEntry3 : idxVarModMap.get(idxArray[k])) {
                                        if (modEntry1.priority + modEn.priority + modEntry3.priority > 1) {
                                            if (Math.abs(modEntry1.mass + modEntry3.mass) >= ptmMassTolerance && Math.abs(modEn.mass + modEntry3.mass) >= ptmMassTolerance) {
                                                if (modEntry1.mass + modEn.mass + modEntry3.mass <= rightMassBound && modEntry1.mass + modEn.mass + modEntry3.mass >= leftMassBound) {
                                                    if (!checkedPtmPattern.contains(idxArray[i] + "-" + Math.round(modEntry1.mass * 1000) + "-" + idxArray[j] + "-" + Math.round(modEn.mass * 1000) + "-" + idxArray[k] + "-" + Math.round(modEntry3.mass * 1000))) {
                                                        checkedPtmPattern.add(idxArray[i] + "-" + Math.round(modEntry1.mass * 1000) + "-" + idxArray[j] + "-" + Math.round(modEn.mass * 1000) + "-" + idxArray[k] + "-" + Math.round(modEntry3.mass * 1000));
                                                        PositionDeltaMassMap positionDeltaMassMap = new PositionDeltaMassMap(ptmFreePeptide.length());
                                                        positionDeltaMassMap.put(new Coordinate(idxArray[i], idxArray[i] + 1), modEntry1.mass);
                                                        positionDeltaMassMap.put(new Coordinate(idxArray[j], idxArray[j] + 1), modEn.mass);
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
                    for (VarModParam modEn : idxVarModMap.get(idxArray[j])) {
                        if (modEntry1.priority + modEn.priority > 0) {
                            if (Math.abs(modEntry1.mass + modEn.mass) >= ptmMassTolerance) {
                                for (int k = j + 1; k < idxArray.length - 1; ++k) {
                                    for (VarModParam modEntry3 : idxVarModMap.get(idxArray[k])) {
                                        if (modEntry1.priority + modEn.priority + modEntry3.priority > 1) {
                                            if (Math.abs(modEntry1.mass + modEntry3.mass) >= ptmMassTolerance && Math.abs(modEn.mass + modEntry3.mass) >= ptmMassTolerance) {
                                                for (int l = k + 1; l < idxArray.length; ++l) {
                                                    for (VarModParam modEntry4 : idxVarModMap.get(idxArray[l])) {
                                                        if (modEntry1.priority + modEn.priority + modEntry3.priority + modEntry4.priority > 2) {
                                                            if (Math.abs(modEntry1.mass + modEntry4.mass) >= ptmMassTolerance && Math.abs(modEn.mass + modEntry4.mass) >= ptmMassTolerance && Math.abs(modEntry3.mass + modEntry4.mass) >= ptmMassTolerance) {
                                                                if (modEntry1.mass + modEn.mass + modEntry3.mass + modEntry4.mass <= rightMassBound && modEntry1.mass + modEn.mass + modEntry3.mass + modEntry4.mass >= leftMassBound) {
                                                                    if (!checkedPtmPattern.contains(idxArray[i] + "-" + Math.round(modEntry1.mass * 1000) + "-" + idxArray[j] + "-" + Math.round(modEn.mass * 1000) + "-" + idxArray[k] + "-" + Math.round(modEntry3.mass * 1000) + "-" + idxArray[l] + "-" + Math.round(modEntry4.mass * 1000))) {
                                                                        checkedPtmPattern.add(idxArray[i] + "-" + Math.round(modEntry1.mass * 1000) + "-" + idxArray[j] + "-" + Math.round(modEn.mass * 1000) + "-" + idxArray[k] + "-" + Math.round(modEntry3.mass * 1000) + "-" + idxArray[l] + "-" + Math.round(modEntry4.mass * 1000));
                                                                        PositionDeltaMassMap positionDeltaMassMap = new PositionDeltaMassMap(ptmFreePeptide.length());
                                                                        positionDeltaMassMap.put(new Coordinate(idxArray[i], idxArray[i] + 1), modEntry1.mass);
                                                                        positionDeltaMassMap.put(new Coordinate(idxArray[j], idxArray[j] + 1), modEn.mass);
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
                    for (VarModParam modEn : idxVarModMap.get(idxArray[j])) {
                        if (modEntry1.priority + modEn.priority > 0) {
                            if (Math.abs(modEntry1.mass + modEn.mass) >= ptmMassTolerance) {
                                for (int k = j + 1; k < idxArray.length - 2; ++k) {
                                    for (VarModParam modEntry3 : idxVarModMap.get(idxArray[k])) {
                                        if (modEntry1.priority + modEn.priority + modEntry3.priority > 1) {
                                            if (Math.abs(modEntry1.mass + modEntry3.mass) >= ptmMassTolerance && Math.abs(modEn.mass + modEntry3.mass) >= ptmMassTolerance) {
                                                for (int l = k + 1; l < idxArray.length - 1; ++l) {
                                                    for (VarModParam modEntry4 : idxVarModMap.get(idxArray[l])) {
                                                        if (modEntry1.priority + modEn.priority + modEntry3.priority + modEntry4.priority > 2) {
                                                            if (Math.abs(modEntry1.mass + modEntry4.mass) >= ptmMassTolerance && Math.abs(modEn.mass + modEntry4.mass) >= ptmMassTolerance && Math.abs(modEntry3.mass + modEntry4.mass) >= ptmMassTolerance) {
                                                                for (int m = l + 1; m < idxArray.length; ++m) {
                                                                    for (VarModParam modEntry5 : idxVarModMap.get(idxArray[m])) {
                                                                        if (modEntry1.priority + modEn.priority + modEntry3.priority + modEntry4.priority + modEntry5.priority > 3) {
                                                                            if (Math.abs(modEntry1.mass + modEntry5.mass) >= ptmMassTolerance && Math.abs(modEn.mass + modEntry5.mass) >= ptmMassTolerance && Math.abs(modEntry3.mass + modEntry5.mass) >= ptmMassTolerance && Math.abs(modEntry4.mass + modEntry5.mass) >= ptmMassTolerance) {
                                                                                if (modEntry1.mass + modEn.mass + modEntry3.mass + modEntry4.mass + modEntry5.mass <= rightMassBound && modEntry1.mass + modEn.mass + modEntry3.mass + modEntry4.mass + modEntry5.mass >= leftMassBound) {
                                                                                    if (!checkedPtmPattern.contains(idxArray[i] + "-" + Math.round(modEntry1.mass * 1000) + "-" + idxArray[j] + "-" + Math.round(modEn.mass * 1000) + "-" + idxArray[k] + "-" + Math.round(modEntry3.mass * 1000) + "-" + idxArray[l] + "-" + Math.round(modEntry4.mass * 1000) + "-" + idxArray[m] + "-" + Math.round(modEntry5.mass * 1000))) {
                                                                                        checkedPtmPattern.add(idxArray[i] + "-" + Math.round(modEntry1.mass * 1000) + "-" + idxArray[j] + "-" + Math.round(modEn.mass * 1000) + "-" + idxArray[k] + "-" + Math.round(modEntry3.mass * 1000) + "-" + idxArray[l] + "-" + Math.round(modEntry4.mass * 1000) + "-" + idxArray[m] + "-" + Math.round(modEntry5.mass * 1000));
                                                                                        PositionDeltaMassMap positionDeltaMassMap = new PositionDeltaMassMap(ptmFreePeptide.length());
                                                                                        positionDeltaMassMap.put(new Coordinate(idxArray[i], idxArray[i] + 1), modEntry1.mass);
                                                                                        positionDeltaMassMap.put(new Coordinate(idxArray[j], idxArray[j] + 1), modEn.mass);
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
