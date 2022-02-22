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

    public PeptidePTMPattern findPTM(GRBEnv env, int scanNum, SparseVector expProcessedPL, TreeMap<Double, Double> plMap, double precursorMass, String ptmFreePeptide, boolean isDecoy, double normalizedCrossCorr, char leftFlank, char rightFlank, int globalRank, int precursorCharge, int localMaxMS2Charge, double localMS1ToleranceL, double localMS1ToleranceR, List<ThreeExpAA> expAaLists) {
        double ptmFreeMass = massTool.calResidueMass(ptmFreePeptide) + massTool.H2O;
        double totalDeltaMass = precursorMass - ptmFreeMass;
        double leftMassBound = totalDeltaMass + localMS1ToleranceL;
        double rightMassBound = totalDeltaMass + localMS1ToleranceR;
        Set<Integer> fixModIdxes = getFixModIdxes(ptmFreePeptide, fixModMap);  // record the indexes of the amino acids with fix mod, e.g. when there is no C, this is empty set.

        PeptidePTMPattern peptidePTMPattern = new PeptidePTMPattern(ptmFreePeptide);
        Peptide peptide = new Peptide(ptmFreePeptide, isDecoy, massTool, localMaxMS2Charge, normalizedCrossCorr, globalRank);
        double[][] theoIonsMatrix = peptide.getTheoIonMatrix();// Note that if you set varMod on it, it returns b y ions with var mod

        Map<Integer, Set<VarModParam>> idxVarModMap = getIdxVarModMap(ptmFreePeptide, fixModIdxes, leftFlank, rightFlank);
        // idxVarModMap records: in this peptide, each amino acid could bear which set of modifications
        Map<Integer, VarModParam[]> idxVarModArrayMap = new HashMap<>();
        for (int id : idxVarModMap.keySet()){
            VarModParam[] modArray = new VarModParam[idxVarModMap.get(id).size()];
            idxVarModMap.get(id).toArray(modArray);
            idxVarModArrayMap.put(id, modArray);
        }
        Set<Integer> modifiedZone = new HashSet<>(idxVarModArrayMap.keySet());// e.g. nABCDEFc, modifiedZone={1,2,3,4,5,6}
//        for (int i = 1; i < ptmFreePeptide.length()-1; i++) {
//            modifiedZone.add(i);
//        }
        Set<Integer> cleanZone = new HashSet<>();
        Set<Integer> tagZone = new HashSet<>();
        double cleanZoneScore = 0d;
        double tagZoneScore = 0d;

        Set<Integer> usedTagIdxes = new HashSet<>();
        Set<Integer> remainedTagIdxes = new HashSet<>(expAaLists.size());
        for (int i = 0; i < expAaLists.size(); i++) {
            remainedTagIdxes.add(i);
        }
        //find feasible zone, clean tags and M_modTags, only for b ions now, wait for Y ions
        for (int tagId : remainedTagIdxes){
            ThreeExpAA tagInfo = expAaLists.get(tagId);
            String tagSeq = tagInfo.getPtmFreeAAString();
            int alignPos = ptmFreePeptide.indexOf(tagSeq);
            if (alignPos == -1) {
                if (!ptmFreePeptide.contains(tagSeq.substring(0,2))){
                    usedTagIdxes.add(tagId);
                }
                continue; // pep does not contain the tag
            }
            int headPeakId = alignPos-1;
            
            double deltaMass = tagInfo.getHeadLocation() - theoIonsMatrix[0][headPeakId];
            if (Math.abs(deltaMass) <= 0.02){
                usedTagIdxes.add(tagId);
                if (alignPos+3 > Collections.max(modifiedZone)){
                    continue; //if the tag is at the C term, it can not be right so skip.
                }
                for (int i = 1; i < alignPos+3; i++){
                    cleanZone.add(i);
                }
                cleanZoneScore += tagInfo.getTotalIntensity();
            } else if(Math.abs(deltaMass - totalDeltaMass) <= 0.02){
                usedTagIdxes.add(tagId);
                if (alignPos <= Collections.min(modifiedZone)) { //skip if it has ptm but it is at the n term
                    continue;
                }
                for (int i = alignPos; i < ptmFreePeptide.length()-1; i++){
                    tagZone.add(i);
                }
                tagZoneScore += tagInfo.getTotalIntensity();
            }
        }

        remainedTagIdxes.removeAll(usedTagIdxes);
        int minPosInTagZone = (tagZone.size() == 0) ? ptmFreePeptide.length()-1 : Collections.min(tagZone);
        int maxPosInCleanZone = (cleanZone.size() == 0) ? 0 : Collections.max(cleanZone);

        if (minPosInTagZone - maxPosInCleanZone >= 2){ // not interfering between tagZone and cleanZone
            modifiedZone.removeAll(tagZone);
            modifiedZone.removeAll(cleanZone);
        } else {
            if (tagZoneScore > cleanZoneScore){
                modifiedZone.removeAll(tagZone);
            } else {
                modifiedZone.removeAll(cleanZone);
            }
        }


        //find feasible zone, clean tags 2 and M_modTags 2, only for b ions now wait for Y ions
        cleanZone.clear();
        tagZone.clear();

        for (int tagId : remainedTagIdxes){
            ThreeExpAA tagInfo = expAaLists.get(tagId);
            String tagSeq = tagInfo.getPtmFreeAAString();
            if (ptmFreePeptide.contains(tagSeq)){
                continue; //this is to skip the m_mod_tag3
            }
            String tagSeqPrefix = tagSeq.substring(0,2);
            int alignPos = ptmFreePeptide.indexOf(tagSeqPrefix);
            if (alignPos == -1) {
                continue; // pep does not contain the tag
            }
            int headPeakId = alignPos - 1;
            double deltaMass = tagInfo.getHeadLocation() - theoIonsMatrix[0][headPeakId];
            if (Math.abs(deltaMass) <= 0.02){
                usedTagIdxes.add(tagId);
                if (alignPos+2 > Collections.max(modifiedZone)){
                    continue; //if the tag is at the C term, it can not be right so skip.
                }
                for (int i = 1; i < alignPos+2; i++){
                    cleanZone.add(i);
                }
            } else if(Math.abs(deltaMass - totalDeltaMass) <= 0.02){
                usedTagIdxes.add(tagId);
                if (alignPos <= Collections.min(modifiedZone)) {
                    continue;
                }
                for (int i = alignPos; i < ptmFreePeptide.length()-1; i++){
                    tagZone.add(i);
                }
//                tagZoneScore += tagInfo.getTotalIntensity();
            }
        }

        remainedTagIdxes.removeAll(usedTagIdxes);
        modifiedZone.removeAll(tagZone);
        modifiedZone.removeAll(cleanZone);


        // find m-mod-tags as proof for multiple PTMs.
        Map<Double, Integer> extraMassMap = new HashMap<>(); // records extra small mass appearance time.
        Set<Integer> extraMassTagIdxes = new HashSet<>();
        for (int tagId : remainedTagIdxes) {
            ThreeExpAA tagInfo = expAaLists.get(tagId);
            String tagSeq = tagInfo.getPtmFreeAAString();
            int alignPos = ptmFreePeptide.indexOf(tagSeq);
            if (alignPos == -1) {
                continue; // pep does not contain the tag
            }
            int headPeakId = alignPos-1;

            double deltaMass = tagInfo.getHeadLocation() - theoIonsMatrix[0][headPeakId];
            if (alignPos <= Collections.min(modifiedZone) || alignPos+3 > Collections.max(modifiedZone)){
                continue; //skip when the tag violates the boundry of modifiedZone. Should I consider when tag has been
                          //judged as clean tags or M_mod_tags?
            }
            boolean isNewMass = true;
            for (double mass : extraMassMap.keySet()) {
                if (Math.abs(deltaMass - mass) <= 0.02){
                    int newMassTimes = 1+extraMassMap.get(mass);
                    double newAveMass = (deltaMass+mass*extraMassMap.get(mass)) / newMassTimes;
                    extraMassMap.remove(mass);
                    extraMassMap.put(newAveMass, newMassTimes);
                    isNewMass = false;
                    break;
                }
            }
            if (isNewMass) {
                extraMassMap.put(deltaMass, 1);
            }
            extraMassTagIdxes.add(tagId);
        }
        double top1ExtraMass = -9999d;
        if (extraMassMap.size() != 0) {
            List<Map.Entry<Double, Integer>> extraMassList = new ArrayList(extraMassMap.entrySet());
            Collections.sort(extraMassList, Comparator.comparingInt(Map.Entry::getValue)); //get the key of the largest value in extraMassMap
            if (extraMassList.get(extraMassList.size()-1).getValue() > 1){
                top1ExtraMass = extraMassList.get(extraMassList.size()-1).getKey();
            }
        }
        Set<Integer> extraMassTagZone = new HashSet<>();
        Set<Integer> constraintZone = new HashSet<>(idxVarModArrayMap.keySet());// e.g. nABCDEFc, modifiedZone={1,2,3,4,5,6}
//        for (int i = 1; i < ptmFreePeptide.length()-1; i++) {
//            constraintZone.add(i);
//        }

        if (top1ExtraMass != -9999d) {
            for (int tagId : extraMassTagIdxes){
                ThreeExpAA tagInfo = expAaLists.get(tagId);
                String tagSeq = tagInfo.getPtmFreeAAString();
                int alignPos = ptmFreePeptide.indexOf(tagSeq);
                int headPeakId = alignPos - 1;
                double deltaMass = tagInfo.getHeadLocation() - theoIonsMatrix[0][headPeakId];
                if (Math.abs(deltaMass - top1ExtraMass) <= 0.02){
                    constraintZone.retainAll( IntStream.range(1, alignPos).boxed().collect(Collectors.toSet()) );
                    for (int i = alignPos; i < alignPos+3; i++) {
                        extraMassTagZone.add(i);
                    }
                }
            }
            modifiedZone.removeAll(extraMassTagZone);
            constraintZone.retainAll(modifiedZone);  //only consider one extraMass now
        }

        if (scanNum == 11557){
            System.out.println("lsz");
        }
        //find all possible ptmCombs and get the bestOne (or best several)
        for (int numPtmsOnPep = 1; numPtmsOnPep <= 4; numPtmsOnPep++){
            collectPtmCombs(scanNum, env, numPtmsOnPep, modifiedZone, idxVarModArrayMap, totalDeltaMass, constraintZone, top1ExtraMass, leftMassBound, rightMassBound, ptmFreePeptide, isDecoy, normalizedCrossCorr, globalRank, peptidePTMPattern, expProcessedPL, plMap, precursorCharge, localMaxMS2Charge);
        }

        return peptidePTMPattern;
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


    private void collectPtmCombs(int scanNum, GRBEnv env, int numPtmsOnPep,  Set<Integer> modifiedZone, Map<Integer, VarModParam[]> idxVarModMap, double totalDeltaMass, Set<Integer> constraintZone, double extraDeltaMass, double leftMassBound, double rightMassBound, String ptmFreePeptide, boolean isDecoy, double normalizedCrossCorr, int globalRank, PeptidePTMPattern peptidePTMPattern, SparseVector expProcessedPL, TreeMap<Double, Double> plMap, int precursorCharge, int localMaxMS2Charge) { // Sometimes, the precursor mass error may affects the digitized spectrum.
        try {
//            GRBEnv env = new GRBEnv(true);
//            env.set(GRB.IntParam.OutputFlag,0);
//            env.start();
            GRBModel model = new GRBModel(env);
//            model.set(GRB.IntParam.OutputFlag, 0);
            double t = 0.02;

            //obj function
            GRBLinExpr objFunction = new GRBLinExpr();
            objFunction.addConstant(t);
            model.setObjective(objFunction, GRB.MINIMIZE);

            //variables
            Map<Integer, GRBVar[]> posVarsMap = new HashMap<>();
            GRBLinExpr totalPtmsOnPepConstr = new GRBLinExpr();
            GRBLinExpr totalMassOnPepConstr = new GRBLinExpr();
            GRBLinExpr totalMassInConstrZoneConstr = new GRBLinExpr();
            for (int pos : modifiedZone){
                if (!idxVarModMap.containsKey(pos)) {
                    System.out.println("lsz Wrong");
                }
                int numPtmsOnAA = idxVarModMap.get(pos).length;
                GRBVar[] varsArray = new GRBVar[numPtmsOnAA];
                for (int i = 0; i < numPtmsOnAA; i++){
                    varsArray[i] = model.addVar(0, 1, 0, GRB.BINARY, "flag_"+pos+"_"+i);
                }
                posVarsMap.put(pos, varsArray);

                GRBLinExpr onePtmOnAaConstr = new GRBLinExpr();
                double[] coeffOneAaArray = new double[numPtmsOnAA];
                Arrays.fill(coeffOneAaArray, 1);
                onePtmOnAaConstr.addTerms(coeffOneAaArray, posVarsMap.get(pos));
                model.addConstr(onePtmOnAaConstr, GRB.LESS_EQUAL, 1, "constr_"+pos);

                totalPtmsOnPepConstr.addTerms(coeffOneAaArray, posVarsMap.get(pos));

                double[] coeffMassAaArray = new double[numPtmsOnAA];
                for (int i = 0; i < numPtmsOnAA; i++){
                    coeffMassAaArray[i] = idxVarModMap.get(pos)[i].mass;
                }
                totalMassOnPepConstr.addTerms(coeffMassAaArray, posVarsMap.get(pos));

                if (constraintZone.contains(pos)) {
                    totalMassInConstrZoneConstr.addTerms(coeffMassAaArray, posVarsMap.get(pos));
                }
            }
            model.addConstr(totalPtmsOnPepConstr, GRB.EQUAL, numPtmsOnPep, "constrTotalNum");
            model.addConstr(totalMassOnPepConstr, GRB.GREATER_EQUAL, totalDeltaMass - t, "constrM1");
            model.addConstr(totalMassOnPepConstr, GRB.LESS_EQUAL, totalDeltaMass + t, "constrM2"); //or put this to constraints as a model.addRange
            if (extraDeltaMass != -9999d) {
                model.addConstr(totalMassInConstrZoneConstr, GRB.GREATER_EQUAL, extraDeltaMass - t, "constrExtraM1");
                model.addConstr(totalMassInConstrZoneConstr, GRB.LESS_EQUAL, extraDeltaMass + t, "constrExtraM2"); //or put this to constraints as a model.addRange
            }
            // solve the model
            model.optimize();
            if (GRB.OPTIMAL != model.get(GRB.IntAttr.Status)){
                model.dispose();
//                env.dispose();
                return;
            }
            int nSolutions = Math.min(model.get(GRB.IntAttr.SolCount), 10);
            if (nSolutions == 0){
                model.dispose();
//                env.dispose();
                return;
            } else {
                for (int solId = 0; solId < nSolutions; solId++){
                    model.set(GRB.IntParam.SolutionNumber, solId);
                    Peptide peptide = new Peptide(ptmFreePeptide, isDecoy, massTool, localMaxMS2Charge, normalizedCrossCorr, globalRank);
                    PositionDeltaMassMap positionDeltaMassMap = new PositionDeltaMassMap(ptmFreePeptide.length());
                    for (int pos : modifiedZone) {
                        int aaHasPtm = 0;
                        double aaMass = 0d;
                        GRBVar[] varsResArray = posVarsMap.get(pos);
                        for (int varId = 0; varId < varsResArray.length; varId++ ) {  //should I be careful about the binary variable to be really 0/1 instead of decimal
                            aaHasPtm += varsResArray[varId].get(GRB.DoubleAttr.Xn);
                            aaMass += varsResArray[varId].get(GRB.DoubleAttr.Xn)*idxVarModMap.get(pos)[varId].mass;
                        }
                        if (aaHasPtm == 1) {
                            positionDeltaMassMap.put(new Coordinate(pos, pos + 1), aaMass);
                        }
                    }
                    peptide.setVarPTM(positionDeltaMassMap); // setVarPTM must be done in one time, otherwise it only get the last PTM
                    double score = massTool.buildVectorAndCalXCorr(peptide.getIonMatrix(), precursorCharge, expProcessedPL);
                    if (score > 0) {
                        peptide.setScore(score);
                        peptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, localMaxMS2Charge, peptide.getIonMatrix(), ms2Tolerance));
                        peptidePTMPattern.update(peptide);
                    }
                }
            }
            model.dispose();
//            env.dispose();
        } catch (GRBException e) {
            System.out.println("Error code: " + e.getErrorCode() + ". " + e.getMessage());
        }
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
                            peptidePTMPattern.update(peptide);
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
                                    peptidePTMPattern.update(peptide);
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
                                                            peptidePTMPattern.update(peptide);
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
                                                                            peptidePTMPattern.update(peptide);
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
                                                                                            peptidePTMPattern.update(peptide);
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
}
