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

import gurobi.*;
import ProteomicsLibrary.*;
import ProteomicsLibrary.MassTool;
import ProteomicsLibrary.Types.*;
//import proteomics.OutputPeff;
import gurobi.GRB;
import gurobi.GRBEnv;
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
    public static final byte N_PART = -2;
    public static final byte C_PART = 2;
    public static final byte N_TERM_PROT = -1;
    public static final byte NON_TERM_PROT = 0;
    public static final byte C_TERM_PROT = 1;
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
            if (modSite == 'M' && modStr[3].contentEquals("Oxidation")) {
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

    private Set<Integer> getInitModZone(String freeSeq, boolean isDecoy, MassTool massTool, int localMaxMS2Charge,
                                        double tagVecScore, int globalRank, SparseVector expProcessedPL, TreeMap<Double, Double> plMap, double cutMass, double ncPart){

        Peptide cleanPep = new Peptide(freeSeq, isDecoy, massTool, localMaxMS2Charge, tagVecScore, globalRank);
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
    private ModPepPool settlePtmOnSide(int scanNum, SparseVector expProcessedPL, TreeMap<Double, Double> plMap, double precursorMass, String partSeq, boolean isDecoy, int tagPosInPep,
                                       Map<Integer, VarPtm[]> posVarPtmArraySrcMap, double cutMass, double deltaMass, int precursorCharge, byte ncPart) throws Exception {
        int localMaxMS2Charge = 1;
        double tagVecScore = -0.99;
        int globalRank = -1;

        TreeMap<Double, Double> unUsedPlMap = new TreeMap<>(plMap);
        Peptide lastPeptide = new Peptide(partSeq, isDecoy, massTool, localMaxMS2Charge, tagVecScore, globalRank);
        lastPeptide.setVarPTM(new PosMassMap(partSeq.length()));
        lastPeptide.shouldPTM = true;

        Set<Integer> modifiedZone = getInitModZone(partSeq, isDecoy, massTool, localMaxMS2Charge, tagVecScore, globalRank, expProcessedPL, plMap, cutMass, ncPart);// e.g. nABCDEFc, modifiedZone={1,2,3,4,5,6}
//        Set<Integer> modifiedZone = IntStream.rangeClosed(0, partSeq.length()-1).boxed().collect(Collectors.toSet());// e.g. nABCDEFc, modifiedZone={1,2,3,4,5,6}
        ModPepPool modPepPoolGood = new ModPepPool(partSeq, 10);
        ModPepPool modPepPoolBad = new ModPepPool(partSeq, 10);

        double massToSettle = deltaMass;
        Set<Integer> toModZone = new HashSet<>(modifiedZone);
        double finalUnsettledMass = 0;
        for (int loop = 1; loop <= 2; loop++) {
            DividedZone dividedZone = dividePepNew(scanNum, toModZone, modPepPoolGood, modPepPoolBad, lastPeptide, posVarPtmArraySrcMap, massToSettle, cutMass, ncPart, partSeq, isDecoy,
                    tagVecScore, globalRank, expProcessedPL, plMap, precursorCharge, 1, unUsedPlMap, 0);

            finalUnsettledMass = massToSettle - dividedZone.ptmMass;
            if (modPepPoolBad.peptideTreeSet.isEmpty() || dividedZone.toModZone.isEmpty() || loop == 2) { //fixme, why the first sentence
                break;
            }
            toModZone = dividedZone.toModZone;
            double massSettled = dividedZone.ptmMass;
            massToSettle -= massSettled;
            lastPeptide = modPepPoolBad.getTopPepPtn();
            modPepPoolBad = new ModPepPool(partSeq, 1);
        }

        // two make up strategy if there is no good pool found
//        if (modPepPoolGood.peptideTreeSet.isEmpty() && ncPart == N_PART && posVarPtmArraySrcMap.containsKey(0)) {
        if (ncPart == N_PART && posVarPtmArraySrcMap.containsKey(0)) { //what if I dont limit only good for dimethyl label dataset
                // Assign pepN or protN varPtm to the nterm aa and try whether it can settle all the mass

            for (VarPtm nVarPtm : posVarPtmArraySrcMap.get(0)) {
                if (nVarPtm.priority != 1) continue;
                ModPepPool tmpModPepPoolBad = new ModPepPool(partSeq, 1);
                double tmpMassToSettle = deltaMass-nVarPtm.mass;
                Set<Integer> tmpToModZone = new HashSet<>(modifiedZone);
                tmpToModZone.remove(0);
                Peptide tmpPeptide = new Peptide(partSeq, isDecoy, massTool, localMaxMS2Charge, tagVecScore, globalRank);
                PosMassMap tmpNewPosMassMap = new PosMassMap(partSeq.length());
                tmpNewPosMassMap.put(0, nVarPtm.mass);
                tmpPeptide.posVarPtmResMap.put(0, nVarPtm);
                tmpPeptide.setVarPTM(tmpNewPosMassMap);
                double[][] ionMatrix = tmpPeptide.getIonMatrixNow();
                updateIonMatrix(ionMatrix, cutMass, ncPart);

                Set<Integer> jRange = IntStream.rangeClosed(0, partSeq.length()-1).boxed().collect(Collectors.toSet());
                double tmpScore = massTool.buildVectorAndCalXCorr(ionMatrix, 1, expProcessedPL, tmpPeptide.matchedBions, tmpPeptide.matchedYions, jRange);
                tmpPeptide.setScore(tmpScore); // init score needs to be bring in the 2-loop

                for (int loop = 1; loop <= 1; loop++) {
                    DividedZone dividedZone = dividePepNew(scanNum, tmpToModZone, modPepPoolGood, tmpModPepPoolBad, tmpPeptide, posVarPtmArraySrcMap, tmpMassToSettle, cutMass, ncPart, partSeq, isDecoy,
                            tagVecScore, globalRank, expProcessedPL, plMap, precursorCharge, 1, unUsedPlMap, 0);
                    if (tmpModPepPoolBad.peptideTreeSet.isEmpty() || dividedZone.toModZone.isEmpty() || loop == 2) { //fixme, why the first sentence
                        break;
                    }
                    tmpToModZone = dividedZone.toModZone;
                    double massSettled = dividedZone.ptmMass;
//                    dividedZone.ptmMass;
                    tmpMassToSettle -= massSettled;
                    tmpPeptide = tmpModPepPoolBad.getTopPepPtn();
                    tmpModPepPoolBad = new ModPepPool(partSeq, 1);
                }
            }
        }
        if (modPepPoolGood.peptideTreeSet.isEmpty() && (massToSettle > -156 && massToSettle < 250)) {
            //  just assign the remaining massToSettle to any aa in toModZone, and record the best one
            double testMass = finalUnsettledMass;
            if (Math.abs(testMass-1*MassTool.PROTON) < 0.02
                    || Math.abs(testMass+1*MassTool.PROTON) < 0.02
                    || Math.abs(testMass-2*MassTool.PROTON) < 0.02 ) {
                modPepPoolGood.push(modPepPoolBad.getTopPepPtn());
                return modPepPoolGood;
            }

            for (int pos : toModZone) { // here must use last toModZone and last massToSettle
                VarPtm fakeVarPtm = new VarPtm(massToSettle, partSeq.charAt(pos), 4, String.format("PIPI_%s", massToSettle), "PIPI_unsettled", -1);
                Peptide fakePeptide = lastPeptide.clone();
//                        Peptide fakePeptide = new Peptide(partSeq, isDecoy, massTool, 1, tagVecScore, globalRank);
                PosMassMap fakeNewPosMassMap = new PosMassMap(partSeq.length());
//                        if (!modPepPoolBad.peptideTreeSet.isEmpty()){
                for (Integer coor : lastPeptide.getVarPTMs().keySet()) {
                    fakeNewPosMassMap.put(coor, lastPeptide.getVarPTMs().get(coor)); // copy the ptms from modPepPoolBad
                }
//                        }
                fakeNewPosMassMap.put(pos, massToSettle);
                fakePeptide.posVarPtmResMap.put(pos, fakeVarPtm);
                fakePeptide.setVarPTM(fakeNewPosMassMap);

                Map<Integer, Double> matchedBions = new HashMap<>();
                Map<Integer, Double> matchedYions = new HashMap<>();
                double[][] ionMatrix = fakePeptide.getIonMatrixNow();
                updateIonMatrix(ionMatrix, cutMass, ncPart);

                int lb = Math.max(0, Collections.min(toModZone));
                int rb = Math.min(partSeq.length() - 1, Collections.max(toModZone));
                Set<Integer> jRange = IntStream.rangeClosed(lb, rb).boxed().collect(Collectors.toSet());
                double fakeScore = massTool.buildVectorAndCalXCorr(ionMatrix, 1, expProcessedPL, matchedBions, matchedYions, jRange);
                fakePeptide.setScore(fakeScore);
                fakePeptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, fakePeptide.getIonMatrix(), ms2Tolerance));
                fakePeptide.matchedBions.putAll(matchedBions);
                fakePeptide.matchedYions.putAll(matchedYions);
                modPepPoolGood.push(fakePeptide);
            }
        }
        return modPepPoolGood;
    }

    public TreeMap<Integer, VarPtm> checkCorrectMass(int scanNum, GRBEnv env,  Map<Integer, VarPtm[]> posVarPtmArraySrcMap, double totalDeltaMass,String partSeq
            , double cutMass, byte isNorC_Part, SparseVector expProcessedPL, double ms1TolAbs) { // Sometimes, the precursor mass error may affects the digitized spectrum.
        int numTimeout = 0;
        int numOfSols = 0;
        Set<Integer> modifiedZone = new HashSet<>(posVarPtmArraySrcMap.keySet());
        modifiedZone.add(-1); // for fakeProton
        ModPepPool modPepPool = new ModPepPool(partSeq, 20);
        int numPtmsOnPep = 2;
        long t1 = 0,t2 = 0,t3 = 0,t4 = 0;
        t1=System.currentTimeMillis();
        try {
            GRBModel model = new GRBModel(env);
            double t = 0.01;
            //obj function
            GRBLinExpr objFunction = new GRBLinExpr();
            int searchMode = 1; //0 min constant obj with 300 sol, 1: max total priority obj with 10 sols

            if (searchMode == 0) {
                objFunction.addConstant(t);
                model.setObjective(objFunction, GRB.MINIMIZE);
            }

            //variables
            Map<Integer, GRBVar[]> posVarsMap = new HashMap<>();
            GRBLinExpr totalPtmsOnPepConstr = new GRBLinExpr();
            GRBLinExpr totalMassOnPepConstr = new GRBLinExpr();
            for (int pos : modifiedZone){
                if (!posVarPtmArraySrcMap.containsKey(pos)) {
                    System.out.println("lsz Wrong");
                }
                if (posVarPtmArraySrcMap.get(pos) == null) {
                    int a = 1;
                }
                int numPtmsOnAA = posVarPtmArraySrcMap.get(pos).length;
                GRBVar[] varsArray = new GRBVar[numPtmsOnAA];
                for (int i = 0; i < numPtmsOnAA; i++){
                    varsArray[i] = model.addVar(0, 1, 0, GRB.BINARY, "isPtmSelected_"+pos+"_"+i);
                }
                posVarsMap.put(pos, varsArray);

                GRBLinExpr onePtmOnAaConstr = new GRBLinExpr();
                double[] coeffOneAaArray = new double[numPtmsOnAA];
                Arrays.fill(coeffOneAaArray, 1);
                onePtmOnAaConstr.addTerms(coeffOneAaArray, posVarsMap.get(pos));
                model.addConstr(onePtmOnAaConstr, GRB.LESS_EQUAL, 1, "max1ptmOn1Pos_"+pos);


                double[] coeffMassAaArray = new double[numPtmsOnAA];
                for (int i = 0; i < numPtmsOnAA; i++){
                    coeffMassAaArray[i] = posVarPtmArraySrcMap.get(pos)[i].mass;
                }

                double[] coeffPriorityAaArray = new double[numPtmsOnAA];
                for (int i = 0; i < numPtmsOnAA; i++){
                    coeffPriorityAaArray[i] = posVarPtmArraySrcMap.get(pos)[i].priority;
                }

                if (pos != -1){  // dont count the fake proton as a ptm , we only need the mass on it
                    totalPtmsOnPepConstr.addTerms(coeffOneAaArray, posVarsMap.get(pos));
                    if (searchMode == 1){
                        objFunction.addTerms(coeffPriorityAaArray, posVarsMap.get(pos));
                    }
                }
                totalMassOnPepConstr.addTerms(coeffMassAaArray, posVarsMap.get(pos));
            }
            if (searchMode == 1) {
                model.setObjective(objFunction, GRB.MAXIMIZE);
            }
            model.addConstr(totalPtmsOnPepConstr, GRB.LESS_EQUAL, 2, "constrTotalNumLT2");
            model.addConstr(totalPtmsOnPepConstr, GRB.GREATER_EQUAL, 1, "constrTotalNumGT1");
            model.addConstr(totalMassOnPepConstr, GRB.GREATER_EQUAL, totalDeltaMass - ms1TolAbs , "constrM1");
            model.addConstr(totalMassOnPepConstr, GRB.LESS_EQUAL, totalDeltaMass +ms1TolAbs, "constrM2"); //or put this to constraints as a model.addRange

            int poolSolutions = 2;
            if (searchMode == 0) {
                poolSolutions = 100;
            }
            model.set(GRB.IntParam.MIPFocus, 1);
            model.set(GRB.IntParam.PoolSearchMode, 1);
            model.set(GRB.IntParam.PoolSolutions, poolSolutions);
            model.set(GRB.DoubleParam.TimeLimit, 1);

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
                    return new TreeMap<>();
            }
            int solCount = model.get(GRB.IntAttr.SolCount);
            if (solCount == 0) return new TreeMap<>();
            for (int solNum = 0; solNum < solCount; solNum++) {
                model.set(GRB.IntParam.SolutionNumber, solNum);
                Map<Integer, Integer> positionOneMap = new HashMap<>();
                Peptide peptide = new Peptide(partSeq, false, massTool, 2, -0.99, 1); // todo abandoned para should be deleted
                PosMassMap posMassMap = new PosMassMap(partSeq.length());
                for (int pos : modifiedZone) {
                    if (pos == -1 ) continue; // skip the fake proton
                    GRBVar[] varsResArray = posVarsMap.get(pos);
                    for (int varId = 0; varId < varsResArray.length; varId++ ) {  //should I be careful about the binary variable to be really 0/1 instead of decimal
                        if (1 == Math.round(varsResArray[varId].get(GRB.DoubleAttr.Xn))) {
                            posMassMap.put(pos, posVarPtmArraySrcMap.get(pos)[varId].mass);
                            peptide.totalPriority += posVarPtmArraySrcMap.get(pos)[varId].priority;
                            peptide.posVarPtmResMap.put(pos, posVarPtmArraySrcMap.get(pos)[varId]);
                            break; // one pos only has at most 1 ptm
                        }
                    }
                }
                peptide.setVarPTM(posMassMap);
                if (posMassMap.containsKey(4) && posMassMap.containsKey(7)){
                    int a = 1;
                }
                double[][] ionMatrix = peptide.getIonMatrixNow();
                updateIonMatrix(ionMatrix, cutMass, isNorC_Part);

                Set<Integer> jRange = IntStream.rangeClosed(0, partSeq.length()-1).boxed().collect(Collectors.toSet());
                double score = massTool.buildVectorAndCalXCorr(ionMatrix, 1, expProcessedPL, peptide.matchedBions, peptide.matchedYions, jRange) ;

                peptide.setScore(score);
//                peptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, peptide.getIonMatrix(), ms2Tolerance));
                modPepPool.push(peptide);
            }

            t4=System.currentTimeMillis();
            System.out.println(scanNum + ", buildTime_Len,solve_Num,verify:" + (t2-t1)+"_"+modifiedZone.size() +","+ (t3-t2)+"_"+solCount+","+ (t4-t3));
            model.dispose();
        } catch (GRBException e) {
            System.out.println("Error code: " + e.getErrorCode() + ". " + e.getMessage());
        }
        if (modPepPool.peptideTreeSet.isEmpty()) {
            return new TreeMap<>();
        }
        return modPepPool.getTopPepPtn().posVarPtmResMap;
    }
    public ModPepPool findPtmNew2(int scanNum, SparseVector expProcessedPL, TreeMap<Double, Double> plMap, double precursorMass, Peptide candiPep, PeptideInfo peptideInfo
            , int precursorCharge, int localMaxMS2Charge, double localMS1ToleranceL, double localMS1ToleranceR) throws Exception  {
//        double ptmFreeMass = massTool.calResidueMass(ptmFreePeptide) + massTool.H2O;
        String freeSeq = candiPep.getFreeSeq();
        boolean isDecoy = candiPep.isDecoy;
        double tagVecScore = candiPep.getTagVecScore();
        int globalRank = candiPep.getGlobalRank();
        char leftFlank = peptideInfo.leftFlank;
        char rightFlank = peptideInfo.rightFlank;

        double nDeltaMass = candiPep.nDeltaMass;
        double cDeltaMass = candiPep.cDeltaMass;

//        if (Math.abs(nDeltaMass) <= 0.2 && Math.abs(cDeltaMass) <= 0.2) { //this must be a candidate borrowed from other scans and with no delta mass both sides but has total delta mass
//            return findPtmNew1(scanNum, expProcessedPL,  plMap, precursorMass, candiPep, peptideInfo, precursorCharge, localMaxMS2Charge, localMS1ToleranceL, localMS1ToleranceR);
//        }

        ModPepPool fullModPepPool = new ModPepPool(freeSeq,5);

        int tagPosInPep = candiPep.tagPosInPep;

        int tagLen = 0;

        double score = 0;
        if (candiPep.finderTag != null) {
            score = candiPep.finderTag.getTotalIntensity();
            tagLen = candiPep.finderTag.size();
        }

//        ModPeptides modPeptides = new ModPeptides(freeSeq);
//        Peptide peptide = new Peptide(freeSeq, isDecoy, massTool, localMaxMS2Charge, tagVecScore, globalRank);
//
//        double totalDeltaMass = precursorMass - peptide.getTheoMass();
//        peptide.absDeltaMass = totalDeltaMass;

        Set<Integer> fixModIdxes = getFixModIdxes(freeSeq, fixModMap);  // positions that has fixed mod on it. Those postions should not bear var mod then. Unless it is at N term

        int pepPosInProt = NON_TERM_PROT; // none of the terms is protein term
        if (leftFlank == '-') {
            pepPosInProt = N_TERM_PROT; //nterm is protein N term
        } else if (rightFlank == '-') {
            pepPosInProt = C_TERM_PROT;  //cterm is protein C term
        }

//        Map<Integer, Set<VarPtm>> idxVarModMap = getIdxVarModMapNew(freeSeq, fixModIdxes,  tagPosInPep, tagLen); //todo no need to generate var mod list for aa again and again, make it stored.
        Map<Integer, VarPtm[]> posVarPtmArrayMap = getIdxVarModMapNew(freeSeq, fixModIdxes,  tagPosInPep, tagLen);
//        for (int id : idxVarModMap.keySet()){
//            VarPtm[] modArray = new VarPtm[idxVarModMap.get(id).size()];
//            idxVarModMap.get(id).toArray(modArray);
//            Arrays.sort(modArray, Comparator.comparingDouble(VarPtm::getMass));
//            posVarPtmArrayMap.put(id, modArray);
//        }


        PosMassMap fullPosMassMap = new PosMassMap(freeSeq.length());
        TreeMap<Integer, VarPtm> posVarPtmResMap = new TreeMap<>();

        if (Math.abs(nDeltaMass) > 0.2) { // nterm has unsettled mass and must has untagged amino acid for var mod
            //use a partial peptide and cterm mass to settle ptm on n side
            String nPartSeq = freeSeq.substring(0,tagPosInPep);
            Map<Integer, VarPtm[]> partPosVarModArrayMap = new HashMap<>();
            for (int i = 0; i < tagPosInPep; i++) {
                if (posVarPtmArrayMap.containsKey(i)){
                    partPosVarModArrayMap.put(i, posVarPtmArrayMap.get(i));
                }
            }
            double cutMass = precursorMass + MassTool.PROTON - candiPep.finderTag.getHeadLocation() - massTool.H2O;
            ModPepPool nPartModPepsSettled = settlePtmOnSide(scanNum, expProcessedPL, plMap, precursorMass, nPartSeq, isDecoy, tagPosInPep,
                    partPosVarModArrayMap, cutMass, nDeltaMass, precursorCharge, N_PART);
            //todo  record the nPartModPepsSettled to the all part for fixed n part
            if (nPartModPepsSettled.peptideTreeSet.isEmpty()) {
                int a = 1;
//                System.out.println(scanNum+", n peptideTreeSet.isEmpty(),"+candiPep.getFreeSeq()+","+String.join("_", peptideInfo.protIdSet));
                return fullModPepPool; //empty
            }

            List<Peptide> nPartPeptideList = new ArrayList<>(nPartModPepsSettled.peptideTreeSet);
            Collections.sort(nPartPeptideList, Comparator.comparing(o -> o.getPriority(), Comparator.reverseOrder()));
            for (Integer coor : nPartPeptideList.get(0).getVarPTMs().keySet()) { //copy the top 1 ptm pattern in n part // whhat if also choose the largeset priority
                fullPosMassMap.put(coor, nPartPeptideList.get(0).getVarPTMs().get(coor)); // copy the ptms from partModPepsUnsettled
                posVarPtmResMap.put(coor, nPartPeptideList.get(0).posVarPtmResMap.get(coor));
            }
//            posVarPtmResMap.putAll(nPartModPepsSettled.getTopPepPtn().posVarPtmResMap);
        }

        if (Math.abs(cDeltaMass) > 0.2) { // nterm has unsettled mass and must has untagged amino acid for var mod
//            TreeMap<Double, Double> unUsedPlMap = new TreeMap<>(plMap);
            String cPartSeq = freeSeq.substring(tagPosInPep+tagLen);
            Map<Integer, VarPtm[]> partPosVarModArrayMap = new HashMap<>();
            for (int i = tagPosInPep+tagLen; i < freeSeq.length(); i++) {
                if (posVarPtmArrayMap.containsKey(i)){
                    partPosVarModArrayMap.put(i-tagPosInPep-tagLen, posVarPtmArrayMap.get(i));
                }
            }
            double cutMass = candiPep.finderTag.getTailLocation() - MassTool.PROTON;
            ModPepPool cPartModPepsSettled = settlePtmOnSide(scanNum, expProcessedPL, plMap, precursorMass, cPartSeq, isDecoy, tagPosInPep,
                    partPosVarModArrayMap, cutMass, cDeltaMass, precursorCharge, C_PART);
            //todo  record the nPartModPepsSettled to the all part for fixed n part
            if (cPartModPepsSettled.peptideTreeSet.isEmpty()) {
                int a = 1;
//                System.out.println(scanNum+", cPeptideTreeSet.isEmpty(),,"+candiPep.getFreeSeq()+","+String.join("_", peptideInfo.protIdSet));
                return fullModPepPool; //empty
            }
            List<Peptide> cPartPeptideList = new ArrayList<>(cPartModPepsSettled.peptideTreeSet);
            Collections.sort(cPartPeptideList, Comparator.comparing(o -> o.getPriority(), Comparator.reverseOrder()));
            for (Integer coor : cPartPeptideList.get(0).getVarPTMs().keySet()) { //copy the top 1 ptm pattern in n part // whhat if also choose the largeset priority
                fullPosMassMap.put(coor+tagPosInPep+tagLen, cPartPeptideList.get(0).getVarPTMs().get(coor)); // copy the ptms from partModPepsUnsettled
                posVarPtmResMap.put(coor+tagPosInPep+tagLen, cPartPeptideList.get(0).posVarPtmResMap.get(coor));
            }
        }


        int idOfAa = -1;
        for (char letter : candiPep.finderTag.getPtmAaString().toCharArray()) {
            if (Character.isUpperCase(letter)) {
                idOfAa += 1;
            } else {
                fullPosMassMap.put(candiPep.tagPosInPep+idOfAa, massTool.labelVarPtmMap.get(letter).mass);
                posVarPtmResMap.put(candiPep.tagPosInPep+idOfAa, massTool.labelVarPtmMap.get(letter));
            }
        }

        Peptide fullPeptide = candiPep.clone();
        fullPeptide.setVarPTM(fullPosMassMap);
        if (fullPeptide.getVarPTMs() != null && fullPeptide.getVarPTMs().isEmpty()) {
            System.out.println("wrong "+scanNum +","+fullPeptide.getFreeSeq()+","+String.join("_", peptideInfo.protIdSet));
        }
        fullPeptide.setScore(massTool.buildVectorAndCalXCorr(fullPeptide.getIonMatrixNow(), 1, expProcessedPL, fullPeptide.matchedBions,  fullPeptide.matchedYions) - 0*0.1);//todo decide the penalty
        fullPeptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, fullPeptide.getIonMatrixNow(), ms2Tolerance));
        fullPeptide.shouldPTM = true;
        fullPeptide.posVarPtmResMap = posVarPtmResMap;
        fullModPepPool.push(fullPeptide);

        return fullModPepPool;
    }

    public TreeMap<Integer, VarPtm> settleSideWithPtm(int scanNum, String partSeq, double deltaMass, byte isNorC_Side, byte isProtNorC_Term, double cutMass, SparseVector expProcessedPL,double ms1TolAbs) throws GRBException {
        Set<Integer> fixModIdxes = getFixModIdxes(partSeq, fixModMap);  // positions that has fixed mod on it. Those postions should not bear var mod then. Unless it is at N term
        Map<Integer, VarPtm[]> partPosVarModArrayMap = getIdxVarModMapNew(partSeq, fixModIdxes, isNorC_Side, isProtNorC_Term); //todo no need to generate var mod list for aa again and again, make it stored.
        GRBEnv env = new GRBEnv(true);
        env.set(GRB.IntParam.OutputFlag,0);
        env.start();
//        Peptide peptide = checkCorrectMass(scanNum, env, partPosVarModArrayMap, deltaMass, partSeq , cutMass, isNorC_Side, expProcessedPL, ms1TolAbs);

        return checkCorrectMass(scanNum, env, partPosVarModArrayMap, deltaMass, partSeq , cutMass, isNorC_Side, expProcessedPL, ms1TolAbs);
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
    private DividedZone dividePep(int scanNum, Set<Integer> modifiedZone, ModPepPool ptmGoodRes, ModPepPool ptmBadRes, ModPepPool ptmTemplate, Map<Integer, VarPtm[]> idxVarModMap, double totalDeltaMass, String ptmFreePeptide, boolean isDecoy, double normalizedCrossCorr, int globalRank, SparseVector expProcessedPL, TreeMap<Double, Double> plMap, int precursorCharge, int localMaxMS2Charge, SortedMap<Double, Double> unUsedPlMap) { // Sometimes, the precursor mass error may affects the digitized spectrum.
        double ptmMass = 0;
        for (int pos : modifiedZone) {
            if (!idxVarModMap.containsKey(pos)) continue;

            for (int ptmId = 0; ptmId < idxVarModMap.get(pos).length; ptmId++){
                Peptide peptide = new Peptide(ptmFreePeptide, isDecoy, massTool, localMaxMS2Charge, normalizedCrossCorr, globalRank);
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

    public Map<Integer, VarPtm[]> getIdxVarModMapNew(String partSeq, Set<Integer> fixModIdxes, int isNorC_Part, int isProtNorC_Term) {
//        partSeq, fixModIdxes, isNorC_Side, isProtNorC_Term
        Map<Integer, VarPtm[]> idxVarModMap = new HashMap<>(partSeq.length() + 1, 1);
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
                Map<String, VarPtm> dstMap = new HashMap<>();
                List<VarPtm> srcSet = finalPtmMap.get(aa);
                if (i == 0) { //aa at seq n term
                    for (VarPtm varPtm : srcSet) {
                        if (massTable.get(partSeq.charAt(i)) + varPtm.mass < ms2Tolerance) continue;//the mass of a modified amino acid cannot be 0 or negative

                        if (varPtm.position == 4 || (varPtm.position == 2 && isNorC_Part == N_PART) || (varPtm.position == 0 && hasProt_N_TermPtm)) { // anywhere or pepN or (protN and pepPos at protN)
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

                        if (varPtm.position == 4 || (varPtm.position == 3 && isNorC_Part == C_PART) || (varPtm.position == 1 && hasProt_C_TermPtm)) { // anywhere or pepC or (protC and pepPos at protC)
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

        // add fake amino acid and possible +- proton to solve the precursor error problem
        VarPtm[] modArray = new VarPtm[2];
        modArray[0] = new VarPtm(-MassTool.PROTON, 'Z', 4, "Fake_Proton_-1", "By_author", 0);
        modArray[1] = new VarPtm(MassTool.PROTON, 'Z', 4, "Fake_Proton_1", "By_author", 0);
//        modArray[0] = new VarPtm(2*MassTool.PROTON, 'Z', 4, "Fake_Proton", "By_author", 0);
        idxVarModMap.put(-1, modArray);
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
