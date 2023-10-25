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
import com.google.common.collect.Sets;
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

import static proteomics.PIPI.*;


public class InferPTM {
    private static final int PoolSolutions = 100;
    private static final int MaxPtmNumInPart = 3;
    private static final Pattern pattern = Pattern.compile("([0-9A-Za-z]+)(\\(([0-9\\-]+)\\))?");
    private static final double ptmMassTolerance = 0.1;
    public static final byte N_PART = 0;
    public static final byte C_PART = 1;
    public static final byte N_TERM_PROT = -1;
    public static final byte NON_TERM_PROT = 0;
    public static final byte C_TERM_PROT = 1;
    public static final byte ANYWHERE = 4;
    public static final byte PEPC = 3;
    public static final byte PEPN = 2;
    public static final byte PROTC = 1;
    public static final byte PROTN = 0;
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
    private Map<Character, List<VarPtm>> finalPtmMap = new HashMap<>(22);
    private Map<Character, Map<Byte, List<VarPtm>>> aaAllVarPtmMap = new HashMap<>(22);
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
            byte modPosition = Byte.valueOf(modStr[2]);
            int priority = 1;
            if (modSite == 'M' && modStr[3].contentEquals("Oxidation") && n_varPtm != 0) {
//                if (modSite == 'M' && modStr[3].contentEquals("Oxidation") && n_varPtm != 1) {
                priority = 0; // oxidation on M is a common phonomenon but not an enriched modification
            }
            if (modPosition == ANYWHERE) {//  position anywhere, highest prority
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

        aaAllVarPtmMap = new HashMap<>(22); //<aachar, <position, varptm>>
        for (Character aa : finalPtmMap.keySet()){
            Map<Byte, List<VarPtm>> positionVarPtmMap = new HashMap<>(5);
            for (VarPtm varPtm : finalPtmMap.get(aa)) {
                byte position = varPtm.position;
                if (positionVarPtmMap.containsKey(position)) {
                    positionVarPtmMap.get(position).add(varPtm);
                } else {
                    List<VarPtm> varPtmList = new ArrayList<>(50);
                    varPtmList.add(varPtm);
                    positionVarPtmMap.put(position, varPtmList);
                }
            }
            aaAllVarPtmMap.put(aa, positionVarPtmMap);
        }
        int a = 1;
    }

    public List<Map<Integer, Integer>> findBestPtmInN_Spec(int scanNum, GRBEnv env,  String partSeq, int refPos, double ms1TolAbs, double flexiableMass,
                                                           Map<Integer, Integer> posYIdMap,Map<Integer, Integer> yIdMinPosMap, byte isProtNorC_Term) {
        Map<Integer, List<VarPtm>> posAllVarPtmMapaa = new HashMap<>();
        List<Map<Integer, Integer>> posPtmIdResList = new LinkedList<>();
        try {
            GRBModel model = new GRBModel(env);
            double t = 0.01;

            //variables
            GRBLinExpr totalNumsOnPep = new GRBLinExpr();
            GRBLinExpr totalFlexiableMass = new GRBLinExpr();

            Map<Integer, GRBVar> yVarMap = new HashMap<>(); // <yId, var>
            Map<Integer, Map<Integer, GRBVar>> xVarMap = new HashMap<>(); // <pos, <ptmId, var>>

            Map<Integer, List<Integer>> xPosFlexiblePtmIdMap = new HashMap<>();
            for (int pos = 0; pos < partSeq.length(); pos++) {
                int absPos = pos + refPos;
                char aa = partSeq.charAt(pos);
                boolean hasFixMod = hasFixMod(aa);
                GRBLinExpr sumX_leq_1orY = new GRBLinExpr();
                if ( ! hasFixMod) {
                    posAllVarPtmMapaa.put(absPos, new ArrayList<>());
                    Map<Integer, GRBVar> ptmIdXVarMap = new HashMap<>();
                    int ptmId = 0;
                    for (VarPtm varPtm : aaAllVarPtmMap.get(aa).get(ANYWHERE)) {
                        GRBVar xVar = model.addVar(0, 1, 0, GRB.BINARY, "x"+absPos+"_"+varPtm.mass);
                        ptmIdXVarMap.put(ptmId, xVar);
                        posAllVarPtmMapaa.get(absPos).add(varPtm);
                        sumX_leq_1orY.addTerm(1, xVar);
                        totalFlexiableMass.addTerm(varPtm.mass, xVar); // + m_i * x_i
                        totalNumsOnPep.addTerm(1,xVar); // + 1 * x_i
                        ptmId++;
                    }
                    if (yIdMinPosMap.containsValue(absPos) || pos == 0) { // is the pos is the max pos (end pos) of each yId poses, they could be c term
                        xPosFlexiblePtmIdMap.put(absPos, new ArrayList<>()); // not including the last pos because no need to.
                        for (VarPtm varPtm : aaAllVarPtmMap.get(aa).get(PEPN)) {
                            GRBVar xVar = model.addVar(0, 1, 0, GRB.BINARY, "x"+absPos+"_"+varPtm.mass);
                            ptmIdXVarMap.put(ptmId, xVar);
                            posAllVarPtmMapaa.get(absPos).add(varPtm);

                            sumX_leq_1orY.addTerm(1, xVar);
                            totalFlexiableMass.addTerm(varPtm.mass, xVar); // + m_i * x_i
                            totalNumsOnPep.addTerm(1,xVar); // + 1 * x_i
                            xPosFlexiblePtmIdMap.get(absPos).add(ptmId);

                            ptmId++;
                        }
                        if (isProtNorC_Term == N_TERM_PROT && absPos <= 1) {
                            for (VarPtm varPtm : aaAllVarPtmMap.get(aa).get(PROTN)) {
                                GRBVar xVar = model.addVar(0, 1, 0, GRB.BINARY, "x"+absPos+"_"+varPtm.mass);
                                ptmIdXVarMap.put(ptmId, xVar);
                                posAllVarPtmMapaa.get(absPos).add(varPtm);

                                sumX_leq_1orY.addTerm(1, xVar);
                                totalFlexiableMass.addTerm(varPtm.mass, xVar); // + m_i * x_i
                                totalNumsOnPep.addTerm(1,xVar); // + 1 * x_i
                                ptmId++;
                            }
                        }
                    }
                    xVarMap.put(absPos, ptmIdXVarMap);
                }

                if (posYIdMap.containsKey(absPos)){
                    int yId = posYIdMap.get(absPos);
                    if (yVarMap.containsKey(yId)) {
                    } else {
                        GRBVar yVar = model.addVar(0, 1, 0, GRB.BINARY, "y"+yId);
                        yVarMap.put(yId, yVar);
                    }
                    if ( ! hasFixMod) {
                        sumX_leq_1orY.addTerm(-1, yVarMap.get(yId));
                    }
//                    sumX_leq_1orY.addTerm(-1, yVar);
                    double aaMass = massTool.getMassTable().get(partSeq.charAt(pos));
                    totalFlexiableMass.addTerm(aaMass, yVarMap.get(yId));

                    if ( ! hasFixMod) {
                        model.addConstr(sumX_leq_1orY, GRB.LESS_EQUAL, 0, "sumX_leq_Y" + absPos);
                    }
                } else {
                    if ( ! hasFixMod) {
                        model.addConstr(sumX_leq_1orY, GRB.LESS_EQUAL, 1, "sumX_leq_1" + absPos);
                    }
                }
            }

            if (yVarMap.size() >= 2) {
                List<Integer> yIdList = new ArrayList<>(yVarMap.keySet());
                yIdList.sort(Comparator.naturalOrder());
                for (int i = 0; i < yIdList.size()-1; i++) {
                    int thisYId = yIdList.get(i);
                    int nextYId = yIdList.get(i+1);
                    GRBLinExpr yiyj = new GRBLinExpr();
                    yiyj.addTerm(1, yVarMap.get(thisYId));
                    yiyj.addTerm(-1, yVarMap.get(nextYId));
                    model.addConstr(yiyj, GRB.GREATER_EQUAL, 0 , "yiyj_"+i);

                    int thisMinPos = yIdMinPosMap.get(thisYId);
                    if (! hasFixMod(partSeq.charAt(thisMinPos-refPos))) {
                        for (int ptmId : xPosFlexiblePtmIdMap.get(thisMinPos)) {
                            GRBLinExpr yiyjxi = new GRBLinExpr();
                            yiyjxi.addTerm(1, yVarMap.get(thisYId));
                            yiyjxi.addTerm(-1, yVarMap.get(nextYId));
                            yiyjxi.addTerm(-1, xVarMap.get(thisMinPos).get(ptmId));
                            model.addConstr(yiyjxi, GRB.GREATER_EQUAL, 0 , "yiyjxi"+i);
                        }
                    }
                }
            }

            model.addConstr(totalFlexiableMass, GRB.GREATER_EQUAL, flexiableMass - ms1TolAbs , "constrM1");
            model.addConstr(totalFlexiableMass, GRB.LESS_EQUAL, flexiableMass + ms1TolAbs, "constrM2"); //or put this to constraints as a model.addRange

            model.addConstr(totalNumsOnPep, GRB.LESS_EQUAL, 3, "totalNumsOnPep"); // this value should not exceed the number of aa in the partSeq

            //obj function
            model.setObjective(totalNumsOnPep, GRB.MINIMIZE); // with this, the solver will find those with smaller number of PTM first then more numbers

//            if (lszDebugScanNum.contains(scanNum) && partSeq.contentEquals("KFGVLSDNFK")){
//                int a = 1;
//            }
            int poolSolutions = 1000;
            model.set(GRB.IntParam.MIPFocus, 1);
            model.set(GRB.IntParam.PoolSearchMode, 2);
            model.set(GRB.IntParam.PoolSolutions, poolSolutions);
            model.set(GRB.DoubleParam.TimeLimit, 1); // second
            model.set(GRB.DoubleParam.PoolGap, 1); // if the best obj is 1, then all sol with 1*(1+poolGap) with be discarded and save time
            model.optimize();
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
                Map<Integer, Integer> thisSol = new HashMap<>();
                for (int absPos = refPos; absPos < partSeq.length() + refPos; absPos++) {
                    if (posYIdMap.containsKey(absPos)) { //just for verification
                        int yId = posYIdMap.get(absPos);
                        GRBVar yVar = yVarMap.get(yId);
                        if ((int) Math.round(yVar.get(GRB.DoubleAttr.Xn)) == 1) {
                            thisSol.put(absPos, -1); //-1 means no ptm yet
                            if (xVarMap.containsKey(absPos)) {
                                Map<Integer, GRBVar> ptmIdXVarMap = xVarMap.get(absPos);
                                for (int ptmId : ptmIdXVarMap.keySet()) {
                                    GRBVar xVar = ptmIdXVarMap.get(ptmId);
                                    if ((int) Math.round(xVar.get(GRB.DoubleAttr.Xn)) == 1) {
                                        thisSol.put(absPos, ptmId);
                                        break;
                                    }
                                }
                            }
                        }
                    } else {
                        thisSol.put(absPos, -1); //-1 means no ptm yet
                        if (xVarMap.containsKey(absPos)){  // when there is a C, the pos wont be in xVarMap but we still record it in thisSol
                            Map<Integer, GRBVar> ptmIdXVarMap = xVarMap.get(absPos);
                            for (int ptmId : ptmIdXVarMap.keySet()) {
                                GRBVar xVar = ptmIdXVarMap.get(ptmId);
                                if ((int) Math.round(xVar.get(GRB.DoubleAttr.Xn)) == 1) {
                                    thisSol.put(absPos, ptmId);
                                    break;
                                }
                            }
                        }
                    }
                }
                posPtmIdResList.add(thisSol);
            }
            model.dispose();
        } catch (GRBException e) {
            System.out.println(scanNum + ","  + " Error code: " + e.getErrorCode() + ". " + e.getMessage());
        }
        return posPtmIdResList;
    }
    public List<Map<Integer, Integer>> findBestPtmInC_Spec(int scanNum, GRBEnv env, String partSeq, int refPos,
                                                           double ms1TolAbs, double flexiableMass, Map<Integer, Integer> posYIdMap,
                                                           Map<Integer, Integer> yIdMaxPosMap, byte isProtNorC_Term) {
        Map<Integer, List<VarPtm>> posAllVarPtmMapaa = new HashMap<>();
        List<Map<Integer, Integer>> posPtmIdResList = new LinkedList<>();
        try {
            GRBModel model = new GRBModel(env);
            double t = 0.01;

            //variables
            GRBLinExpr totalNumsOnPep = new GRBLinExpr();
            GRBLinExpr totalFlexiableMass = new GRBLinExpr();

            Map<Integer, GRBVar> yVarMap = new HashMap<>(); // <yId, var>
            Map<Integer, Map<Integer, GRBVar>> xVarMap = new HashMap<>(); // <pos, <ptmId, var>>
            Map<Integer, List<Integer>> xPosFlexiblePtmIdMap = new HashMap<>();
            for (int pos = 0; pos < partSeq.length(); pos++) {
                char aa = partSeq.charAt(pos);
                int absPos = refPos + pos;
                boolean hasFixMod = hasFixMod(aa);
                GRBLinExpr sumX_leq_1orY = new GRBLinExpr();
                if ( ! hasFixMod) {
                    posAllVarPtmMapaa.put(absPos, new ArrayList<>());
                    Map<Integer, GRBVar> ptmIdXVarMap = new HashMap<>();
                    int ptmId = 0;
                    for (VarPtm varPtm : aaAllVarPtmMap.get(aa).get(ANYWHERE)) {
                        GRBVar xVar = model.addVar(0, 1, 0, GRB.BINARY, "x"+absPos+"_"+varPtm.mass);
                        ptmIdXVarMap.put(ptmId, xVar);
                        posAllVarPtmMapaa.get(absPos).add(varPtm);
                        sumX_leq_1orY.addTerm(1, xVar);
                        totalFlexiableMass.addTerm(varPtm.mass, xVar); // + m_i * x_i
                        totalNumsOnPep.addTerm(1,xVar); // + 1 * x_i
                        ptmId++;
                    }
                    if (yIdMaxPosMap.containsValue(absPos) || pos == partSeq.length()-1) { // is the pos is the max pos (end pos) of each yId poses, they could be c term
                        xPosFlexiblePtmIdMap.put(absPos, new ArrayList<>()); // not including the last pos because no need to.
                        for (VarPtm varPtm : aaAllVarPtmMap.get(aa).get(PEPC)) {
                            GRBVar xVar = model.addVar(0, 1, 0, GRB.BINARY, "x"+absPos+"_"+varPtm.mass);
                            ptmIdXVarMap.put(ptmId, xVar);
                            posAllVarPtmMapaa.get(absPos).add(varPtm);

                            sumX_leq_1orY.addTerm(1, xVar);
                            totalFlexiableMass.addTerm(varPtm.mass, xVar); // + m_i * x_i
                            totalNumsOnPep.addTerm(1,xVar); // + 1 * x_i
                            xPosFlexiblePtmIdMap.get(absPos).add(ptmId);

                            ptmId++;
                        }
                        if (isProtNorC_Term == C_TERM_PROT && pos == partSeq.length()-1) {
                            for (VarPtm varPtm : aaAllVarPtmMap.get(aa).get(PROTC)) {
                                GRBVar xVar = model.addVar(0, 1, 0, GRB.BINARY, "x"+absPos+"_"+varPtm.mass);
                                ptmIdXVarMap.put(ptmId, xVar);
                                posAllVarPtmMapaa.get(absPos).add(varPtm);

                                sumX_leq_1orY.addTerm(1, xVar);
                                totalFlexiableMass.addTerm(varPtm.mass, xVar); // + m_i * x_i
                                totalNumsOnPep.addTerm(1,xVar); // + 1 * x_i
                                ptmId++;
                            }
                        }
                    }
                    xVarMap.put(absPos, ptmIdXVarMap);
                }

                if (posYIdMap.containsKey(absPos)){
                    int yId = posYIdMap.get(absPos);
                    if (yVarMap.containsKey(yId)) {
                    } else {
                        GRBVar yVar = model.addVar(0, 1, 0, GRB.BINARY, "y"+yId);
                        yVarMap.put(yId, yVar);
                    }
                    if ( ! hasFixMod) {
                        sumX_leq_1orY.addTerm(-1, yVarMap.get(yId));
                    }
//                    sumX_leq_1orY.addTerm(-1, yVar);
                    double aaMass = massTool.getMassTable().get(partSeq.charAt(pos));
                    totalFlexiableMass.addTerm(aaMass, yVarMap.get(yId));

                    if ( ! hasFixMod) {
                        model.addConstr(sumX_leq_1orY, GRB.LESS_EQUAL, 0, "sumX_leq_Y" + absPos);
                    }
                } else {
                    if ( ! hasFixMod) {
                        model.addConstr(sumX_leq_1orY, GRB.LESS_EQUAL, 1, "sumX_leq_1" + absPos);
                    }
                }
            }
            if (yVarMap.size() >= 2) {
                List<Integer> yIdList = new ArrayList<>(yVarMap.keySet());
                yIdList.sort(Comparator.naturalOrder());
                for (int i = 0; i < yIdList.size()-1; i++) {
                    int thisYId = yIdList.get(i);
                    int nextYId = yIdList.get(i+1);
                    GRBLinExpr yiyj = new GRBLinExpr();
                    yiyj.addTerm(1, yVarMap.get(thisYId));
                    yiyj.addTerm(-1, yVarMap.get(nextYId));
                    model.addConstr(yiyj, GRB.GREATER_EQUAL, 0 , "yiyj_"+i);

                    int thisMaxPos = yIdMaxPosMap.get(thisYId);
                    if (! hasFixMod(partSeq.charAt(thisMaxPos - refPos))) {
                        for (int ptmId : xPosFlexiblePtmIdMap.get(thisMaxPos)) {
                            GRBLinExpr yiyjxi = new GRBLinExpr();
                            yiyjxi.addTerm(1, yVarMap.get(thisYId));
                            yiyjxi.addTerm(-1, yVarMap.get(nextYId));
                            yiyjxi.addTerm(-1, xVarMap.get(thisMaxPos).get(ptmId));
                            model.addConstr(yiyjxi, GRB.GREATER_EQUAL, 0, "yiyjxi" + i);
                        }
                    }
                }
            }
            model.addConstr(totalFlexiableMass, GRB.GREATER_EQUAL, flexiableMass - ms1TolAbs , "constrM1");
            model.addConstr(totalFlexiableMass, GRB.LESS_EQUAL, flexiableMass + ms1TolAbs, "constrM2"); //or put this to constraints as a model.addRange

            model.addConstr(totalNumsOnPep, GRB.LESS_EQUAL, 3.1, "totalNumsOnPep"); // this value should not exceed the number of aa in the partSeq

            //obj function
            model.setObjective(totalNumsOnPep, GRB.MINIMIZE); // with this, the solver will find those with smaller number of PTM first then more numbers

//            if (lszDebugScanNum.contains(scanNum) && partSeq.contentEquals("KFGVLSDNFK")){
//                int a = 1;
//            }
            int poolSolutions = 1000;
            model.set(GRB.IntParam.MIPFocus, 1); // 2 seems better than 1 but dont know why
            model.set(GRB.IntParam.PoolSearchMode, 2 );  //0 for only one sol, 1 for possible more but not guaranteed = poolSolutions, 2 for guaranteed but long time
            model.set(GRB.IntParam.PoolSolutions, poolSolutions);
            model.set(GRB.DoubleParam.TimeLimit, 1); // second
//            model.set(GRB.DoubleParam.Heuristics, 0.05); //
//            model.set(GRB.DoubleParam.PoolGap, 1); // if the best obj is 1, then all sol with 1*(1+poolGap) with be discarded and save time

            model.optimize();
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
                Map<Integer, Integer> thisSol = new HashMap<>();
                for (int absPos = refPos; absPos < partSeq.length()+refPos; absPos++) {
                    if (posYIdMap.containsKey(absPos)) { //just for verification
                        int yId = posYIdMap.get(absPos);
                        GRBVar yVar = yVarMap.get(yId);
                        if ((int) Math.round(yVar.get(GRB.DoubleAttr.Xn)) == 1) {
                            thisSol.put(absPos, -1); //-1 means no ptm yet
                            if (xVarMap.containsKey(absPos)) {
                                Map<Integer, GRBVar> ptmIdXVarMap = xVarMap.get(absPos);
                                for (int ptmId : ptmIdXVarMap.keySet()) {
                                    GRBVar xVar = ptmIdXVarMap.get(ptmId);
                                    if ((int) Math.round(xVar.get(GRB.DoubleAttr.Xn)) == 1) {
                                        thisSol.put(absPos, ptmId);
                                        break;
                                    }
                                }
                            }
                        }
                    } else {
                        thisSol.put(absPos, -1); //-1 means no ptm yet
                        if (xVarMap.containsKey(absPos)) {  // when there is a C, the pos wont be in xVarMap but we still record it in thisSol
                            Map<Integer, GRBVar> ptmIdXVarMap = xVarMap.get(absPos);
                            for (int ptmId : ptmIdXVarMap.keySet()) {
                                GRBVar xVar = ptmIdXVarMap.get(ptmId);
                                if ((int) Math.round(xVar.get(GRB.DoubleAttr.Xn)) == 1) {
                                    thisSol.put(absPos, ptmId);
                                    break;
                                }
                            }
                        }
                    }
                }
                posPtmIdResList.add(thisSol);
            }
            if (lszDebugScanNum.contains(scanNum) && partSeq.contentEquals("EGPDAKL")){
                int a = 1;
            }
            if (lszDebugScanNum.contains(scanNum)){
                int a = 1;
            }
            model.dispose();
        } catch (GRBException e) {
            System.out.println("Error code: " + e.getErrorCode() + ". " + e.getMessage());
        }
        return posPtmIdResList;
    }

    public List<Map<Double, Integer>> findBestPtmMIP(int scanNum, GRBEnv env, Map<Double, List<Integer>> allMassAllPosesMap, double totalDeltaMass, String partSeq,
                                                     double ms1TolAbs, Map<Integer, Set<Double>> oneTimeMassGroups, List<Pair<Integer, Set<Double>>> multiTimeMassGroups) throws GRBException {

        List<Map<Double, Integer>> massTimeMapList = new LinkedList<>();

        try {
            GRBModel model = new GRBModel(env);
            double t = 0.01;

            //variables
            GRBLinExpr totalNumsOnPepConstr = new GRBLinExpr();
            GRBLinExpr totalMassOnPepConstr = new GRBLinExpr();


            Map<Double, GRBVar> massVarMap = new HashMap<>();
            for (double mass : allMassAllPosesMap.keySet()) {
                GRBVar tmpVar = model.addVar(0, allMassAllPosesMap.get(mass).size(), 0, GRB.INTEGER, "isPtmSelected_"+mass);
                massVarMap.put(mass, tmpVar);
                totalMassOnPepConstr.addTerm(mass, tmpVar); // + m_i * x_i
                totalNumsOnPepConstr.addTerm(1,tmpVar); // + 1 * x_i
            }
            model.addConstr(totalMassOnPepConstr, GRB.GREATER_EQUAL, totalDeltaMass - ms1TolAbs , "constrM1");
            model.addConstr(totalMassOnPepConstr, GRB.LESS_EQUAL, totalDeltaMass +ms1TolAbs, "constrM2"); //or put this to constraints as a model.addRange

            model.addConstr(totalNumsOnPepConstr, GRB.LESS_EQUAL, Math.min(partSeq.length(),3), "totalNumsOnPepConstr"); // this value should not exceed the number of aa in the partSeq

            // two extra constraints on the numebr of PTMs.
            for (int pos : oneTimeMassGroups.keySet()) {
                Set<Double> massSet = oneTimeMassGroups.get(pos);
                GRBLinExpr oneTimeMassConstr = new GRBLinExpr();
                for (double mass : massSet) {
                    oneTimeMassConstr.addTerm(1, massVarMap.get(mass));
                }
                model.addConstr(oneTimeMassConstr, GRB.LESS_EQUAL, 1, "oneTimeMassConstr_"+pos);
            }
            int dummyI = 0;
            for (Pair<Integer, Set<Double>> timeMassSetPair : multiTimeMassGroups) {
                int time = timeMassSetPair.getKey();
                Set<Double> massSet = timeMassSetPair.getValue();
                GRBLinExpr multiTimeMassConstr = new GRBLinExpr();
                for (double mass : massSet) {
                    multiTimeMassConstr.addTerm(1, massVarMap.get(mass));
                }
                model.addConstr(multiTimeMassConstr, GRB.LESS_EQUAL, time, "multiTimeMassConstr_"+time+"_"+dummyI);
                dummyI++;
            }
            //obj function
            GRBLinExpr objFunction = new GRBLinExpr();
            objFunction.addConstant(t);
            model.setObjective(totalNumsOnPepConstr, GRB.MINIMIZE); // with this, the solver will find those with smaller number of PTM first then more numbers

            Set<Double> allMassSet = new HashSet<>(allMassAllPosesMap.keySet());
//            allMassSet.add(MassTool.PROTON);
            if (lszDebugScanNum.contains(scanNum) && partSeq.contentEquals("KFGVLSDNFK")){
                int a = 1;
            }
            int poolSolutions = 100;
            model.set(GRB.IntParam.MIPFocus, 1);
            model.set(GRB.IntParam.PoolSearchMode, 1);
            model.set(GRB.IntParam.PoolSolutions, poolSolutions);
            model.set(GRB.DoubleParam.TimeLimit, 1); // second

            model.optimize();
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
//                    if (mass == MassTool.PROTON) continue; //dont record the mass of proton because it should not be placed on any aa as a ptm

                    GRBVar var = massVarMap.get(mass);
                    int time = (int) Math.round(var.get(GRB.DoubleAttr.Xn));
                    if (time > 0) {
                        massTimeMap.put(mass, time);
                    }
                }
//                if ( ! massTimeMap.isEmpty()) {   // when there is only one proton ptm, it will not be count,
                massTimeMapList.add(massTimeMap);
//                }
            }
            model.dispose();
        } catch (GRBException e) {
            System.out.println("Error code: " + e.getErrorCode() + ". " + e.getMessage());
        }
        return massTimeMapList;
    }

    public void findBestPtmMIPExtC(int scanNum, GRBEnv env, Map<Double, Set<Integer>> allMassAllPosesMap, double totalDeltaMass, int refPos, String partSeq,
                                   double ms1TolAbs, Map<Integer, Set<Double>> oneTimeMassGroups, final Map<Set<Integer>, Set<Double>> posComb_multiMassSet_Map,
                                   Map<Integer, Integer> posYIdMap, Map<Integer, Map<Double, VarPtm>> absPos_MassVarPtm_Map, List<Pair<Integer, Map<Double, Integer>>> resList) {

        Map<Integer, List<Integer>> yIdAllPosesMap = new HashMap<>(posYIdMap.values().size());
        for (int pos : posYIdMap.keySet()) {
            int yId = posYIdMap.get(pos);
            List<Integer> allPoses = yIdAllPosesMap.get(yId);
            if ( allPoses != null) {
                allPoses.add(pos);
            } else {
                allPoses = new ArrayList<>();
                allPoses.add(pos);
                yIdAllPosesMap.put(yId, allPoses);
            }
        }
        try {
            GRBModel model = new GRBModel(env);
            //// Constraints
            GRBLinExpr totalNumsOnPepConstr = new GRBLinExpr();
            GRBLinExpr totalFlexiableMass = new GRBLinExpr();
            GRBLinExpr totalFlexiableMass1 = new GRBLinExpr();
            GRBVar t = model.addVar(0, ms1TolAbs/10, 0, GRB.CONTINUOUS, "t");
            totalFlexiableMass.addTerm(-1,t);
            totalFlexiableMass1.addTerm(1,t);
            // x y variables
            Map<Integer, GRBVar> yVarMap = new HashMap<>( posYIdMap.values().size() ); // <yId, var>
            for (int yId : posYIdMap.values()) {
                GRBVar yVar = model.addVar(0, 1, 0, GRB.BINARY, "y"+yId);
                yVarMap.put(yId, yVar);

                //constraints
                double aaMass = 0;
                for (int absPos : yIdAllPosesMap.get(yId)) {
                    aaMass += massTool.getMassTable().get(partSeq.charAt(absPos-refPos));
                }
                totalFlexiableMass.addTerm(aaMass, yVarMap.get(yId));
                totalFlexiableMass1.addTerm(aaMass, yVarMap.get(yId));

            }
            List<Integer> yIdList = new ArrayList<>(yVarMap.keySet());
            yIdList.sort(Comparator.naturalOrder());
            Map<Double, GRBVar> xVarMap = new HashMap<>(allMassAllPosesMap.size());
            for (double mass : allMassAllPosesMap.keySet()) {
                Set<Integer> allAbsPoses = allMassAllPosesMap.get(mass);
                GRBVar xVar = model.addVar(0, allAbsPoses.size(), 0, GRB.INTEGER, "x_"+mass);
                xVarMap.put(mass, xVar);

                //constraints
                totalFlexiableMass.addTerm(mass, xVar); // += m_i * x_i
                totalFlexiableMass1.addTerm(mass, xVar); // += m_i * x_i

                totalNumsOnPepConstr.addTerm(1,xVar); // + 1 * x_i

                //constraints
                GRBLinExpr massOccurence = new GRBLinExpr();
                massOccurence.addTerm(1, xVar);
                int fixPosNum = 0;

                for (int absPos : allAbsPoses) {
                    if ( ! posYIdMap.containsKey(absPos)) {
                        fixPosNum++;
                        continue;
                    }

                    int yId = posYIdMap.get(absPos);
                    int position = absPos_MassVarPtm_Map.get(absPos).get(mass).position;

                    if (position == PEPC
                            && absPos-refPos != partSeq.length()-1
                            && yId != yIdList.get(yIdList.size()-1)) {
                        massOccurence.addTerm(-1, yVarMap.get(yId)); // thisYId
                        massOccurence.addTerm(1, yVarMap.get(yId+1));// nextYId
                    } else {
                        massOccurence.addTerm(-1, yVarMap.get(yId));
                    }
                }
                if (massOccurence.size() > 1) { //if it is just x <= 1, no any y, then no need to add this constraint
                    model.addConstr(massOccurence, GRB.LESS_EQUAL, fixPosNum, "massOccurence_"+mass); // xi <= fixPosNum + y1+y2+...
                }
            }

            //// add constraints
            model.addConstr(totalFlexiableMass1, GRB.GREATER_EQUAL, totalDeltaMass , "totalFlexiableMassGe");
            model.addConstr(totalFlexiableMass, GRB.LESS_EQUAL, totalDeltaMass, "totalFlexiableMassLe");
            model.addConstr(totalNumsOnPepConstr, GRB.LESS_EQUAL, Math.min(partSeq.length(), MaxPtmNumInPart), "totalNumsOnPepConstr"); // this value should not exceed the number of aa in the partSeq

            // constraints, Sum(oneMass on certain aa) < 1 or y
            for (int absPos : oneTimeMassGroups.keySet()) {
                GRBLinExpr sumX_leq_1orY = new GRBLinExpr();
//                int absPos = refPos + pos;
                for (double mass : oneTimeMassGroups.get(absPos)) {
                    sumX_leq_1orY.addTerm(1, xVarMap.get(mass));
                }
                if (posYIdMap.containsKey(absPos)){
                    sumX_leq_1orY.addTerm(-1, yVarMap.get(posYIdMap.get(absPos)));
                    model.addConstr(sumX_leq_1orY, GRB.LESS_EQUAL, 0, "sumX_leq_Y" + absPos);
                } else {
                    model.addConstr(sumX_leq_1orY, GRB.LESS_EQUAL, 1, "sumX_leq_1" + absPos);
                }
            }

            int dummyI = 0;
            for (Set<Integer> posComb : posComb_multiMassSet_Map.keySet()) {
                GRBLinExpr multiTimeMassConstr = new GRBLinExpr();
                int fixPosNum = 0;
                for (int absPos : posComb) {
                    if (posYIdMap.containsKey(absPos)) {
                        multiTimeMassConstr.addTerm(-1, yVarMap.get(posYIdMap.get(absPos)));
                    } else {
                        fixPosNum++;
                    }
                }
                Set<Double> massSet = posComb_multiMassSet_Map.get(posComb);
                for (double mass : massSet) {
                    multiTimeMassConstr.addTerm(1, xVarMap.get(mass));
                }
                model.addConstr(multiTimeMassConstr, GRB.LESS_EQUAL, fixPosNum, "multiTimeMassConstr_"+fixPosNum+"_"+dummyI);
                dummyI++;
            }

            if (yVarMap.size() >= 2) {
                for (int i = 0; i < yIdList.size()-1; i++) {
                    int thisYId = yIdList.get(i);
                    int nextYId = yIdList.get(i+1);
                    GRBLinExpr yiyj = new GRBLinExpr();
                    yiyj.addTerm(1, yVarMap.get(thisYId));
                    yiyj.addTerm(-1, yVarMap.get(nextYId));
                    model.addConstr(yiyj, GRB.GREATER_EQUAL, 0 , "yiyj_"+i);
                }
            }

            //obj function
//            model.setObjective(totalNumsOnPepConstr, GRB.MINIMIZE); // with this, the solver will find those with smaller number of PTM first then more numbers
            GRBLinExpr tConstr = new GRBLinExpr();
            tConstr.addTerm(1, t);
            model.setObjective(tConstr, GRB.MINIMIZE); // with this, the solver will find those with smaller number of PTM first then more numbers

            if (lszDebugScanNum.contains(scanNum) && partSeq.contentEquals("KFGVLSDNFK")){
                int a = 1;
            }
            model.set(GRB.IntParam.MIPFocus, 1); // 2 seems better than 1 but dont know why
            model.set(GRB.IntParam.PoolSearchMode, 2 );  //0 for only one sol, 1 for possible more but not guaranteed = poolSolutions, 2 for guaranteed but long time
            model.set(GRB.IntParam.PoolSolutions, PoolSolutions);
            model.set(GRB.DoubleParam.TimeLimit, 1); // second
            model.set(GRB.IntParam.ConcurrentMIP, 32); // second
            model.set(GRB.IntParam.Threads, 32); // second
//            model.set(GRB.DoubleParam.Heuristics, 0.05); //
//            model.set(GRB.DoubleParam.PoolGap, 1); // if the best obj is 1, then all sol with 1*(1+poolGap) with be discarded and save time

            model.optimize();
            switch (model.get(GRB.IntAttr.Status)) {
                case GRB.OPTIMAL:
                    break;
                case GRB.TIME_LIMIT:
                    break;
                default:
                    model.dispose();
                    return; // if it is disposed here (i.e. infeasible, unbounded...), dont collect solutions, just return
            }
            int solCount = model.get(GRB.IntAttr.SolCount);
            for (int solNum = 0; solNum < solCount; solNum++) {
                model.set(GRB.IntParam.SolutionNumber, solNum);

                int maxYId = -99;
                for (int yId : yVarMap.keySet()) {
                    GRBVar yVar = yVarMap.get(yId);
                    int time = (int) Math.round(yVar.get(GRB.DoubleAttr.Xn));
                    if (time > 0) {
                        if (yId > maxYId) {
                            maxYId = yId;
                        }
                    }
                }
                Map<Double, Integer> massTimeMap = new HashMap<>();
                for (double mass : xVarMap.keySet()){
                    GRBVar xVar = xVarMap.get(mass);
                    int time = (int) Math.round(xVar.get(GRB.DoubleAttr.Xn));
                    if (time > 0) {
                        massTimeMap.put(mass, time);
                    }
                }
                resList.add(new Pair<>(maxYId, massTimeMap));
            }
            model.dispose();
        } catch (GRBException e) {
            System.out.println(scanNum + "," + partSeq+ " , Error code: " + e.getErrorCode() + ". " + e.getMessage());
        }
    }

    public void findBestPtmMIPExtN(int scanNum, GRBEnv env, Map<Double, Set<Integer>> allMassAllPosesMap, double totalDeltaMass, int refPos, String partSeq,
                                   double ms1TolAbs, Map<Integer, Set<Double>> oneTimeMassGroups, final Map<Set<Integer>, Set<Double>> posComb_multiMassSet_Map,
                                   Map<Integer, Integer> posYIdMap, Map<Integer, Map<Double, VarPtm>> absPos_MassVarPtm_Map, List<Pair<Integer, Map<Double, Integer>>> resList) {

        Map<Integer, List<Integer>> yIdAllPosesMap = new HashMap<>();
        for (int pos : posYIdMap.keySet()) {
            int yId = posYIdMap.get(pos);
            List<Integer> allPoses = yIdAllPosesMap.get(yId);
            if ( allPoses != null) {
                allPoses.add(pos);
            } else {
                allPoses = new ArrayList<>();
                allPoses.add(pos);
                yIdAllPosesMap.put(yId, allPoses);
            }
        }

        try {
            GRBModel model = new GRBModel(env);
            //// Constraints
            GRBLinExpr totalNumsOnPepConstr = new GRBLinExpr();
            GRBLinExpr totalFlexiableMass = new GRBLinExpr();
            GRBLinExpr totalFlexiableMass1 = new GRBLinExpr();
            GRBVar t = model.addVar(0, ms1TolAbs, 0, GRB.CONTINUOUS, "t");
            totalFlexiableMass.addTerm(-1,t);
            totalFlexiableMass1.addTerm(1,t);
            // x y variables
            Map<Integer, GRBVar> yVarMap = new HashMap<>( posYIdMap.values().size() ); // <yId, var>
            for (int yId : posYIdMap.values()) {
                GRBVar yVar = model.addVar(0, 1, 0, GRB.BINARY, "y"+yId);
                yVarMap.put(yId, yVar);

                //constraints
                double aaMass = 0;
                for (int absPos : yIdAllPosesMap.get(yId)) {
                    aaMass += massTool.getMassTable().get(partSeq.charAt(absPos + partSeq.length()-refPos));
                }
                totalFlexiableMass.addTerm(aaMass, yVarMap.get(yId));
                totalFlexiableMass1.addTerm(aaMass, yVarMap.get(yId));

            }
            List<Integer> yIdList = new ArrayList<>(yVarMap.keySet());
            yIdList.sort(Comparator.naturalOrder());
            Map<Double, GRBVar> xVarMap = new HashMap<>(allMassAllPosesMap.size());
            for (double mass : allMassAllPosesMap.keySet()) {
                Set<Integer> allAbsPoses = allMassAllPosesMap.get(mass);
                GRBVar xVar = model.addVar(0, allAbsPoses.size(), 0, GRB.INTEGER, "x_"+mass);
                xVarMap.put(mass, xVar);

                //constraints
                totalFlexiableMass.addTerm(mass, xVar); // += m_i * x_i
                totalFlexiableMass1.addTerm(mass, xVar); // += m_i * x_i

                totalNumsOnPepConstr.addTerm(1,xVar); // + 1 * x_i

                //constraints
                GRBLinExpr massOccurence = new GRBLinExpr();
                massOccurence.addTerm(1, xVar);
                int fixPosNum = 0;

                for (int absPos : allAbsPoses) {
                    if ( ! posYIdMap.containsKey(absPos)) {
                        fixPosNum++;
                        continue;
                    }

                    int yId = posYIdMap.get(absPos);
                    int position = absPos_MassVarPtm_Map.get(absPos).get(mass).position;

                    if (position == PEPN
                            && absPos+partSeq.length() != refPos
                            && yId != yIdList.get(yIdList.size()-1)) { // NC the same
                        massOccurence.addTerm(-1, yVarMap.get(yId)); // thisYId
                        massOccurence.addTerm(1, yVarMap.get(yId+1));// nextYId // NC the same
                    } else {
                        massOccurence.addTerm(-1, yVarMap.get(yId));
                    }
                }
                if (massOccurence.size() > 1) { //if it is just x <= 1, no any y, then no need to add this constraint
                    model.addConstr(massOccurence, GRB.LESS_EQUAL, fixPosNum, "massOccurence_"+mass); // xi <= fixPosNum + y1+y2+...
                }
            }

            //// add constraints
            model.addConstr(totalFlexiableMass1, GRB.GREATER_EQUAL, totalDeltaMass , "totalFlexiableMassGe");
            model.addConstr(totalFlexiableMass, GRB.LESS_EQUAL, totalDeltaMass, "totalFlexiableMassLe");
            model.addConstr(totalNumsOnPepConstr, GRB.LESS_EQUAL, Math.min(partSeq.length(), MaxPtmNumInPart), "totalNumsOnPepConstr"); // this value should not exceed the number of aa in the partSeq

            // constraints, Sum(oneMass on certain aa) < 1 or y
            for (int absPos : oneTimeMassGroups.keySet()) {
                GRBLinExpr sumX_leq_1orY = new GRBLinExpr();
                for (double mass : oneTimeMassGroups.get(absPos)) {
                    sumX_leq_1orY.addTerm(1, xVarMap.get(mass));
                }
                if (posYIdMap.containsKey(absPos)){
                    sumX_leq_1orY.addTerm(-1, yVarMap.get(posYIdMap.get(absPos)));
                    model.addConstr(sumX_leq_1orY, GRB.LESS_EQUAL, 0, "sumX_leq_Y" + absPos);
                } else {
                    model.addConstr(sumX_leq_1orY, GRB.LESS_EQUAL, 1, "sumX_leq_1" + absPos);
                }
            }

            int dummyI = 0;
            for (Set<Integer> posComb : posComb_multiMassSet_Map.keySet()) {
                GRBLinExpr multiTimeMassConstr = new GRBLinExpr();
                int fixPosNum = 0;
                for (int absPos : posComb) {
                    if (posYIdMap.containsKey(absPos)) {
                        multiTimeMassConstr.addTerm(-1, yVarMap.get(posYIdMap.get(absPos)));
                    } else {
                        fixPosNum++;
                    }
                }
                Set<Double> massSet = posComb_multiMassSet_Map.get(posComb);
                for (double mass : massSet) {
                    multiTimeMassConstr.addTerm(1, xVarMap.get(mass));
                }
                model.addConstr(multiTimeMassConstr, GRB.LESS_EQUAL, fixPosNum, "multiTimeMassConstr_"+fixPosNum+"_"+dummyI);
                dummyI++;
            }

            if (yVarMap.size() >= 2) {
                for (int i = 0; i < yIdList.size()-1; i++) {
                    int thisYId = yIdList.get(i);
                    int nextYId = yIdList.get(i+1);
                    GRBLinExpr yiyj = new GRBLinExpr();
                    yiyj.addTerm(1, yVarMap.get(thisYId));
                    yiyj.addTerm(-1, yVarMap.get(nextYId));
                    model.addConstr(yiyj, GRB.GREATER_EQUAL, 0 , "yiyj_"+i);
                }
            }

            //obj function
            GRBLinExpr tConstr = new GRBLinExpr();
            tConstr.addTerm(1, t);
            model.setObjective(tConstr, GRB.MINIMIZE); // with this, the solver will find those with smaller number of PTM first then more numbers

            if (lszDebugScanNum.contains(scanNum) && partSeq.contentEquals("KFGVLSDNFK")){
                int a = 1;
            }
            model.set(GRB.IntParam.MIPFocus, 1); // 2 seems better than 1 but dont know why
            model.set(GRB.IntParam.PoolSearchMode, 2 );  //0 for only one sol, 1 for possible more but not guaranteed = poolSolutions, 2 for guaranteed but long time
            model.set(GRB.IntParam.PoolSolutions, PoolSolutions);
            model.set(GRB.DoubleParam.TimeLimit, 1); // second
            model.set(GRB.IntParam.ConcurrentMIP, 32); // second
            model.set(GRB.IntParam.Threads, 32); // second
//            model.set(GRB.DoubleParam.Heuristics, 0.05); //
//            model.set(GRB.DoubleParam.PoolGap, 1); // if the best obj is 1, then all sol with 1*(1+poolGap) with be discarded and save time

            model.optimize();
            switch (model.get(GRB.IntAttr.Status)) {
                case GRB.OPTIMAL:
                    break;
                case GRB.TIME_LIMIT:
                    break;
                default:
                    model.dispose();
                    return;
            }
            int solCount = model.get(GRB.IntAttr.SolCount);
            for (int solNum = 0; solNum < solCount; solNum++) {
                model.set(GRB.IntParam.SolutionNumber, solNum);

                int maxYId = -99;
                for (int yId : yVarMap.keySet()) {
                    GRBVar yVar = yVarMap.get(yId);
                    int time = (int) Math.round(yVar.get(GRB.DoubleAttr.Xn));
                    if (time > 0) {
                        if (yId > maxYId) {
                            maxYId = yId;
                        }
                    }
                }
                Map<Double, Integer> massTimeMap = new HashMap<>();
                for (double mass : xVarMap.keySet()){
                    GRBVar xVar = xVarMap.get(mass);
                    int time = (int) Math.round(xVar.get(GRB.DoubleAttr.Xn));
                    if (time > 0) {
                        massTimeMap.put(mass, time);
                    }
                }
                resList.add(new Pair<>(maxYId, massTimeMap));
            }
            model.dispose();
        } catch (GRBException e) {
            System.out.println(scanNum + "," + partSeq+ " , Error code: " + e.getErrorCode() + ". " + e.getMessage());
        }
    }

    public void updateIonMatrix(double[][] ionMatrix, double cutMass, byte ncPart){
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

    public void getFeasibleMassPosMapC(int scanNum, List<Pair<Integer, Map<Double, Integer>>> massTimeResList, TreeMap<Double, Double> plMap, String fullPartSeq,
                                       double cutMass, byte ncPart, SparseVector expProcessedPL, boolean isDecoy, Map<Double, Set<Integer>> allMassAllPosesMap,
                                       Map<Integer, Map<Double, VarPtm>> absPos_MassVarPtm_Map, TreeSet<Peptide> cModPepsSet, int startRefPos, Map<Integer, Integer> yIdMaxAbsPosMap,
                                       int optStartPos) {

        Map<Integer, Double> singleAbsPosMassMap;
        Map<Double, Integer> massTimeMap;
        for (Pair<Integer, Map<Double, Integer>> y_x_Res : massTimeResList) {
            int maxYId = y_x_Res.getFirst();
            int seqEndPos;

            if (maxYId == -99) {
                seqEndPos = optStartPos;
            } else {
                seqEndPos = yIdMaxAbsPosMap.get(maxYId)+1;
            }

            String partSeq = fullPartSeq.substring(0, seqEndPos-startRefPos);

            massTimeMap = y_x_Res.getSecond();
            singleAbsPosMassMap = new HashMap<>();
            List<Double> massMaxMultiTime = new ArrayList<>(); // the mass order in this list is fixed. should be reused to match there poses.
            for (double mass : massTimeMap.keySet()){
                Set<Integer> allTheoAbsPoses = allMassAllPosesMap.get(mass);
                if (allTheoAbsPoses.size() == 1) { // max poses num is 1. First use them to occupy some AAs
                    singleAbsPosMassMap.put(Collections.min(allTheoAbsPoses), mass); // not finding minimum just get the element
                } else {
                    for (int i = 0; i < massTimeMap.get(mass); i++) {  // if in res, one mass has xi = 2, then add it twice for correct cartesianProduct later
                        massMaxMultiTime.add(mass);
                    }                }
            } // finish checking all res mass whose MaxPosSize==1

            List<Set<Integer>> updatedFeasiblePosesList = new ArrayList<>(massMaxMultiTime.size());
            for (double mass : massMaxMultiTime) { // for mass whose MaxPosSize > 1
                Set<Integer> feasibleAbsPoses = new HashSet<>();
                for (int absPos : allMassAllPosesMap.get(mass)) {
                    if (absPos < seqEndPos && !singleAbsPosMassMap.containsKey(absPos)) {
                        feasibleAbsPoses.add(absPos);
                    }
                }
                updatedFeasiblePosesList.add(feasibleAbsPoses);
            }

            // then the updatedFeasiblePosesList is empty, suprisingly it still works to find the only-single-time PTM pattern.
            for (List<Integer> cartesianList : Sets.cartesianProduct(updatedFeasiblePosesList)){
                Set<Integer> cartesianSet = new HashSet<>(cartesianList);
                if (cartesianSet.size() < cartesianList.size()) continue; // if the cartesian list has repeated elements, i.e., more than one ptm on one pos, skip

                Peptide tmpPeptide = new Peptide(partSeq, isDecoy, massTool);
                PosMassMap posMassMap = new PosMassMap();

                for (int absPos : singleAbsPosMassMap.keySet()) {
                    int relPos = absPos - startRefPos;
                    posMassMap.put(relPos, singleAbsPosMassMap.get(absPos)); // set every the single Pos Ptms
                    tmpPeptide.posVarPtmResMap.put(relPos, absPos_MassVarPtm_Map.get(absPos).get(singleAbsPosMassMap.get(absPos)));
                }
                for (int i = 0; i < massMaxMultiTime.size(); i++) {
                    double mass = massMaxMultiTime.get(i);
                    int absPos = cartesianList.get(i);
                    int relPos = absPos - startRefPos;
                    posMassMap.put(relPos, mass); // set every the single Pos Ptms
                    tmpPeptide.posVarPtmResMap.put(relPos, absPos_MassVarPtm_Map.get(absPos).get(mass));
                }

                double[][] ionMatrix = tmpPeptide.getIonMatrix();
                if ( ! posMassMap.isEmpty()) { // if it is empty, the ptm must be just proton and ignored, so we just use the original ionMatrix to calculate score
                    tmpPeptide.setVarPTM(posMassMap);
                    ionMatrix = tmpPeptide.getIonMatrixNow();
                    updateIonMatrix(ionMatrix, cutMass, ncPart);
                }

                Set<Integer> jRange = IntStream.rangeClosed(0, partSeq.length()-1).boxed().collect(Collectors.toSet());
                double score = massTool.buildVectorAndCalXCorr(ionMatrix, 1, expProcessedPL, tmpPeptide.matchedBions, tmpPeptide.matchedYions, jRange) ;
//                if (score > 0) {
                tmpPeptide.setScore(score*(1-tmpPeptide.posVarPtmResMap.size()*0.05));
                tmpPeptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, tmpPeptide.getIonMatrix(), ms2Tolerance));
                if (cModPepsSet.size() < 2) { //max restore 2 patterns for one peptide  //todo make this avaliable to differentiate priority
                    cModPepsSet.add(tmpPeptide);
                } else if (cModPepsSet.last().compareTo(tmpPeptide) < 0) {
                    cModPepsSet.pollLast();
                    cModPepsSet.add(tmpPeptide);
                }
//                }
            }
        }
        int a = 1;
    }

    public void getFeasibleMassPosMapN(int scanNum, List<Pair<Integer, Map<Double, Integer>>> massTimeResList, TreeMap<Double, Double> plMap, String fullPartSeq,
                                       double cutMass, byte ncPart, SparseVector expProcessedPL, boolean isDecoy, Map<Double, Set<Integer>> allMassAllPosesMap,
                                       Map<Integer, Map<Double, VarPtm>> absPos_MassVarPtm_Map, TreeSet<Peptide> nModPepsSet, int startRefPos, Map<Integer, Integer> yIdMinAbsPosMap,
                                       int maxNPos) {

        int maxYId;
        int seqStartPos;
        String partSeq;
        Map<Double, Integer> massTimeMap;
        Map<Integer, Double> singleAbsPosMassMap;
        List<Double> massMaxMultiTime;
        List<Set<Integer>> updatedFeasiblePosesList;
        for (Pair<Integer, Map<Double, Integer>> y_x_Res : massTimeResList) {
            maxYId = y_x_Res.getFirst();
            if (maxYId == -99) {
                seqStartPos = maxNPos;
            } else {
                seqStartPos = yIdMinAbsPosMap.get(maxYId);
            }

            partSeq = fullPartSeq.substring(seqStartPos-startRefPos + fullPartSeq.length());

            massTimeMap = y_x_Res.getSecond();
            singleAbsPosMassMap = new HashMap<>();
            massMaxMultiTime = new ArrayList<>(); // the mass order in this list is fixed. should be reused to match there poses.
            Set<Integer> allTheoAbsPoses;
            for (double mass : massTimeMap.keySet()){
                allTheoAbsPoses = allMassAllPosesMap.get(mass);
                if (allTheoAbsPoses.size() == 1) { // max poses num is 1. First use them to occupy some AAs
                    singleAbsPosMassMap.put(Collections.min(allTheoAbsPoses), mass);
                } else {
                    for (int i = 0; i < massTimeMap.get(mass); i++) {  // if in res, one mass has xi = 2, then add it twice for correct cartesianProduct later
                        massMaxMultiTime.add(mass);
                    }                }
            } // finish checking all res mass whose MaxPosSize==1

            updatedFeasiblePosesList = new ArrayList<>(massMaxMultiTime.size());
            for (double mass : massMaxMultiTime) { // for mass whose MaxPosSize > 1
                Set<Integer> feasibleAbsPoses = new HashSet<>();
                for (int absPos : allMassAllPosesMap.get(mass)) {
                    if (absPos >= seqStartPos && !singleAbsPosMassMap.containsKey(absPos)) {
                        feasibleAbsPoses.add(absPos);
                    }
                }
                updatedFeasiblePosesList.add(feasibleAbsPoses);
            }

            // then the updatedFeasiblePosesList is empty, suprisingly it still works to find the only-single-time PTM pattern.
            for (List<Integer> cartesianList : Sets.cartesianProduct(updatedFeasiblePosesList)){
                Set<Integer> cartesianSet = new HashSet<>(cartesianList);
                if (cartesianSet.size() < cartesianList.size()) continue; // if the cartesian list has repeated elements, i.e., more than one ptm on one pos, skip

                Peptide tmpPeptide = new Peptide(partSeq, isDecoy, massTool);
                PosMassMap posMassMap = new PosMassMap();

                for (int absPos : singleAbsPosMassMap.keySet()) {
                    int relPos = absPos - startRefPos + partSeq.length();
                    posMassMap.put(relPos, singleAbsPosMassMap.get(absPos)); // set every the single Pos Ptms
                    tmpPeptide.posVarPtmResMap.put(relPos, absPos_MassVarPtm_Map.get(absPos).get(singleAbsPosMassMap.get(absPos)));
                }
                for (int i = 0; i < massMaxMultiTime.size(); i++) {
                    double mass = massMaxMultiTime.get(i);
                    int absPos = cartesianList.get(i);
                    int relPos = absPos - startRefPos + partSeq.length();
                    posMassMap.put(relPos, mass); // set every the single Pos Ptms
                    tmpPeptide.posVarPtmResMap.put(relPos, absPos_MassVarPtm_Map.get(absPos).get(mass));
                }

                double[][] ionMatrix = tmpPeptide.getIonMatrix();
                if ( ! posMassMap.isEmpty()) { // if it is empty, the ptm must be just proton and ignored, so we just use the original ionMatrix to calculate score
                    tmpPeptide.setVarPTM(posMassMap);
                    ionMatrix = tmpPeptide.getIonMatrixNow();
                    updateIonMatrix(ionMatrix, cutMass, ncPart);
                }

                Set<Integer> jRange = IntStream.rangeClosed(0, partSeq.length()-1).boxed().collect(Collectors.toSet());
                double score = massTool.buildVectorAndCalXCorr(ionMatrix, 1, expProcessedPL, tmpPeptide.matchedBions, tmpPeptide.matchedYions, jRange) ;
//                if (score > 0) {
                tmpPeptide.setScore(score*(1-tmpPeptide.posVarPtmResMap.size()*0.05));
                tmpPeptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, tmpPeptide.getIonMatrix(), ms2Tolerance));
                if (nModPepsSet.size() < 2) { //max restore 5 patterns for one peptide
                    nModPepsSet.add(tmpPeptide);
                } else if (nModPepsSet.last().compareTo(tmpPeptide) < 0) {
                    nModPepsSet.pollLast();
                    nModPepsSet.add(tmpPeptide);
                }
//                }
            }
        }
        int a = 1;
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
                PosMassMap newPosMassMap = new PosMassMap();
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
                PosMassMap newPosMassMap = new PosMassMap();
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
                PosMassMap newPtmPtn = new PosMassMap();
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
                    byte position = 0;
                    switch (positionStr) {
                        case "Protein N-term":
                            position = PROTN;
                            break;
                        case "Protein C-term":
                            position = PROTC;
                            break;
                        case "Any N-term":
                            position = PEPN;
                            break;
                        case "Any C-term":
                            position = PEPC;
                            break;
                        case "Anywhere":
                            position = ANYWHERE;
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

    private boolean hasFixMod(char aa) {
        if (fixModMap.containsKey(aa) && Math.abs(fixModMap.get(aa)) > 0.1) {
            return true;
        }
        return false;
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

    public Map<Integer, VarPtm[]> getIdxVarModMapNew(String partSeq, Set<Integer> fixModIdxes, byte isNorC_Part, byte isProtNorC_Term) {
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

        return idxVarModMap;
    }

    public Map<Integer, VarPtm[]> getPosInProtVarModMap(String partSeq, Set<Integer> fixModIdxes, byte isNorC_Part, byte isProtNorC_Term, int startPos) {
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
                    idxVarModMap.put(i+startPos, modArray);
                }
            }
        }

        return idxVarModMap;
    }

    public void prepareInfoCTerm(int scanNum, String partSeq,
                                 Map<Integer, Set<Double>> oneTimeMassSetsMap,
                                 Map<Set<Integer>, Set<Double>> posComb_multiMassSet_Map,
                                 Map<Integer, Map<Double, VarPtm>> pos_MassVarPtm_Map,
                                 Map<Double, Set<Integer>> allMassAllPosesMap,
                                 final Map<Integer, Integer> yIdMaxAbsPosMap,
                                 final boolean couldBeProtC,
                                 final int optStartPos,
                                 final int startRefPos) {

        int partSeqLen = partSeq.length();
        char aa;
        int absPos;
        for (int relPos = 0; relPos < partSeqLen; ++relPos) {
            aa = partSeq.charAt(relPos);
            absPos = relPos + startRefPos;
            if (aaWithFixModSet.contains(aa)) continue;

            List<Byte> positionsToTry = new ArrayList<>(3);
            positionsToTry.add(ANYWHERE);
            if (yIdMaxAbsPosMap.containsValue(absPos)) {
                positionsToTry.add(PEPC);
            }
            if (absPos == optStartPos-1) {
                if ((aa == 'K' || aa == 'R') || !cTermSpecific) {
                    positionsToTry.add(PEPC);
                }
            }
            if (couldBeProtC && relPos == partSeqLen-1) {
                positionsToTry.add(PROTC);
            }

            if (aaAllVarPtmMap.containsKey(aa)) {
                Map<Byte, List<VarPtm>> allVarPtmMap = aaAllVarPtmMap.get(aa);
                Map<String, VarPtm> dstMap = new HashMap<>(128);
                for (Byte position : positionsToTry) {
                    for (VarPtm varPtm : allVarPtmMap.get(position)) {
                        double mass = varPtm.mass;
                        if (massTable.get(partSeq.charAt(relPos)) + mass < ms2Tolerance) continue;
                        String varPtmStr = varPtm.getStr();

                        VarPtm oldVarPtm = dstMap.get(varPtmStr);
                        if (oldVarPtm != null) {
                            if (varPtm.priority > oldVarPtm.priority){
                                dstMap.put(varPtmStr, varPtm);
                            }
                        } else {
                            dstMap.put(varPtmStr, varPtm);
                        }

                        Set<Integer> allPoses = allMassAllPosesMap.get(mass);
                        if (allPoses != null) {
                            allPoses.add(absPos);
                        } else {
                            allPoses = new HashSet<>();
                            allPoses.add(absPos);
                            allMassAllPosesMap.put(varPtm.mass, allPoses);
                        }
                    }
                }

                if (!dstMap.isEmpty()) {
                    Map<Double, VarPtm> massVarPtmMap = new HashMap<>(dstMap.size());
                    for (VarPtm varPtm : dstMap.values()){
                        massVarPtmMap.put(varPtm.mass, varPtm);
                    }
                    pos_MassVarPtm_Map.put(absPos, massVarPtmMap);
                }
            }
        } // end for (int i = 0; i < partSeq.length(); ++i) {

        for (double mass : allMassAllPosesMap.keySet()){
            int times = allMassAllPosesMap.get(mass).size();
            if (times == 1) {
                absPos = Collections.min(allMassAllPosesMap.get(mass));

                Set<Double> massSet = oneTimeMassSetsMap.get(absPos);
                if (massSet != null) {
                    massSet.add(mass);
                } else {
                    massSet = new HashSet<>();
                    massSet.add(mass);
                    oneTimeMassSetsMap.put(absPos, massSet);
                }
            } else {
                Set<Integer> posComb = new HashSet<>(allMassAllPosesMap.get(mass));
                Set<Double> multiMassSet = posComb_multiMassSet_Map.get(posComb);
                if (multiMassSet != null) {
                    multiMassSet.add(mass);
                } else {
                    multiMassSet = new HashSet<>();
                    multiMassSet.add(mass);
                    for (int absPosTmp : posComb) {
                        if (oneTimeMassSetsMap.containsKey(absPosTmp)) {  // just think about PEPWD, pos 0 2 wont be in oneTimeMassSetsMap.keyset(), but P is in massWithMultiTime, so oneTimeMassSetsMap.get(pos) will fail.
                            multiMassSet.addAll(oneTimeMassSetsMap.get(absPosTmp));
                        }
                    }
                    posComb_multiMassSet_Map.put(posComb, multiMassSet);
                }
            }
        }
//        if (lszDebugScanNum.contains(scanNum)){
//            int a = 1;
//        }
    }

    public void prepareInfoNTerm(int scanNum, String partSeq,
                                 Map<Integer, Set<Double>> oneTimeMassSetsMap,
                                 Map<Set<Integer>, Set<Double>> posComb_multiMassSet_Map,
                                 Map<Integer, Map<Double, VarPtm>> pos_MassVarPtm_Map,
                                 Map<Double, Set<Integer>> allMassAllPosesMap,
                                 Map<Integer, Integer> yIdMinAbsPosMap,
                                 boolean couldBeProtN,
                                 int maxAbsNPos,
                                 int endRefPos,
                                 String protSeq) {

        int partSeqLen = partSeq.length();
        char aa;
        int absPos;
        for (int relPos = 0; relPos < partSeqLen; ++relPos) {
            aa = partSeq.charAt(relPos);
            absPos = relPos + endRefPos - partSeqLen;
            if (aaWithFixModSet.contains(aa) && absPos != 0) continue;  // if that pos has fix mod but is not N term, dont let it

            List<Byte> positionsToTry = new ArrayList<>(3);
            positionsToTry.add(ANYWHERE);
            if (yIdMinAbsPosMap.containsValue(absPos)) {
                positionsToTry.add(PEPN);
            }
            if (absPos == maxAbsNPos) {
                if (! nTermSpecific || absPos == 0 || protSeq.charAt(absPos-1) == 'K' || protSeq.charAt(absPos-1) == 'R') {
                    positionsToTry.add(PEPN);
                }
            }
            if (couldBeProtN && relPos == 0) {
                positionsToTry.add(PROTN);
            }

            if (aaAllVarPtmMap.containsKey(aa)) {
                Map<Byte, List<VarPtm>> allVarPtmMap = aaAllVarPtmMap.get(aa);
                Map<String, VarPtm> dstMap = new HashMap<>();
                for (Byte position : positionsToTry) {
                    for (VarPtm varPtm : allVarPtmMap.get(position)) {
                        double mass = varPtm.mass;
                        if (massTable.get(partSeq.charAt(relPos)) + mass < ms2Tolerance) continue;
                        String varPtmStr = varPtm.getStr();

                        VarPtm oldVarPtm = dstMap.get(varPtmStr);
                        if (oldVarPtm != null) {
                            if (varPtm.priority > oldVarPtm.priority){
                                dstMap.put(varPtmStr, varPtm);
                            }
                        } else {
                            dstMap.put(varPtmStr, varPtm);
                        }

                        Set<Integer> allPoses = allMassAllPosesMap.get(mass);
                        if (allPoses != null) {
                            allPoses.add(absPos);
                        } else {
                            allPoses = new HashSet<>();
                            allPoses.add(absPos);
                            allMassAllPosesMap.put(varPtm.mass, allPoses);
                        }
                    }
                }

                if (!dstMap.isEmpty()) {
                    Map<Double, VarPtm> massVarPtmMap = new HashMap<>();
                    for (VarPtm varPtm : dstMap.values()){
                        massVarPtmMap.put(varPtm.mass, varPtm);
                    }
                    pos_MassVarPtm_Map.put(absPos, massVarPtmMap);
                }
            }
        } // end for (int i = 0; i < partSeq.length(); ++i) {

        for (double mass : allMassAllPosesMap.keySet()){
            int times = allMassAllPosesMap.get(mass).size();
            if (times == 1) {
                absPos = Collections.min(allMassAllPosesMap.get(mass));

                Set<Double> massSet = oneTimeMassSetsMap.get(absPos);
                if (massSet != null) {
                    massSet.add(mass);
                } else {
                    massSet = new HashSet<>();
                    massSet.add(mass);
                    oneTimeMassSetsMap.put(absPos, massSet);
                }
            } else {
                Set<Integer> posComb = new HashSet<>(allMassAllPosesMap.get(mass));
                Set<Double> multiMassSet = posComb_multiMassSet_Map.get(posComb);
                if (multiMassSet != null) {
                    multiMassSet.add(mass);
                } else {
                    multiMassSet = new HashSet<>();
                    multiMassSet.add(mass);
                    for (int absPosTmp : posComb) {
                        if (oneTimeMassSetsMap.containsKey(absPosTmp)) {  // just think about PEPWD, pos 0 2 wont be in oneTimeMassSetsMap.keyset(), but P is in massWithMultiTime, so oneTimeMassSetsMap.get(pos) will fail.
                            multiMassSet.addAll(oneTimeMassSetsMap.get(absPosTmp));
                        }
                    }
                    posComb_multiMassSet_Map.put(posComb, multiMassSet);
                }
            }
        }
//        if (lszDebugScanNum.contains(scanNum)){
//            int a = 1;
//        }
    }

    public Map<Double, List<Integer>> getMassPosInfoMap(int scanNum, String partSeq, Set<Integer> fixModIdxes, int isNorC_Part, int isProtNorC_Term,
                                                        Map<Integer, Set<Double>> oneTimeMassSetsMap, List<Pair<Integer, Set<Double>>> multiTimeMassGroups,
                                                        Map<Integer, Map<Double, VarPtm>> pos_MassVarPtm_Map) {
//        partSeq, fixModIdxes, isNorC_Side, isProtNorC_Term
        Map<Double, List<Integer>> allMassAllPosesMap = new HashMap<>(partSeq.length() * 10, 1); //Set<pos>
        boolean hasProt_N_TermPtm = false;
        boolean hasProt_C_TermPtm = false;
        if (isNorC_Part == N_PART && isProtNorC_Term == N_TERM_PROT) {
            hasProt_N_TermPtm = true;
        }
        if (isNorC_Part == C_PART && isProtNorC_Term == C_TERM_PROT) {
            hasProt_C_TermPtm = true;
        }
        Map<Integer, Set<Double>> posMassSetMap = new HashMap<>(partSeq.length());
        for (int i = 0; i < partSeq.length(); ++i) {
            if (fixModIdxes.contains(i) && i != 0) continue;  // if that pos has fix mod but is not N term, dont let it
            char aa = partSeq.charAt(i);

            if (finalPtmMap.containsKey(aa)) {
                Map<String, VarPtm> dstMap = new HashMap<>();
                Set<Double> massesOnAa = new HashSet<>();
                List<VarPtm> srcSet = finalPtmMap.get(aa);
                if (i == 0) { //aa at seq n term
                    for (VarPtm varPtm : srcSet) {
                        if (massTable.get(partSeq.charAt(i)) + varPtm.mass < ms2Tolerance) continue;//the mass of a modified amino acid cannot be 0 or negative

                        if (varPtm.position == 4 || (varPtm.position == 2 && isNorC_Part == N_PART) || (varPtm.position == 0 && hasProt_N_TermPtm)) { // anywhere or pepN or (protN and pepPos at protN)
//                            massesOnAa.add(Double.parseDouble(df2.format(varPtm.mass)));
                            String varPtmStr = varPtm.getStr();
                            if (dstMap.containsKey(varPtmStr)) {
                                VarPtm oldVarPtm = dstMap.get(varPtmStr);
                                if (varPtm.priority > oldVarPtm.priority){
                                    dstMap.put(varPtmStr, varPtm);
                                }
                            } else {
                                dstMap.put(varPtmStr, varPtm);
                            }
                            massesOnAa.add(varPtm.mass);
                        }
                    }
                } else if (i == partSeq.length()-1) { //aa at seq c term
                    for (VarPtm varPtm : srcSet) {
                        if (massTable.get(partSeq.charAt(i)) + varPtm.mass < ms2Tolerance) continue;//the mass of a modified amino acid cannot be 0 or negative

                        if (varPtm.position == 4 || (varPtm.position == 3 && isNorC_Part == C_PART) || (varPtm.position == 1 && hasProt_C_TermPtm)) { // anywhere or pepC or (protC and pepPos at protC)
//                            massesOnAa.add(Double.parseDouble(df2.format(varPtm.mass)));
                            String varPtmStr = varPtm.getStr();
                            if (dstMap.containsKey(varPtmStr)) {
                                VarPtm oldVarPtm = dstMap.get(varPtmStr);
                                if (varPtm.priority > oldVarPtm.priority){
                                    dstMap.put(varPtmStr, varPtm);
                                }
                            } else {
                                dstMap.put(varPtmStr, varPtm);
                            }
                            massesOnAa.add(varPtm.mass);
                        }
                    }
                } else {//aa at middle
                    for (VarPtm varPtm : srcSet) {
                        if (massTable.get(partSeq.charAt(i)) + varPtm.mass < ms2Tolerance) continue;//the mass of a modified amino acid cannot be 0 or negative

                        if (varPtm.position == 4) { // anywhere
//                            massesOnAa.add(Double.parseDouble(df2.format(varPtm.mass)));
                            String varPtmStr = varPtm.getStr();
                            if (dstMap.containsKey(varPtmStr)) {
                                VarPtm oldVarPtm = dstMap.get(varPtmStr);
                                if (varPtm.priority > oldVarPtm.priority){
                                    dstMap.put(varPtmStr, varPtm);
                                }
                            } else {
                                dstMap.put(varPtmStr, varPtm);
                            }
                            massesOnAa.add(varPtm.mass);
                        }
                    }
                }
                if (!dstMap.isEmpty()) {
                    Map<Double, VarPtm> massVarPtmMap = new HashMap<>();
                    for (VarPtm varPtm : dstMap.values()){
                        massVarPtmMap.put(varPtm.mass, varPtm);
                    }
                    pos_MassVarPtm_Map.put(i, massVarPtmMap);
                }
                for (double mass : massesOnAa) {
                    if (allMassAllPosesMap.containsKey(mass)) {
                        allMassAllPosesMap.get(mass).add(i);
//                        allMassMaxTimeMap.put(mass, allMassMaxTimeMap.get(mass)+1);
                    } else {
                        List<Integer> tmpPoses = new LinkedList<>();
                        tmpPoses.add(i);
                        allMassAllPosesMap.put(mass, tmpPoses);
//                        allMassMaxTimeMap.put(mass, 1);
                    }
                }
                posMassSetMap.put(i, massesOnAa);
            }
        } // end for (int i = 0; i < partSeq.length(); ++i) {

        Set<Double> massWithMultiTime = new HashSet<>();
        for (double mass : allMassAllPosesMap.keySet()){
            int times = allMassAllPosesMap.get(mass).size();
            if (times == 1) {
                int pos = allMassAllPosesMap.get(mass).get(0);
                if (oneTimeMassSetsMap.containsKey(pos)) {
                    oneTimeMassSetsMap.get(pos).add(mass);
                } else {
                    Set<Double> massSet = new HashSet<>();
                    massSet.add(mass);
                    oneTimeMassSetsMap.put(pos, massSet);
                }
            } else {
                massWithMultiTime.add(mass);
            }
        }

        if (lszDebugScanNum.contains(scanNum)){
            int a = 1;
        }
        Map<Set<Integer>, Set<Double>> posComb_multiMassSet_Map = new HashMap<>();
        for (double mass : massWithMultiTime) {
            Set<Integer> posComb = new HashSet<>(allMassAllPosesMap.get(mass));
            if (posComb_multiMassSet_Map.containsKey(posComb)) {
                posComb_multiMassSet_Map.get(posComb).add(mass);
            } else {
                Set<Double> multiMassSet = new HashSet<>();
                multiMassSet.add(mass);
                posComb_multiMassSet_Map.put(posComb, multiMassSet);
            }
        }

        for (Set<Integer> posComb: posComb_multiMassSet_Map.keySet()) {
            Set<Double> multiMassSet = new HashSet<>();
            multiMassSet.addAll(posComb_multiMassSet_Map.get(posComb));
            for (int pos : posComb) {
                if (oneTimeMassSetsMap.containsKey(pos)) {  // just think about PEPWD, pos 0 2 wont be in oneTimeMassSetsMap.keyset(), but P is in massWithMultiTime, so oneTimeMassSetsMap.get(pos) will fail.
                    multiMassSet.addAll(oneTimeMassSetsMap.get(pos));
                }
            }
            Pair<Integer, Set<Double>> maxTime_multiMassSetPair = new Pair<>(posComb.size(), multiMassSet);
            multiTimeMassGroups.add(maxTime_multiMassSetPair);
        }
        //testlsz
//        for (double mass : massWithMultiTime) {
//            Set<Double> multiMassSet = new HashSet<>();
//            multiMassSet.add(mass);
//            for (int pos : allMassAllPosesMap.get(mass)) {
//                if (oneTimeMassSetsMap.containsKey(pos)) {  // just think about PEPWD, pos 0 2 wont be in oneTimeMassSetsMap.keyset(), but P is in massWithMultiTime, so oneTimeMassSetsMap.get(pos) will fail.
//                    multiMassSet.addAll(oneTimeMassSetsMap.get(pos));
//                }
//            }
//            Pair<Integer, Set<Double>> maxTime_multiMassSetPair = new Pair<>(allMassAllPosesMap.get(mass).size(), multiMassSet);
//            multiTimeMassGroups.add(maxTime_multiMassSetPair);
//        }
        return allMassAllPosesMap;
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
