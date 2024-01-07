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
    //    private static final int PoolSolutions = 100;
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
//    public final static DecimalFormat df2 = new DecimalFormat("0.00");
    private final MassTool massTool;
    private final Map<String, Double> elementTable;
    private final Map<Character, Double> massTable;
    private final Map<Character, Double> fixModMap;
    private Set<VarPtm> varPtmSet = new HashSet<>();
    private final double minPtmMass;
    private final double maxPtmMass;
//    public double minUserPtmMass = 0;
//    public double maxUserPtmMass = 530;
    private final double ms2Tol;
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
        this.ms2Tol = Double.valueOf(parameterMap.get("ms2_tolerance"));

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
                        GRBVar yVar = model.addVar(0, 1, -0.1, GRB.BINARY, "y"+yId);
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
                        GRBVar yVar = model.addVar(0, 1, -0.1, GRB.BINARY, "y"+yId);
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

    public void findBestPtmMIPExtCPWL(int scanNum, GRBEnv env, double totalDeltaMass, int refPos, String partSeq,
                                   double ms1TolAbs, Map<Integer, Integer> absPosYIdMap, Map<Integer, TreeMap<Double, VarPtm>> absPos_MassVarPtm_Map,
                                   int optStartPos,
                                   boolean couldBeProtC, Map<Integer, Integer> yIdMaxAbsPosMap,
                                   TreeMap<Double, Double> plMap, TreeSet<Peptide> cModPepsSet,Map<Integer, List<Byte>> absPos_ptmPositions_Map) {


        Map<Integer, List<Integer>> yIdAllPosesMap = new HashMap<>(absPosYIdMap.values().size());
        for (int pos : absPosYIdMap.keySet()) {
            int yId = absPosYIdMap.get(pos);
            List<Integer> allPoses = yIdAllPosesMap.get(yId);
            if (allPoses != null) {
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
            GRBLinExpr totalNumsOnPep = new GRBLinExpr();
            GRBLinExpr totalFlexiableMass = new GRBLinExpr();
            // y variables
            Map<Integer, GRBVar> yVarMap = new HashMap<>(yIdAllPosesMap.size()); // <yId, var>
            for (int yId : absPosYIdMap.values()) {

                GRBVar yVar = model.addVar(0, 1, -0.1, GRB.BINARY, "y" + yId);
                yVarMap.put(yId, yVar);
                //constraints
                double aaMass = 0;
                for (int absPos : yIdAllPosesMap.get(yId)) {
                    aaMass += massTool.getMassTable().get(partSeq.charAt(absPos - refPos));
                }
                totalFlexiableMass.addTerm(aaMass, yVarMap.get(yId));
            }

            // x variables
            Map<Integer, Map<Double, GRBVar>> xVarMap = new HashMap<>(); // <pos, <mass, GRBVar>>
            Map<Integer, List<Double>> xPosPEPC_MassMap = new HashMap<>();
            int partSeqLen = partSeq.length();
            for (int relPos = 0; relPos < partSeqLen; relPos++) {
                char aa = partSeq.charAt(relPos);
                int absPos = refPos + relPos;

                Map<Double, GRBVar> massXVarMap = new HashMap<>();
                xVarMap.put(absPos, massXVarMap);

//                boolean hasFixMod = aaWithFixModSet.contains(aa);
//                if (hasFixMod) continue;
                List<Byte> positionsToTry = absPos_ptmPositions_Map.get(absPos);

                GRBLinExpr sumX_leq_1orY = new GRBLinExpr();
                Map<Byte, List<VarPtm>> allVarPtmMap = aaAllVarPtmMap.get(aa);
                Set<Double> usedMassForThisAa = new HashSet<>();
                for (Byte position : positionsToTry) { // in the order of anywhere-pepC-protC
                    List<VarPtm> varPtmList = allVarPtmMap.get(position);
                    if (varPtmList == null) continue;
                    double ptmMass;
                    for (VarPtm varPtm : varPtmList) {
                        ptmMass = varPtm.mass;
                        if (massTable.get(partSeq.charAt(relPos)) + ptmMass < ms2Tol) continue;

                        if (usedMassForThisAa.contains(ptmMass)) continue;

                        GRBVar xVar = model.addVar(0, 1, -0.1, GRB.BINARY, "x" + absPos + "_" + ptmMass);
                        massXVarMap.put(ptmMass, xVar);
                        sumX_leq_1orY.addTerm(1, xVar);
                        totalFlexiableMass.addTerm(ptmMass, xVar); // + m_i * x_i
                        totalNumsOnPep.addTerm(1, xVar); // + 1 * x_i
//                        ptmId++;
                        usedMassForThisAa.add(ptmMass);
                        if (position == PEPC) {
                            List<Double> PEPC_MassList = xPosPEPC_MassMap.get(absPos);
                            if (PEPC_MassList == null) {
                                PEPC_MassList = new ArrayList<>();
                                PEPC_MassList.add(ptmMass);
                                xPosPEPC_MassMap.put(absPos, PEPC_MassList);
                            } else {
                                PEPC_MassList.add(ptmMass); // for later PEPC ptm and yi yj constraint
                            }
                        }
                    }
                }

                if (sumX_leq_1orY.size() != 0) {
                    if (absPosYIdMap.containsKey(absPos)) {
                        int yId = absPosYIdMap.get(absPos);
                        sumX_leq_1orY.addTerm(-1, yVarMap.get(yId));
                        model.addConstr(sumX_leq_1orY, GRB.LESS_EQUAL, 0, "sumX_leq_Y" + absPos);
                    } else {
                        model.addConstr(sumX_leq_1orY, GRB.LESS_EQUAL, 1, "sumX_leq_1" + absPos);
                    }
                }
            }
            // all poses x Var finished

            List<Integer> yIdList = new ArrayList<>(yVarMap.keySet());
            yIdList.sort(Comparator.naturalOrder());
            int yIdNum = yIdList.size();
            if (yIdNum >= 2) {
                for (int i = 0; i < yIdNum - 1; i++) {
                    int thisYId = yIdList.get(i);
                    int nextYId = yIdList.get(i + 1);
                    GRBLinExpr yiyj = new GRBLinExpr();
                    yiyj.addTerm(1, yVarMap.get(thisYId));
                    yiyj.addTerm(-1, yVarMap.get(nextYId));
                    model.addConstr(yiyj, GRB.GREATER_EQUAL, 0, "yiyj_" + i);

                    int thisMaxPos = yIdMaxAbsPosMap.get(thisYId);
                    if (!hasFixMod(partSeq.charAt(thisMaxPos - refPos))) {
                        for (double ptmMass : xPosPEPC_MassMap.get(thisMaxPos)) {
                            GRBLinExpr yiyjxi = new GRBLinExpr();
                            yiyjxi.addTerm(1, yVarMap.get(thisYId));
                            yiyjxi.addTerm(-1, yVarMap.get(nextYId));
                            yiyjxi.addTerm(-1, xVarMap.get(thisMaxPos).get(ptmMass));
                            model.addConstr(yiyjxi, GRB.GREATER_EQUAL, 0, "yiyjxi" + i);
                        }
                    }
                }
            } else if (yIdNum == 1) { // the PEPC ptm on the optStartPos-1 relies on the value of this yId
                int onlyOneYId = yIdList.get(0);
                if (!hasFixMod(partSeq.charAt(optStartPos - 1 - refPos))) {
                    for (double ptmMass : xPosPEPC_MassMap.get(optStartPos - 1)) {
                        GRBLinExpr yixi_1 = new GRBLinExpr();
                        yixi_1.addTerm(1, yVarMap.get(onlyOneYId));
                        yixi_1.addTerm(1, xVarMap.get(optStartPos - 1).get(ptmMass));
                        model.addConstr(yixi_1, GRB.LESS_EQUAL, 1, "yixi_1" + (optStartPos - 1));
                    }
                }
            }

            //// add constraints
            model.addConstr(totalFlexiableMass, GRB.LESS_EQUAL, totalDeltaMass + 1.2*ms1TolAbs, "constrTotalMassLE");
            model.addConstr(totalFlexiableMass, GRB.GREATER_EQUAL, totalDeltaMass - 1.2*ms1TolAbs, "constrTotalMassGE");
            model.addConstr(totalNumsOnPep, GRB.LESS_EQUAL, Math.min(partSeq.length(), MaxPtmNumInPart), "totalNumsOnPepConstr"); // this value should not exceed the number of aa in the partSeq

            // peaks
            Map<Integer, TreeMap<Double, GRBVar>> z_yIonVarMap = new HashMap<>(partSeqLen);// <tpId, <eMass, zVar>>
            int numOfTps = Math.min(partSeqLen-1, 12);
            for (int tpId = 1; tpId <= numOfTps; tpId++) { // tpId = theoretical peak Id, one less than aa num
                int absPos = refPos + tpId;
                //y ions
                double minYmz;
                double maxYmz;
                int absPosOfNextStart;
                if (absPosYIdMap.containsKey(absPos)) { // this pos is in opt zone , it has a y id
                    int thisYId = absPosYIdMap.get(absPos);
                    absPosOfNextStart = yIdMaxAbsPosMap.get(thisYId)+1;
                } else { // it does not has a y id , it is in fixed zone
                    absPosOfNextStart = optStartPos;
                }

                //for minYmz
                if (absPosYIdMap.containsKey(absPos)) {
                    minYmz = 0;
                } else {
                    int aaNumC = Math.min(MaxPtmNumInPart-1, absPosOfNextStart-absPos);
                    List<Double> bigPtmMassesOnN = new ArrayList<>(absPos - refPos);
                    for (int tmpPos = refPos; tmpPos < absPos; tmpPos++) {  //only one PTM is supposed to be here to estimate maxPtmMassOnC
                        if (absPos_MassVarPtm_Map.containsKey(tmpPos) && absPos_MassVarPtm_Map.get(tmpPos).lastKey() > 0) {
                            bigPtmMassesOnN.add(absPos_MassVarPtm_Map.get(tmpPos).lastKey());
                        } else {
                            bigPtmMassesOnN.add(0d);
                        }
                    }
                    bigPtmMassesOnN.sort(Comparator.reverseOrder());
                    double maxPtmMassOnN = 0;
                    for (int tmpPos = 0; tmpPos < Math.min(MaxPtmNumInPart-aaNumC, bigPtmMassesOnN.size()); ++tmpPos) {
                        maxPtmMassOnN += bigPtmMassesOnN.get(tmpPos);
                    }

                    List<Double> smallPtmMassesHere = new ArrayList<>(absPosOfNextStart-absPos);
                    for (int tmpPos = absPos; tmpPos < absPosOfNextStart; tmpPos++) {
                        if (absPos_MassVarPtm_Map.containsKey(tmpPos) && absPos_MassVarPtm_Map.get(tmpPos).firstKey() < 0) {
                            smallPtmMassesHere.add(absPos_MassVarPtm_Map.get(tmpPos).firstKey());
                        }else {
                            smallPtmMassesHere.add(0d);
                        }
                    }
                    smallPtmMassesHere.sort(Comparator.naturalOrder());
                    double minPtmMassHere = 0;
                    for (int i = 0; i < aaNumC; i++) { ////MaxPtmNumInPart-1 PTM here  to estimate minPtmMassHere
                        minPtmMassHere += smallPtmMassesHere.get(i);
                    }
                    double seqMass_TpId_nextStart = massTool.calResidueMass(partSeq.substring(tpId, absPosOfNextStart-refPos));
                    minYmz = Math.max(0, seqMass_TpId_nextStart + Math.max( totalDeltaMass - maxPtmMassOnN, minPtmMassHere ) );
                }

                //for maxYmz
                int aaNumC = Math.min(MaxPtmNumInPart-1, partSeqLen - tpId);
                List<Double> smallPtmMassesC = new ArrayList<>();
                double minPtmMassOnN = 0;
                for (int tmpPos = refPos; tmpPos < absPos; tmpPos++) {
                    if (absPos_MassVarPtm_Map.containsKey(tmpPos) && absPos_MassVarPtm_Map.get(tmpPos).firstKey() < 0) {
                        smallPtmMassesC.add(absPos_MassVarPtm_Map.get(tmpPos).firstKey());
                    }else {
                        smallPtmMassesC.add(0d);
                    }
                }
                smallPtmMassesC.sort(Comparator.naturalOrder());
                for (int tmpPos = 0; tmpPos < Math.min(MaxPtmNumInPart-aaNumC, smallPtmMassesC.size()); tmpPos++) {
                    minPtmMassOnN += smallPtmMassesC.get(tmpPos);
                }

                List<Double> bigPtmMassesHere = new ArrayList<>(partSeqLen - tpId);
                for (int tmpPos = absPos; tmpPos < partSeqLen+refPos; tmpPos++) {
                    if (absPos_MassVarPtm_Map.containsKey(tmpPos) && absPos_MassVarPtm_Map.get(tmpPos).lastKey() > 0) {
                        bigPtmMassesHere.add(absPos_MassVarPtm_Map.get(tmpPos).lastKey());
                    }else {
                        bigPtmMassesHere.add(0d);
                    }
                }
                bigPtmMassesHere.sort(Comparator.reverseOrder());
                double maxPtmMassHere = 0;
                for (int i = 0; i < aaNumC; i++) { ////MaxPtmNumInPart-1 PTM here  to estimate minPtmMassHere
                    maxPtmMassHere += bigPtmMassesHere.get(i);
                }

                maxYmz = Math.min(3000, massTool.calResidueMass(partSeq.substring(tpId)) + Math.min(totalDeltaMass-massTool.calResidueMass(partSeq.substring(absPosOfNextStart-refPos))-minPtmMassOnN, maxPtmMassHere) );
                // this 1000 is the max mz that I expect it should consider. No need to consider large values because the normally are not fragmented.
                if(minYmz >= maxYmz) {
                    continue;
                }

                GRBVar yMassVar = model.addVar(minYmz, maxYmz, 0, GRB.CONTINUOUS, "yMassVar_"+tpId);
                GRBLinExpr yMassDef = new GRBLinExpr();
                yMassDef.addTerm(1, yMassVar);
                for (int tmpAbsPos = absPos; tmpAbsPos < refPos+partSeqLen; tmpAbsPos++) {
                    for (double tmpMass : xVarMap.get(tmpAbsPos).keySet()){
                        yMassDef.addTerm(-tmpMass, xVarMap.get(tmpAbsPos).get(tmpMass));
                    }
                    if (absPosYIdMap.containsKey(tmpAbsPos)) {
                        yMassDef.addTerm(-massTable.get(partSeq.charAt(tmpAbsPos-refPos)), yVarMap.get(absPosYIdMap.get(tmpAbsPos)));
                    } else {
                        yMassDef.addConstant(-massTable.get(partSeq.charAt(tmpAbsPos-refPos)));
                    }
                }
                model.addConstr(yMassDef, GRB.EQUAL, MassTool.PROTON+massTool.H2O, "yMassDef_"+tpId);
                NavigableMap<Double, Double> feasiblePeaks = plMap.subMap(minYmz, true, maxYmz,true);
//                List<Map.Entry<Double, Double>> peakEntries = new ArrayList<>(feasiblePeaks.entrySet());
                int nPeaks = feasiblePeaks.size();
                if (nPeaks != 0) {
                    List<Double> pwX_yIonList = new ArrayList<>(4*nPeaks+2);
                    List<Double> pwY_yIonList = new ArrayList<>(4*nPeaks+2);
                    generatePWLpts(feasiblePeaks, pwX_yIonList, pwY_yIonList, minYmz, maxYmz);
                    double[] pwX_yIon = pwX_yIonList.stream().mapToDouble(Double::doubleValue).toArray();
                    double[] pwY_yIon = pwY_yIonList.stream().mapToDouble(Double::doubleValue).toArray();
                    model.setPWLObj(yMassVar, pwX_yIon, pwY_yIon);
                }
            }

            model.set(GRB.IntAttr.ModelSense, GRB.MAXIMIZE);
            //obj function
            model.set(GRB.DoubleParam.TimeLimit, 5); // second
            model.set(GRB.IntParam.ConcurrentMIP, 4); // second
            model.set(GRB.IntParam.Threads, 32); // second

            model.set(GRB.DoubleParam.MIPGap, 1e-1); // second
//            model.set(GRB.IntParam.CutPasses, 3); // 2: slower
//            model.set(GRB.IntParam.PrePasses, 2); // second
            model.set(GRB.IntParam.PreSparsify, 1); // second
//            model.set(GRB.DoubleParam.Heuristics, 0.5); // second

//            model.set(GRB.DoubleParam.TuneTimeLimit, 43200);
//            model.set(GRB.IntParam.TuneTrials, 5);
//            model.tune();

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
            if (lszDebugScanNum.contains(scanNum) ){ // SLEEEGAA
                int a = 1;
            }
            double optObj = model.get(GRB.DoubleAttr.ObjVal);
            double objValThres;
            if (optObj > 0) {
                objValThres = 0.8*optObj;
            } else {
                objValThres = optObj-0.2;
            }
//            double objValThres = model.get(GRB.DoubleAttr.ObjVal) - 0.2*Math.abs(model.get(GRB.DoubleAttr.ObjVal)); // becareful for negative objval
            for (int solNum = 0; solNum < solCount; solNum++) {
                model.set(GRB.IntParam.SolutionNumber, solNum);
                double poolObjVal = model.get(GRB.DoubleAttr.PoolObjVal);
                if (poolObjVal < objValThres){
                    break;
                }
                Map<Integer, Double> thisSol = new HashMap<>();
                for (int absPos = refPos; absPos < partSeq.length() + refPos; absPos++) {
                    if (absPosYIdMap.containsKey(absPos)) { //just for verification
                        int yId = absPosYIdMap.get(absPos);
                        GRBVar yVar = yVarMap.get(yId);
                        if ((int) Math.round(yVar.get(GRB.DoubleAttr.Xn)) == 1) {
                            thisSol.put(absPos, -1.0); //-1 means no ptm yet
                            if (xVarMap.containsKey(absPos)) {
                                Map<Double, GRBVar> ptmIdXVarMap = xVarMap.get(absPos);
                                for (double ptmMass : ptmIdXVarMap.keySet()) {
                                    GRBVar xVar = ptmIdXVarMap.get(ptmMass);
                                    if ((int) Math.round(xVar.get(GRB.DoubleAttr.Xn)) == 1) {
                                        thisSol.put(absPos, ptmMass);
                                        break;
                                    }
                                }
                            }
                        }
                    } else {
                        thisSol.put(absPos, -1.0); //-1 means no ptm yet
                        if (xVarMap.containsKey(absPos)) {  // when there is a C, the pos wont be in xVarMap but we still record it in thisSol
                            Map<Double, GRBVar> ptmIdXVarMap = xVarMap.get(absPos);
                            for (double ptmMass : ptmIdXVarMap.keySet()) {
                                GRBVar xVar = ptmIdXVarMap.get(ptmMass);
                                if ((int) Math.round(xVar.get(GRB.DoubleAttr.Xn)) == 1) {
                                    thisSol.put(absPos, ptmMass);
                                    break;
                                }
                            }
                        }
                    }
                }
                int trueSeqEndAbsPos = Collections.max(thisSol.keySet());
                Peptide tmpPeptide = new Peptide(partSeq.substring(0, trueSeqEndAbsPos+1-refPos), false, massTool);
                PosMassMap posMassMap = new PosMassMap();
                for (int absPos : thisSol.keySet()) {
                    if (thisSol.get(absPos) == -1) continue;
                    posMassMap.put(absPos-refPos, thisSol.get(absPos));
                    tmpPeptide.posVarPtmResMap.put(absPos-refPos, absPos_MassVarPtm_Map.get(absPos).get(thisSol.get(absPos)));
                }
                if (!posMassMap.isEmpty()) tmpPeptide.setVarPTM(posMassMap);
                tmpPeptide.setScore(poolObjVal);
                cModPepsSet.add(tmpPeptide);
            }
            model.dispose();
        } catch (GRBException e) {
            System.out.println("Error code: " + e.getErrorCode() + ". " + e.getMessage());
        }
    }

    public void findBestPtmInGapPWL(int scanNum, GRBEnv env, double totalDeltaMass, int refAbsPos, int rTagPosInProt, String gapSeq,
                                 double ms1TolAbs, Map<Integer, TreeMap<Double, VarPtm>> absPos_MassVarPtm_Map,
                                 TreeMap<Double, Double> plMap, TreeSet<Peptide> midModPepsSet, double bIonRefMass, double yIonRefMass, double precursorMass) {

        try {
            GRBModel model = new GRBModel(env);
            //// Constraints
            GRBLinExpr totalNumsInGap = new GRBLinExpr();
            GRBLinExpr deltaMass = new GRBLinExpr();

            // x variables
            Map<Integer, Map<Double, GRBVar>> xVarMap = new HashMap<>(); // <pos, <mass, GRBVar>>
            int gapLen = gapSeq.length();
            for (int relPos = 0; relPos < gapLen; relPos++) {
                char aa = gapSeq.charAt(relPos);
                int absPos = refAbsPos + relPos;
//                boolean hasFixMod = aaWithFixModSet.contains(aa);
//                if (hasFixMod) continue;

                Map<Double, GRBVar> massXVarMap = new HashMap<>();
                xVarMap.put(absPos, massXVarMap);
                GRBLinExpr sumX_leq_1 = new GRBLinExpr();
                TreeMap<Double, VarPtm> massVarPtmMap = absPos_MassVarPtm_Map.get(absPos);
                if (massVarPtmMap == null) continue;
                for (double ptmMass : massVarPtmMap.keySet()) {
                    if (massTable.get(gapSeq.charAt(relPos)) + ptmMass < ms2Tol) continue;
//                    if(Math.abs(ptmMass-68)>1 && Math.abs(ptmMass-16)>1)continue;
                    GRBVar xVar = model.addVar(0, 1, -0.1, GRB.BINARY, "x_" + absPos + "_" + ptmMass);
                    massXVarMap.put(ptmMass, xVar);
                    sumX_leq_1.addTerm(1, xVar);
                    deltaMass.addTerm(ptmMass, xVar); // + m_i * x_i
                    totalNumsInGap.addTerm(1, xVar); // + 1 * x_i
                }
                if (sumX_leq_1.size() != 0) {
                    model.addConstr(sumX_leq_1, GRB.LESS_EQUAL, 1, "sumX_leq_1" + absPos);
                }
            }// all poses x Var finished

            //// add constraints
            model.addConstr(deltaMass, GRB.LESS_EQUAL, totalDeltaMass + ms1TolAbs, "constrTotalMassLE");
            model.addConstr(deltaMass, GRB.GREATER_EQUAL, totalDeltaMass - ms1TolAbs, "constrTotalMassGE");
            model.addConstr(totalNumsInGap, GRB.LESS_EQUAL, Math.min(gapSeq.length(), MaxPtmNumInPart), "totalNumsOnPepConstr"); // this value should not exceed the number of aa in the gapSeq

            // peaks
//            Map<Integer, GRBVar> yMassVarMap = new HashMap<>(gapLen-1);
//            Map<Integer, GRBVar> bMassVarMap = new HashMap<>(gapLen-1);
            int numOfTps = Math.min(gapLen-1, 12);
            NavigableMap<Double, Double> feasiblePeaks;
            Pair<Double, Double> massRange;
            for (int tpId = 1; tpId <= numOfTps; tpId++) { // tpId = theoretical peak Id, one less than aa num
                GRBLinExpr bySumMass = new GRBLinExpr();

                int absPos = refAbsPos + tpId;
                //y ions
                massRange = getPotentialYmzRange(rTagPosInProt, absPos, absPos_MassVarPtm_Map, gapSeq, refAbsPos, totalDeltaMass);
                double minYmz = massRange.getFirst() + yIonRefMass;
                double maxYmz = massRange.getSecond() + yIonRefMass;
                if (minYmz >= maxYmz) {
                    continue;
                }
                GRBVar yMassVar = model.addVar(minYmz, maxYmz, 0, GRB.CONTINUOUS, "yMassVar_"+tpId);
                GRBLinExpr yMassDef = new GRBLinExpr();
                yMassDef.addTerm(1, yMassVar);
                for (int tmpAbsPos = absPos; tmpAbsPos < rTagPosInProt; tmpAbsPos++) {
                    for (double tmpMass : xVarMap.get(tmpAbsPos).keySet()){
                        yMassDef.addTerm(-tmpMass, xVarMap.get(tmpAbsPos).get(tmpMass));
                    }
                    yMassDef.addConstant(-massTable.get(gapSeq.charAt(tmpAbsPos-refAbsPos)));
                }
                model.addConstr(yMassDef, GRB.EQUAL, yIonRefMass, "yMassDef_"+tpId);
                feasiblePeaks = plMap.subMap(minYmz, true, maxYmz,true);
//                List<Map.Entry<Double, Double>> peakEntries = new ArrayList<>(feasiblePeaks.entrySet());
                int nPeaks = feasiblePeaks.size();
                if (nPeaks != 0) {
                    List<Double> pwX_yIonList = new ArrayList<>(4 * nPeaks + 2);
                    List<Double> pwY_yIonList = new ArrayList<>(4 * nPeaks + 2);
                    generatePWLpts(feasiblePeaks, pwX_yIonList, pwY_yIonList, minYmz, maxYmz);
                    double[] pwX_yIon = pwX_yIonList.stream().mapToDouble(Double::doubleValue).toArray();
                    double[] pwY_yIon = pwY_yIonList.stream().mapToDouble(Double::doubleValue).toArray();
                    model.setPWLObj(yMassVar, pwX_yIon, pwY_yIon);
                }
//                yMassVarMap.put(tpId, yMassVar);
                bySumMass.addTerm(1, yMassVar);
                //b ions
                massRange = getPotentialBmzRange(rTagPosInProt, absPos, absPos_MassVarPtm_Map, gapSeq, refAbsPos, totalDeltaMass);
                double minBmz = massRange.getFirst() + bIonRefMass;
                double maxBmz = massRange.getSecond() + bIonRefMass;
                if (minBmz >= maxBmz) {
                    continue;
                }
                GRBVar bMassVar = model.addVar(minBmz, maxBmz, 0, GRB.CONTINUOUS, "bMassVar_"+tpId);
                GRBLinExpr bMassDef = new GRBLinExpr();
                bMassDef.addTerm(1, bMassVar);
                for (int tmpAbsPos = refAbsPos; tmpAbsPos < absPos; tmpAbsPos++) {
                    for (double tmpMass : xVarMap.get(tmpAbsPos).keySet()){
                        bMassDef.addTerm(-tmpMass, xVarMap.get(tmpAbsPos).get(tmpMass));
                    }
                    bMassDef.addConstant(-massTable.get(gapSeq.charAt(tmpAbsPos-refAbsPos)));
                }
                model.addConstr(bMassDef, GRB.EQUAL, bIonRefMass, "bMassDef_"+tpId);
                feasiblePeaks = plMap.subMap(minBmz, true, maxBmz,true);
//                peakEntries = new ArrayList<>(feasiblePeaks.entrySet());
                nPeaks = feasiblePeaks.size();
                if (nPeaks != 0) {
                    List<Double> pwX_bIonList = new ArrayList<>(4 * nPeaks + 2);
                    List<Double> pwY_bIonList = new ArrayList<>(4 * nPeaks + 2);
                    generatePWLpts(feasiblePeaks, pwX_bIonList, pwY_bIonList, minBmz, maxBmz);
                    double[] pwX_bIon = pwX_bIonList.stream().mapToDouble(Double::doubleValue).toArray();
                    double[] pwY_bIon = pwY_bIonList.stream().mapToDouble(Double::doubleValue).toArray();
                    model.setPWLObj(bMassVar, pwX_bIon, pwY_bIon);
                }
//                bMassVarMap.put(tpId, bMassVar);
                bySumMass.addTerm(1, bMassVar);

                model.addConstr(bySumMass, GRB.LESS_EQUAL, precursorMass+2*MassTool.PROTON + ms1TolAbs, "bySumMassLEQ_"+tpId);
                model.addConstr(bySumMass, GRB.GREATER_EQUAL, precursorMass+2*MassTool.PROTON - ms1TolAbs, "bySumMassGEQ_"+tpId);
            }

//            for (int tpId = 1; tpId <= numOfTps-1; tpId++) {
//                GRBVar thisYVar = yMassVarMap.get(tpId);
//                GRBVar nextYVar = yMassVarMap.get(tpId+1);
//                GRBLinExpr thisYgtNextY = new GRBLinExpr();
//                thisYgtNextY.addTerm(1, thisYVar);
//                thisYgtNextY.addTerm(-1, nextYVar);
//                model.addConstr(thisYgtNextY, GRB.GREATER_EQUAL, 57, "thisYgtNextY_"+tpId);
//
//                GRBVar thisBVar = bMassVarMap.get(tpId);
//                GRBVar nextBVar = bMassVarMap.get(tpId+1);
//                GRBLinExpr thisBgtNextB = new GRBLinExpr();
//                thisBgtNextB.addTerm(-1, thisBVar);
//                thisBgtNextB.addTerm(1, nextBVar);
//                model.addConstr(thisBgtNextB, GRB.GREATER_EQUAL, 57, "thisBgtNextB_"+tpId);
//            }

            model.set(GRB.IntAttr.ModelSense, GRB.MAXIMIZE);
            //obj function
            model.set(GRB.DoubleParam.TimeLimit, 5); // second
            model.set(GRB.IntParam.ConcurrentMIP, 4); // second
            model.set(GRB.IntParam.Threads, 32); // second
//            model.set(GRB.IntParam.AggFill, 1000); // second

//            model.set(GRB.IntParam.MIPFocus, 3); // second
//            model.set(GRB.IntParam.CutPasses, 3); // 2: slowe r
//            model.set(GRB.IntParam.Presolve, 1); // second
            model.set(GRB.IntParam.PreSparsify, 1); // second
//            model.set(GRB.DoubleParam.Heuristics, 0.5); // second

//            model.set(GRB.DoubleParam.TuneTimeLimit, 3600);
//            model.set(GRB.IntParam.TuneTrials, 5);
//            model.tune();
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
            if (lszDebugScanNum.contains(scanNum) ){ // SLEEEGAA
                int a = 1;
            }
            double optObj = model.get(GRB.DoubleAttr.ObjVal);
            double objValThres;
            if (optObj > 0) {
                objValThres = 0.8*optObj;
            } else {
                objValThres = optObj-0.2;
            }
            for (int solNum = 0; solNum < solCount; solNum++) {
                model.set(GRB.IntParam.SolutionNumber, solNum);
                double poolObjVal = model.get(GRB.DoubleAttr.PoolObjVal);
                if (poolObjVal < objValThres){
                    break;
                }
                Map<Integer, Double> thisSol = new HashMap<>();
                for (int absPos = refAbsPos; absPos < gapSeq.length() + refAbsPos; absPos++) {
                    thisSol.put(absPos, -1.0); //-1 means no ptm yet
                    if (xVarMap.containsKey(absPos)) {  // when there is a C, the pos wont be in xVarMap but we still record it in thisSol
                        Map<Double, GRBVar> ptmIdXVarMap = xVarMap.get(absPos);
                        for (double ptmMass : ptmIdXVarMap.keySet()) {
                            GRBVar xVar = ptmIdXVarMap.get(ptmMass);
//                            System.out.println(absPos+ ", "+gapSeq.substring(absPos-refAbsPos, absPos-refAbsPos+1)+ ","+ptmMass + ","+xVar.get(GRB.DoubleAttr.Xn));
                            if ((int) Math.round(xVar.get(GRB.DoubleAttr.Xn)) == 1) {
                                thisSol.put(absPos, ptmMass);
                                break;
                            }
                        }
                    }
                }
//                System.out.println("up is sol + " + solNum+ ","+thisSol);

//                int trueSeqEndAbsPos = Collections.max(thisSol.keySet());
//                posPtmIdResSet.add(thisSol);
                Peptide tmpPeptide = new Peptide(gapSeq, false, massTool);
                PosMassMap posMassMap = new PosMassMap();
                for (int absPos : thisSol.keySet()) {
                    if (thisSol.get(absPos) == -1) continue;
                    posMassMap.put(absPos-refAbsPos, thisSol.get(absPos));
                    tmpPeptide.posVarPtmResMap.put(absPos-refAbsPos, absPos_MassVarPtm_Map.get(absPos).get(thisSol.get(absPos)));
                }
                if ( ! posMassMap.isEmpty()) {
                    tmpPeptide.setVarPTM(posMassMap);
                }
                tmpPeptide.setScore(poolObjVal);
//                try {
                midModPepsSet.add(tmpPeptide);
//                } catch (Exception e ){
//                    System.out.println(scanNum +",wrong");
//                }
            }
            model.dispose();
        } catch (GRBException e) {
            System.out.println("Error code: " + e.getErrorCode() + ". " + e.getMessage());
        }
    }

//    private void generatePWLpts(List<Map.Entry<Double, Double>> peakEntries, List<Double> pwX_List, List<Double> pwY_List){
//        int nPeaks = peakEntries.size();
//        pwX_List.set(0, 0d); pwY_List.set(0, 0d);
//        pwX_List.set(4*nPeaks+1, peakEntries.get(nPeaks-1).getKey()+1000); pwY_List.set(4*nPeaks+1, 0d);
//        for (int pId = 0; pId < nPeaks; ++pId) {
//            double tmpMz = peakEntries.get(pId).getKey();
//            double tmpIntes = peakEntries.get(pId).getValue();
//            pwX_List.set(pId*4+1, tmpMz-2*ms2Tol); pwY_List.set(pId*4+1, 0d);
//            pwX_List.set(pId*4+2, tmpMz-2*ms2Tol); pwY_List.set(pId*4+2, tmpIntes);
//            pwX_List.set(pId*4+3, tmpMz+2*ms2Tol); pwY_List.set(pId*4+3, tmpIntes);
//            pwX_List.set(pId*4+4, tmpMz+2*ms2Tol); pwY_List.set(pId*4+4, 0d);
//        }
//
//        List<Integer> pointsToDel = new LinkedList<>();
//        // adjust and delete
//        for (int pIdThis = 0; pIdThis < nPeaks-1; ++pIdThis) {
//            double thisMass = peakEntries.get(pIdThis).getKey();
//            double thisInte = peakEntries.get(pIdThis).getValue();
//            int pIdNext = pIdThis+1;
//            double nextMass = peakEntries.get(pIdNext).getKey();
//            double nextInte = peakEntries.get(pIdNext).getValue();
//            double massDiff = nextMass - thisMass;
//            if (massDiff < 4*ms2Tol) {
//                if (thisInte < nextInte) {
//                    pointsToDel.add(pIdThis*4+3);
//                    pointsToDel.add(pIdThis*4+4);
//                    pwY_List.set(pIdNext*4+1, thisInte);
//                } else if (thisInte > nextInte) {
//                    pointsToDel.add(pIdNext*4+1);
//                    pointsToDel.add(pIdNext*4+2);
//                    pwY_List.set(pIdThis*4+4, nextInte);
//                } else {
//                    pointsToDel.add(pIdThis*4+3);
//                    pointsToDel.add(pIdThis*4+4);
//                    pointsToDel.add(pIdNext*4+1);
//                    pointsToDel.add(pIdNext*4+2);
//                }
//            } else if (massDiff == 4*ms2Tol) {
//                if (thisInte != nextInte) {
//                    pointsToDel.add(pIdThis*4+4);
//                    pointsToDel.add(pIdNext*4+1);
//                } else {
//                    pointsToDel.add(pIdThis*4+3);
//                    pointsToDel.add(pIdThis*4+4);
//                    pointsToDel.add(pIdNext*4+1);
//                    pointsToDel.add(pIdNext*4+2);
//                }
//            }
//        }
//
//        pointsToDel.sort(Comparator.reverseOrder());
//        for (int delId : pointsToDel) {
//            pwX_List.remove(delId);
//            pwY_List.remove(delId);
//        }
//    }

    private void generatePWLpts(NavigableMap<Double, Double> feasiblePeaks, List<Double> pwX_List, List<Double> pwY_List, double minMz, double maxMz){
        int nPeaks = feasiblePeaks.size();
        List<Dash> dashList = new ArrayList<>(nPeaks);
        for (double mz : feasiblePeaks.keySet()) {
            dashList.add(new Dash(mz-1.1*ms2Tol, mz+1.1*ms2Tol, feasiblePeaks.get(mz)));
        }

        Set<Integer> uselessDid = new HashSet<>();
        for (int dId = 0; dId < nPeaks-1; ++dId) {
            Dash thisDash = dashList.get(dId);
            if (uselessDid.contains(dId)) continue; //useless
            for (int nextdId = dId+1; nextdId < nPeaks; ++nextdId) {
                Dash nextDash = dashList.get(nextdId);
                if (uselessDid.contains(nextdId)) continue; //useless
                if (thisDash.y < nextDash.x) {
                    break;
                } else { // two dashes overlaps or connects
                    if (thisDash.h == nextDash.h) { // two dashses same height
                        thisDash.y = Math.max(thisDash.y, nextDash.y);
                        uselessDid.add(nextdId);
//                        nextDash.y = nextDash.x; // marked as useless
                    } else if (thisDash.h < nextDash.h) {
                        if (thisDash.x < nextDash.x) {
                            thisDash.y = nextDash.x;
                        } else { //thisDash.x < nextDash.x only, thisDash.x > nextDash.x is impossible
//                            thisDash.y = thisDash.x; // marked as useless
                            uselessDid.add(dId);
                        }
                    } else { //thisDash.h > nextDash.h
                        nextDash.x = thisDash.y;
                    }

                }
            }
        }
        for (int dId : uselessDid) {
            dashList.remove(dId);
        }
        dashList.sort(Comparator.comparingDouble(o->o.x));
        int dashNum = dashList.size();
        pwX_List.add(minMz-400); pwY_List.add(0d);
        pwX_List.add(dashList.get(0).x); pwY_List.add(0d);
        pwX_List.add(dashList.get(0).x); pwY_List.add(dashList.get(0).h);
        for (int dId = 0; dId < dashNum-1; ++dId) {
            Dash thisDash = dashList.get(dId);
            Dash nextDash = dashList.get(dId+1);
            if (thisDash.y < nextDash.x) {
                pwX_List.add(thisDash.y); pwY_List.add(thisDash.h);
                pwX_List.add(thisDash.y); pwY_List.add(0d);
                pwX_List.add(nextDash.x); pwY_List.add(0d);
                pwX_List.add(nextDash.x); pwY_List.add(nextDash.h);
            } else { //thisDash.y == nextDash.x
                pwX_List.add(thisDash.y); pwY_List.add(thisDash.h);
                pwX_List.add(nextDash.x); pwY_List.add(nextDash.h);
            }
        }
        pwX_List.add(dashList.get(dashNum-1).y); pwY_List.add(dashList.get(dashNum-1).h);
        pwX_List.add(dashList.get(dashNum-1).y); pwY_List.add(0d);
        pwX_List.add(maxMz+400); pwY_List.add(0d);
    }

    private void fill_zYionVarMap(double D, double yIonRefMass, int absPos, int gapLen, int refAbsPos, NavigableMap<Double, Double> feasiblePeaks,
                                  GRBModel model, String gapSeq, Map<Integer, Map<Double, GRBVar>> xVarMap, Map<Integer, TreeMap<Double, GRBVar>> z_yIonVarMap) throws GRBException {
        TreeMap<Double, GRBVar> eMassVarMap = new TreeMap<>();
        GRBLinExpr tp_eMass_match_Base = new GRBLinExpr();
        tp_eMass_match_Base.addConstant(yIonRefMass); // yIonRefMass has contained ( MassTool.PROTON + massTool.H2O), no need to add them again
        for (int tmpPos = absPos; tmpPos < gapLen+refAbsPos; tmpPos++) {
            tp_eMass_match_Base.addConstant(massTable.get(gapSeq.charAt(tmpPos - refAbsPos)));
            Map<Double, GRBVar> ptmMassVarMap = xVarMap.get(tmpPos);
            if (ptmMassVarMap == null) continue;
            for (double ptmMass : ptmMassVarMap.keySet()) {   // Sum (xi * mi) , all needed ptm masses
                tp_eMass_match_Base.addTerm(ptmMass, ptmMassVarMap.get(ptmMass));
            }
        }
        GRBLinExpr one_tp_for_one_ep = new GRBLinExpr();
        for (double eMass : feasiblePeaks.keySet()) {
            GRBVar zVar = model.addVar(0, 1, feasiblePeaks.get(eMass), GRB.BINARY, "z_yion" + absPos + "_" + eMass); // obj coeff, directly build the obj func
            eMassVarMap.put(eMass, zVar);
            one_tp_for_one_ep.addTerm(1, zVar);
            GRBLinExpr tp_eMass_match_Actual_1 = new GRBLinExpr(tp_eMass_match_Base);
            tp_eMass_match_Actual_1.addTerm(D, zVar);
            model.addConstr(tp_eMass_match_Actual_1, GRB.LESS_EQUAL, D+ ms2Tol +eMass, "tp_eMass_match_Actual_1_yion"+absPos+"_"+eMass);

            GRBLinExpr tp_eMass_match_Actual_2 = new GRBLinExpr(tp_eMass_match_Base);
            tp_eMass_match_Actual_2.addTerm(-D, zVar);
            model.addConstr(tp_eMass_match_Actual_2, GRB.GREATER_EQUAL, -D-ms2Tol+eMass, "tp_eMass_match_Actual_2_yion"+absPos+"_"+eMass);
        }
        if (!eMassVarMap.isEmpty()) {
            z_yIonVarMap.put(absPos-refAbsPos, eMassVarMap);
            model.addConstr(one_tp_for_one_ep, GRB.LESS_EQUAL, 1, "one_tp_for_one_ep_yion"+absPos);
        }
    }

    private void fill_zBionVarMap(double D, double bIonRefMass, int absPos, int gapLen, int refAbsPos, NavigableMap<Double, Double> feasiblePeaks,
                                  GRBModel model, String gapSeq, Map<Integer, Map<Double, GRBVar>> xVarMap, Map<Integer, TreeMap<Double, GRBVar>> z_bIonVarMap) throws GRBException {
        TreeMap<Double, GRBVar> eMassVarMap = new TreeMap<>();
        GRBLinExpr tp_eMass_match_Base = new GRBLinExpr();
        tp_eMass_match_Base.addConstant(bIonRefMass); // bIonRefMass has contained ( MassTool.PROTON + massTool.H2O), no need to add them again
        for (int tmpPos = refAbsPos; tmpPos < absPos; tmpPos++) { // diff from C
            tp_eMass_match_Base.addConstant(massTable.get(gapSeq.charAt(tmpPos - refAbsPos)));
            Map<Double, GRBVar> ptmMassVarMap = xVarMap.get(tmpPos);
            if (ptmMassVarMap == null) continue;
            for (double ptmMass : ptmMassVarMap.keySet()) {   // Sum (xi * mi) , all needed ptm masses
                tp_eMass_match_Base.addTerm(ptmMass, ptmMassVarMap.get(ptmMass));
            }
        }
        GRBLinExpr one_tp_for_one_ep = new GRBLinExpr();
        for (double eMass : feasiblePeaks.keySet()) {
            GRBVar zVar = model.addVar(0, 1, feasiblePeaks.get(eMass), GRB.BINARY, "z_bion" + absPos + "_" + eMass); // obj coeff, directly build the obj func
            eMassVarMap.put(eMass, zVar);
            one_tp_for_one_ep.addTerm(1, zVar);
            GRBLinExpr tp_eMass_match_Actual_1 = new GRBLinExpr(tp_eMass_match_Base);
            tp_eMass_match_Actual_1.addTerm(D, zVar);
            model.addConstr(tp_eMass_match_Actual_1, GRB.LESS_EQUAL, D+ 2*ms2Tol +eMass, "tp_eMass_match_Actual_1_bion"+absPos+"_"+eMass);

            GRBLinExpr tp_eMass_match_Actual_2 = new GRBLinExpr(tp_eMass_match_Base);
            tp_eMass_match_Actual_2.addTerm(-D, zVar);
            model.addConstr(tp_eMass_match_Actual_2, GRB.GREATER_EQUAL, -D-2*ms2Tol+eMass, "tp_eMass_match_Actual_2_bion"+absPos+"_"+eMass);
        }
        if (!eMassVarMap.isEmpty()) {
            z_bIonVarMap.put(absPos-refAbsPos, eMassVarMap);
            model.addConstr(one_tp_for_one_ep, GRB.LESS_EQUAL, 1, "one_tp_for_one_ep_bion"+absPos);
        }
    }

    private void addZYionNoCrossConstr(int numOfTps, int gapLen, Map<Integer, TreeMap<Double, GRBVar>> z_yIonVarMap, GRBModel model) throws GRBException {
        for (int thisTpId = 1; thisTpId < numOfTps; thisTpId++) {
            if (!z_yIonVarMap.containsKey(thisTpId)) continue;
            double thisMinYmz = z_yIonVarMap.get(thisTpId).firstKey();
            int nextTpId;
            for (nextTpId = thisTpId+1; nextTpId < gapLen; nextTpId++) {
                if (!z_yIonVarMap.containsKey(nextTpId)) continue;
                double nextMaxYmz = z_yIonVarMap.get(nextTpId).lastKey();
                if (nextMaxYmz < thisMinYmz) break;
                for (double thisYmz : z_yIonVarMap.get(thisTpId).keySet()){
                    if (thisYmz > nextMaxYmz) break;
                    for (double nextYmz : z_yIonVarMap.get(nextTpId).subMap(thisYmz, true, nextMaxYmz, true).keySet()) {
                        GRBLinExpr zThisNextNoCross = new GRBLinExpr();
                        zThisNextNoCross.addTerm(1, z_yIonVarMap.get(thisTpId).get(thisYmz));
                        zThisNextNoCross.addTerm(1, z_yIonVarMap.get(nextTpId).get(nextYmz));
                        model.addConstr(zThisNextNoCross, GRB.LESS_EQUAL,1, String.format("zYionThisNextNoCross_%d_%f_%d_%f", thisTpId, thisYmz, nextTpId, nextYmz));
                    }
                }
            }
        }
    }

    private void addZBionNoCrossConstr(int numOfTps, int gapLen, Map<Integer, TreeMap<Double, GRBVar>> z_bIonVarMap, GRBModel model) throws GRBException {
        for (int thisTpId = numOfTps-1; thisTpId >= 1; thisTpId--) {
            if (!z_bIonVarMap.containsKey(thisTpId)) continue;
            double thisMinYmz = z_bIonVarMap.get(thisTpId).firstKey();
            int nextTpId;
            for (nextTpId = thisTpId-1; nextTpId >= 0; nextTpId--) {
                if (!z_bIonVarMap.containsKey(nextTpId)) continue;
                double nextMaxYmz = z_bIonVarMap.get(nextTpId).lastKey();
                if (nextMaxYmz < thisMinYmz) break;
                for (double thisYmz : z_bIonVarMap.get(thisTpId).keySet()){
                    if (thisYmz > nextMaxYmz) break;
                    for (double nextYmz : z_bIonVarMap.get(nextTpId).subMap(thisYmz, true, nextMaxYmz, true).keySet()) {
                        GRBLinExpr zThisNextNoCross = new GRBLinExpr();
                        zThisNextNoCross.addTerm(1, z_bIonVarMap.get(thisTpId).get(thisYmz));
                        zThisNextNoCross.addTerm(1, z_bIonVarMap.get(nextTpId).get(nextYmz));
                        model.addConstr(zThisNextNoCross, GRB.LESS_EQUAL,1, String.format("zBionThisNextNoCross_%d_%f_%d_%f", thisTpId, thisYmz, nextTpId, nextYmz));
                    }
                }
            }
        }
    }
    private Pair<Double, Double> getPotentialYmzRange(int rTagPosInProt, int absPos, Map<Integer, TreeMap<Double, VarPtm>> absPos_MassVarPtm_Map, String gapSeq, int refAbsPos, double totalDeltaMass) {
        double minYmz;
        double maxYmz;
        int cLen = rTagPosInProt-absPos;
        int nLen = absPos - refAbsPos;
        int ptmNumOnC = Math.min(MaxPtmNumInPart-1, cLen);
        int ptmNumOnN = Math.min(MaxPtmNumInPart-ptmNumOnC, nLen);

        // prepare small and big PTM masses for N
        List<Double> allSmallPtmsOnN = new ArrayList<>(nLen);
        List<Double> allBigPtmsOnN = new ArrayList<>(nLen);
        for (int tmpPos = refAbsPos; tmpPos < absPos; tmpPos++) {  //only one PTM is supposed to be here to estimate maxPtmMassOnC
            if (absPos_MassVarPtm_Map.containsKey(tmpPos) && absPos_MassVarPtm_Map.get(tmpPos).lastKey() > 0) {
                allBigPtmsOnN.add(absPos_MassVarPtm_Map.get(tmpPos).lastKey());
            } else {
                allBigPtmsOnN.add(0d);
            }
            if (absPos_MassVarPtm_Map.containsKey(tmpPos) && absPos_MassVarPtm_Map.get(tmpPos).firstKey() < 0) {
                allSmallPtmsOnN.add(absPos_MassVarPtm_Map.get(tmpPos).firstKey());
            }else {
                allSmallPtmsOnN.add(0d);
            }
        }
        allBigPtmsOnN.sort(Comparator.reverseOrder());
        allSmallPtmsOnN.sort(Comparator.naturalOrder());
        double maxPtmMassOnN = 0;
        double minPtmMassOnN = 0;
        for (int tmpPos = 0; tmpPos < ptmNumOnN; ++tmpPos) {
            minPtmMassOnN += allSmallPtmsOnN.get(tmpPos);
            maxPtmMassOnN += allBigPtmsOnN.get(tmpPos);
        }

        // prepare small and big PTM masses for C
        List<Double> allBigPtmsOnC = new ArrayList<>(cLen);
        List<Double> allSmallPtmsOnC = new ArrayList<>(cLen);
        double sumAaMassOnC = 0;
        for (int tmpPos = absPos; tmpPos < rTagPosInProt; tmpPos++) {
            if (absPos_MassVarPtm_Map.containsKey(tmpPos) && absPos_MassVarPtm_Map.get(tmpPos).firstKey() < 0) {
                allSmallPtmsOnC.add(absPos_MassVarPtm_Map.get(tmpPos).firstKey());
            }else {
                allSmallPtmsOnC.add(0d);
            }
            if (absPos_MassVarPtm_Map.containsKey(tmpPos) && absPos_MassVarPtm_Map.get(tmpPos).lastKey() > 0) {
                allBigPtmsOnC.add(absPos_MassVarPtm_Map.get(tmpPos).lastKey());
            }else {
                allBigPtmsOnC.add(0d);
            }
            sumAaMassOnC += massTable.get(gapSeq.charAt(tmpPos-refAbsPos));
        }
        allSmallPtmsOnC.sort(Comparator.naturalOrder());
        allBigPtmsOnC.sort(Comparator.reverseOrder());
        double minPtmMassOnC = 0;
        double maxPtmMassOnC = 0;
        for (int i = 0; i < ptmNumOnC; i++) { ////MaxPtmNumInPart-1 PTM here  to estimate minPtmMassOnC
            minPtmMassOnC += allSmallPtmsOnC.get(i);
            maxPtmMassOnC += allBigPtmsOnC.get(i);
        }
//        GRBVar bMassVar = new GRBVar();
        minYmz = Math.max(0   , sumAaMassOnC + Math.max( totalDeltaMass-maxPtmMassOnN, minPtmMassOnC ) );
        maxYmz = Math.min(3000, sumAaMassOnC + Math.min( totalDeltaMass-minPtmMassOnN, maxPtmMassOnC ) );
        // this 1000 is the max mz that I expect it should consider. No need to consider large values because the normally are not fragmented.
        return new Pair<>(minYmz, maxYmz);
    }

    private Pair<Double, Double> getPotentialBmzRange(int rTagPosInProt, int absPos, Map<Integer, TreeMap<Double, VarPtm>> absPos_MassVarPtm_Map, String gapSeq, int refAbsPos, double totalDeltaMass) {
        double minBmz;
        double maxBmz;
        int nLen = absPos - refAbsPos;
        int cLen = rTagPosInProt-absPos;
        int ptmNumOnN = Math.min(MaxPtmNumInPart-1, nLen); // diff from C
        int ptmNumOnC = Math.min(MaxPtmNumInPart-ptmNumOnN, cLen);// diff from C

        // prepare small and big PTM masses for N
        List<Double> allSmallPtmsOnN = new ArrayList<>(nLen);
        List<Double> allBigPtmsOnN = new ArrayList<>(nLen);
        double sumAaMassOnN = 0;
        for (int tmpPos = refAbsPos; tmpPos < absPos; tmpPos++) {  //only one PTM is supposed to be here to estimate maxPtmMassOnC
            if (absPos_MassVarPtm_Map.containsKey(tmpPos) && absPos_MassVarPtm_Map.get(tmpPos).lastKey() > 0) {
                allBigPtmsOnN.add(absPos_MassVarPtm_Map.get(tmpPos).lastKey());
            } else {
                allBigPtmsOnN.add(0d);
            }
            if (absPos_MassVarPtm_Map.containsKey(tmpPos) && absPos_MassVarPtm_Map.get(tmpPos).firstKey() < 0) {
                allSmallPtmsOnN.add(absPos_MassVarPtm_Map.get(tmpPos).firstKey());
            }else {
                allSmallPtmsOnN.add(0d);
            }
            sumAaMassOnN += massTable.get(gapSeq.charAt(tmpPos-refAbsPos));
        }
        allBigPtmsOnN.sort(Comparator.reverseOrder());
        allSmallPtmsOnN.sort(Comparator.naturalOrder());
        double maxPtmMassOnN = 0;
        double minPtmMassOnN = 0;
        for (int tmpPos = 0; tmpPos < ptmNumOnN; ++tmpPos) {
            minPtmMassOnN += allSmallPtmsOnN.get(tmpPos);
            maxPtmMassOnN += allBigPtmsOnN.get(tmpPos);
        }

        // prepare small and big PTM masses for C
        List<Double> allBigPtmsOnC = new ArrayList<>(cLen);
        List<Double> allSmallPtmsOnC = new ArrayList<>(cLen);
        for (int tmpPos = absPos; tmpPos < rTagPosInProt; tmpPos++) {
            if (absPos_MassVarPtm_Map.containsKey(tmpPos) && absPos_MassVarPtm_Map.get(tmpPos).firstKey() < 0) {
                allSmallPtmsOnC.add(absPos_MassVarPtm_Map.get(tmpPos).firstKey());
            }else {
                allSmallPtmsOnC.add(0d);
            }
            if (absPos_MassVarPtm_Map.containsKey(tmpPos) && absPos_MassVarPtm_Map.get(tmpPos).lastKey() > 0) {
                allBigPtmsOnC.add(absPos_MassVarPtm_Map.get(tmpPos).lastKey());
            }else {
                allBigPtmsOnC.add(0d);
            }
        }
        allSmallPtmsOnC.sort(Comparator.naturalOrder());
        allBigPtmsOnC.sort(Comparator.reverseOrder());
        double minPtmMassOnC = 0;
        double maxPtmMassOnC = 0;
        for (int i = 0; i < ptmNumOnC; i++) { ////MaxPtmNumInPart-1 PTM here  to estimate minPtmMassOnC
            minPtmMassOnC += allSmallPtmsOnC.get(i);
            maxPtmMassOnC += allBigPtmsOnC.get(i);
        }

        minBmz = Math.max(0   , sumAaMassOnN + Math.max( totalDeltaMass-maxPtmMassOnN, minPtmMassOnC ) );
        maxBmz = Math.min(3000, sumAaMassOnN + Math.min( totalDeltaMass-minPtmMassOnN, maxPtmMassOnC ) );
        // this 1000 is the max mz that I expect it should consider. No need to consider large values because the normally are not fragmented.
        return new Pair<>(minBmz, maxBmz);
    }
    public SparseVector TmpdigitizePL(TreeMap<Double, Double> plMap) {
        SparseVector digitizedPL = new SparseVector();
        for (double mz : plMap.keySet()) {
            int idx = massTool.mzToBin(mz);
            if (Math.abs(plMap.get(mz)) > 1e-6) {
                digitizedPL.put(idx, Math.max(plMap.get(mz), digitizedPL.get(idx)));
            }
        }
        return digitizedPL;
    }

    public void findBestPtmMIPExtNPWL(int scanNum, GRBEnv env, double totalDeltaMass, int refPos, String partSeq,
                                   double ms1TolAbs, Map<Integer, Integer> absPosYIdMap,
                                   Map<Integer, TreeMap<Double, VarPtm>> absPos_MassVarPtm_Map, int optEndPosP1,
                                   boolean couldBeProtN, Map<Integer, Integer> yIdMinAbsPosMap,
                                   String protSeq, TreeMap<Double, Double> plMap,
                                   TreeSet<Peptide> nModPepsSet,
                                   Map<Integer, List<Byte>> absPos_ptmPositions_Map) {

        Map<Integer, List<Integer>> yIdAllPosesMap = new HashMap<>();
        for (int pos : absPosYIdMap.keySet()) {
            int yId = absPosYIdMap.get(pos);
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
            model.set(GRB.IntAttr.ModelSense, GRB.MAXIMIZE);

            //// Constraints
            GRBLinExpr totalNumsOnPep = new GRBLinExpr();
            GRBLinExpr totalFlexiableMass = new GRBLinExpr();
            int partSeqLen = partSeq.length();
            // x y variables
            Map<Integer, GRBVar> yVarMap = new HashMap<>( yIdAllPosesMap.size() ); // <yId, var>
            for (int yId : yIdAllPosesMap.keySet()) {
                GRBVar yVar = model.addVar(0, 1, -0.1, GRB.BINARY, "y"+yId);
                yVarMap.put(yId, yVar);
                //constraints
                double aaMass = 0;
                for (int absPos : yIdAllPosesMap.get(yId)) {
                    aaMass += massTool.getMassTable().get(partSeq.charAt(absPos+partSeqLen - refPos));
                }
                totalFlexiableMass.addTerm(aaMass, yVarMap.get(yId));
            }

            // x variables
            Map<Integer, Map<Double, GRBVar>> xVarMap = new HashMap<>(); // <pos, <mass, GRBVar>>
            Map<Integer, List<Double>> xPosPEPN_MassMap = new HashMap<>();
            for (int relPos = 0; relPos < partSeqLen; relPos++) {
                char aa = partSeq.charAt(relPos);
                int absPos = refPos + relPos - partSeqLen;

//                boolean hasFixMod = aaWithFixModSet.contains(aa);
//                if (hasFixMod && absPos != 0) continue;  // if that pos has fix mod but is not N term, dont let it

                GRBLinExpr sumX_leq_1orY = new GRBLinExpr();
                Map<Byte, List<VarPtm>> allVarPtmMap = aaAllVarPtmMap.get(aa);
                Map<Double, GRBVar> massXVarMap = new HashMap<>();
                xVarMap.put(absPos, massXVarMap);
                Set<Double> usedMassForThisAa = new HashSet<>();
                List<Byte> positionsToTry = absPos_ptmPositions_Map.get(absPos);
                for (Byte position : positionsToTry) { // in the order of anywhere-pepC-protC
                    List<VarPtm> varPtmList = allVarPtmMap.get(position);
                    if (varPtmList == null) continue;

                    double ptmMass;
                    for (VarPtm varPtm : varPtmList) {
                        ptmMass = varPtm.mass;
                        if (massTable.get(partSeq.charAt(relPos)) + ptmMass < ms2Tol) continue;
                        if (usedMassForThisAa.contains(ptmMass)) continue;
//                        if (Math.abs(86-ptmMass) >3 && Math.abs(14-ptmMass) >3) continue;// debug

                        GRBVar xVar = model.addVar(0, 1, -0.1, GRB.BINARY, "x" + absPos + "_" + ptmMass);
                        massXVarMap.put(ptmMass, xVar);
                        sumX_leq_1orY.addTerm(1, xVar);
                        totalFlexiableMass.addTerm(ptmMass, xVar); // + m_i * x_i
                        totalNumsOnPep.addTerm(1, xVar); // + 1 * x_i

                        usedMassForThisAa.add(ptmMass);
                        if (position == PEPN) {
                            List<Double> PEPN_MassList = xPosPEPN_MassMap.get(absPos);
                            if (PEPN_MassList == null) {
                                PEPN_MassList = new ArrayList<>();
                                PEPN_MassList.add(ptmMass);
                                xPosPEPN_MassMap.put(absPos, PEPN_MassList);
                            } else {
                                PEPN_MassList.add(ptmMass); // for later PEPC ptm and yi yj constraint
                            }
                        }
                    }
                }
                if (sumX_leq_1orY.size() != 0) {
                    if (absPosYIdMap.containsKey(absPos)) {
                        int yId = absPosYIdMap.get(absPos);
                        sumX_leq_1orY.addTerm(-1, yVarMap.get(yId));
                        model.addConstr(sumX_leq_1orY, GRB.LESS_EQUAL, 0, "sumX_leq_Y" + absPos);
                    } else {
                        model.addConstr(sumX_leq_1orY, GRB.LESS_EQUAL, 1, "sumX_leq_1" + absPos);
                    }
                }
            }
            // all poses x Var finished

            List<Integer> yIdList = new ArrayList<>(yVarMap.keySet());
            yIdList.sort(Comparator.naturalOrder());
            int yIdNum = yIdList.size();
            if (yIdNum >= 2) {
                for (int i = 0; i < yIdNum - 1; i++) {
                    int thisYId = yIdList.get(i);
                    int nextYId = yIdList.get(i + 1);
                    GRBLinExpr yiyj = new GRBLinExpr();
                    yiyj.addTerm(1, yVarMap.get(thisYId));
                    yiyj.addTerm(-1, yVarMap.get(nextYId));
                    model.addConstr(yiyj, GRB.GREATER_EQUAL, 0, "yiyj_" + i);

                    int thisMinPos = yIdMinAbsPosMap.get(thisYId);
                    if (!hasFixMod(partSeq.charAt(thisMinPos+partSeqLen - refPos))) {
                        for (double ptmMass : xPosPEPN_MassMap.get(thisMinPos)) {
                            GRBLinExpr yiyjxi = new GRBLinExpr();
                            yiyjxi.addTerm(1, yVarMap.get(thisYId));
                            yiyjxi.addTerm(-1, yVarMap.get(nextYId));
                            yiyjxi.addTerm(-1, xVarMap.get(thisMinPos).get(ptmMass));
                            model.addConstr(yiyjxi, GRB.GREATER_EQUAL, 0, "yiyjxi" + i);
                        }
                    }
                }
            } else if (yIdNum == 1) { // the PEPC ptm on the optStartPos-1 relies on the value of this yId
                int onlyOneYId = yIdList.get(0);
                if (!hasFixMod(partSeq.charAt(optEndPosP1+partSeqLen - refPos))) {
                    for (double ptmMass : xPosPEPN_MassMap.get(optEndPosP1)) {
                        GRBLinExpr yixi_1 = new GRBLinExpr();
                        yixi_1.addTerm(1, yVarMap.get(onlyOneYId));
                        yixi_1.addTerm(1, xVarMap.get(optEndPosP1).get(ptmMass));
                        model.addConstr(yixi_1, GRB.LESS_EQUAL, 1, "yixi_1" + (optEndPosP1));
                    }
                }
            }

            //// add constraints
            model.addConstr(totalFlexiableMass, GRB.LESS_EQUAL, totalDeltaMass + 1.2*ms1TolAbs, "constrTotalMassLE");
            model.addConstr(totalFlexiableMass, GRB.GREATER_EQUAL, totalDeltaMass - 1.2*ms1TolAbs, "constrTotalMassGE");
            model.addConstr(totalNumsOnPep, GRB.LESS_EQUAL, Math.min(partSeqLen, MaxPtmNumInPart), "totalNumsOnPep"); // this value should not exceed the number of aa in the partSeq

            // peaks
//            Map<Integer, TreeMap<Double, GRBVar>> z_bIonVarMap = new HashMap<>(partSeqLen);// <tpId, <eMass, zVar>>
            int numOfTps = Math.min(partSeqLen-1, 12);
            for (int tpId = 1; tpId <= numOfTps; tpId++) { // tpId = theoretical peak Id, one less than aa num
                int absPos = refPos + tpId - partSeqLen;
                //b ions
                double minBmz;
                double maxBmz;
                int absPosOfNextEnd;
                if (absPosYIdMap.containsKey(absPos)) { // this pos is in opt zone , it has a y id
                    int thisYId = absPosYIdMap.get(absPos);
                    absPosOfNextEnd = yIdMinAbsPosMap.get(thisYId)-1;
                } else { // it does not has a y id , it is in fixed zone
                    absPosOfNextEnd = optEndPosP1-1;
                }
                //for minBmz

                if (absPosYIdMap.containsKey(absPos) || absPosYIdMap.containsKey(absPos-1)) {
                    minBmz = 0;// acutually it should be proton mass but no difference todo
                } else {
                    int aaNumN = Math.min(MaxPtmNumInPart-1, absPos-absPosOfNextEnd-1);
                    List<Double> bigPtmMassesOnC = new ArrayList<>(refPos - absPos -1);
                    for (int tmpPos = absPos; tmpPos < refPos; tmpPos++) {  //only one PTM is supposed to be here to estimate maxPtmMassOnC
                        if (absPos_MassVarPtm_Map.containsKey(tmpPos) && absPos_MassVarPtm_Map.get(tmpPos).lastKey() > 0) {
                            bigPtmMassesOnC.add(absPos_MassVarPtm_Map.get(tmpPos).lastKey());
                        }else {
                            bigPtmMassesOnC.add(0d);
                        }
                    }
                    bigPtmMassesOnC.sort(Comparator.reverseOrder());
                    double maxPtmMassOnC = 0;
                    for (int tmpPos = 0; tmpPos < Math.min(MaxPtmNumInPart-aaNumN, bigPtmMassesOnC.size()); ++tmpPos) {
                        maxPtmMassOnC += bigPtmMassesOnC.get(tmpPos);
                    }

                    List<Double> smallPtmMassesHere = new ArrayList<>(absPos - absPosOfNextEnd);
                    for (int tmpPos = absPosOfNextEnd+1; tmpPos < absPos; tmpPos++) {
                        if (absPos_MassVarPtm_Map.containsKey(tmpPos) && absPos_MassVarPtm_Map.get(tmpPos).firstKey() < 0) {
                            smallPtmMassesHere.add(absPos_MassVarPtm_Map.get(tmpPos).firstKey());
                        }else {
                            smallPtmMassesHere.add(0d);
                        }
                    }
                    smallPtmMassesHere.sort(Comparator.naturalOrder());
                    double minPtmMassHere = 0;
                    for (int i = 0; i < aaNumN; i++) { ////MaxPtmNumInPart-1 PTM here  to estimate minPtmMassHere
                        minPtmMassHere += smallPtmMassesHere.get(i);
                    }
                    double seqMass_nextEnd_TpId = massTool.calResidueMass(partSeq.substring(absPosOfNextEnd+1-refPos+partSeqLen, tpId));
                    minBmz = Math.max(0, seqMass_nextEnd_TpId + Math.max( totalDeltaMass - maxPtmMassOnC, minPtmMassHere ) );
                }

                //for maxBmz
                int aaNumN = Math.min(MaxPtmNumInPart-1, tpId+1);
                List<Double> smallPtmMassesC = new ArrayList<>();
                double minPtmMassOnC = 0;
                for (int tmpPos = absPos+1; tmpPos < refPos; tmpPos++) {
                    if (absPos_MassVarPtm_Map.containsKey(tmpPos) && absPos_MassVarPtm_Map.get(tmpPos).firstKey() < 0) {
                        smallPtmMassesC.add(absPos_MassVarPtm_Map.get(tmpPos).firstKey());
                    }else {
                        smallPtmMassesC.add(0d);
                    }
                }
                smallPtmMassesC.sort(Comparator.naturalOrder());
                for (int tmpPos = 0; tmpPos < Math.min(MaxPtmNumInPart-aaNumN, smallPtmMassesC.size()); tmpPos++) {
                    minPtmMassOnC += smallPtmMassesC.get(tmpPos);
                }
                List<Double> bigPtmMassesHere = new ArrayList<>(tpId+1);
                for (int tmpPos = refPos-partSeqLen; tmpPos < absPos+1; tmpPos++) {
                    if (absPos_MassVarPtm_Map.containsKey(tmpPos) && absPos_MassVarPtm_Map.get(tmpPos).lastKey() > 0) {
                        bigPtmMassesHere.add(absPos_MassVarPtm_Map.get(tmpPos).lastKey());
                    }else {
                        bigPtmMassesHere.add(0d);
                    }
                }
                bigPtmMassesHere.sort(Comparator.reverseOrder());
                double maxPtmMassHere = 0;
                for (int i = 0; i < aaNumN; i++) { ////MaxPtmNumInPart-1 PTM here  to estimate minPtmMassHere
                    maxPtmMassHere += bigPtmMassesHere.get(i);
                }
                maxBmz = Math.min(3000, massTool.calResidueMass(partSeq.substring(0, tpId+1)) + Math.min(totalDeltaMass-massTool.calResidueMass(partSeq.substring(0, absPosOfNextEnd+1-refPos+partSeqLen))-minPtmMassOnC, maxPtmMassHere) );
                // this 1000 is the max mz that I expect it should consider. No need to consider large values because the normally are not fragmented.
                if (minBmz > maxBmz) {
//                    System.out.println(scanNum + ",minBmz > maxBmz,");
                    continue;
                }

                GRBVar bMassVar = model.addVar(minBmz, maxBmz, 0, GRB.CONTINUOUS, "bMassVar_"+tpId);
                GRBLinExpr bMassDef = new GRBLinExpr();
                bMassDef.addTerm(1, bMassVar);
                for (int tmpAbsPos = refPos-partSeqLen; tmpAbsPos < absPos; tmpAbsPos++) {
                    for (double tmpMass : xVarMap.get(tmpAbsPos).keySet()){
                        bMassDef.addTerm(-tmpMass, xVarMap.get(tmpAbsPos).get(tmpMass));
                    }
                    if (absPosYIdMap.containsKey(tmpAbsPos)) {
                        bMassDef.addTerm(-massTable.get(partSeq.charAt(tmpAbsPos-refPos+partSeqLen)), yVarMap.get(absPosYIdMap.get(tmpAbsPos)));
                    } else {
                        bMassDef.addConstant(-massTable.get(partSeq.charAt(tmpAbsPos-refPos+partSeqLen)));
                    }
                }
                model.addConstr(bMassDef, GRB.EQUAL, MassTool.PROTON, "bMassDef_"+tpId);
                NavigableMap<Double, Double> feasiblePeaks = plMap.subMap(minBmz, true, maxBmz,true);
                int nPeaks = feasiblePeaks.size();
                if (nPeaks != 0) {
                    List<Double> pwX_bIonList = new ArrayList<>(4*nPeaks+2);
                    List<Double> pwY_bIonList = new ArrayList<>(4*nPeaks+2);
                    generatePWLpts(feasiblePeaks, pwX_bIonList, pwY_bIonList, minBmz, maxBmz);
                    double[] pwX_bIon = pwX_bIonList.stream().mapToDouble(Double::doubleValue).toArray();
                    double[] pwY_bIon = pwY_bIonList.stream().mapToDouble(Double::doubleValue).toArray();
                    model.setPWLObj(bMassVar, pwX_bIon, pwY_bIon);
                }
//                bySumMass.addTerm(1, bMassVar);
            }

            //obj function
            model.set(GRB.DoubleParam.TimeLimit, 5); // timeLimit
            model.set(GRB.IntParam.ConcurrentMIP, 4); // second
            model.set(GRB.IntParam.Threads, 32); // second

//            model.set(GRB.DoubleParam.MIPGap, 1e-1); // second
//            model.set(GRB.IntParam.CutPasses, 3); // 2: slower
//            model.set(GRB.IntParam.PrePasses, 2); // second
            model.set(GRB.IntParam.PreSparsify, 1); // second
//            model.set(GRB.DoubleParam.Heuristics, 0.5); // second

//            model.set(GRB.IntParam.AggFill, 1000); // second

//            model.set(GRB.DoubleParam.TuneTimeLimit, 3600);//43200
//            model.set(GRB.IntParam.TuneTrials, 5);
////            model.set(GRB.IntParam.TuneCriterion, 3);
//            model.tune();

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
            double optObj = model.get(GRB.DoubleAttr.ObjVal);
            double objValThres;
            if (optObj > 0) {
                objValThres = 0.8*optObj;
            } else {
                objValThres = optObj-0.2;
            }
//            double objValThres = model.get(GRB.DoubleAttr.ObjVal) - 0.2*Math.abs(model.get(GRB.DoubleAttr.ObjVal)); // becareful for negative objval
            for (int solNum = 0; solNum < solCount; solNum++) {
                model.set(GRB.IntParam.SolutionNumber, solNum);
                double poolObjVal = model.get(GRB.DoubleAttr.PoolObjVal);
                if (model.get(GRB.DoubleAttr.PoolObjVal) < objValThres){
                    break;
                }
                Map<Integer, Double> thisSol = new HashMap<>();
                for (int absPos = refPos-partSeqLen; absPos < refPos; absPos++) {
                    if (absPosYIdMap.containsKey(absPos)) { //just for verification
                        int yId = absPosYIdMap.get(absPos);
                        GRBVar yVar = yVarMap.get(yId);
//                        if (yVar == null) {
//                            int a = 1;
//                        }
                        if ((int) Math.round(yVar.get(GRB.DoubleAttr.Xn)) == 1) {
                            thisSol.put(absPos, -1.0); //-1 means no ptm yet
                            if (xVarMap.containsKey(absPos)) {
                                Map<Double, GRBVar> ptmIdXVarMap = xVarMap.get(absPos);
                                for (double ptmMass : ptmIdXVarMap.keySet()) {
                                    GRBVar xVar = ptmIdXVarMap.get(ptmMass);
                                    if ((int) Math.round(xVar.get(GRB.DoubleAttr.Xn)) == 1) {
                                        thisSol.put(absPos, ptmMass);
                                        break;
                                    }
                                }
                            }
                        }
                    } else {
                        thisSol.put(absPos, -1.0); //-1 means no ptm yet
                        if (xVarMap.containsKey(absPos)) {  // when there is a C, the pos wont be in xVarMap but we still record it in thisSol
                            Map<Double, GRBVar> ptmIdXVarMap = xVarMap.get(absPos);
                            for (double ptmMass : ptmIdXVarMap.keySet()) {
                                GRBVar xVar = ptmIdXVarMap.get(ptmMass);
                                if ((int) Math.round(xVar.get(GRB.DoubleAttr.Xn)) == 1) {
                                    thisSol.put(absPos, ptmMass);
                                    break;
                                }
                            }
                        }
                    }
                }
//                nModPepsSet
//                posPtmIdResSet.add(thisSol);

                int trueSeqStartAbsPos = Collections.min(thisSol.keySet());
                Peptide tmpPeptide = new Peptide(partSeq.substring(trueSeqStartAbsPos-refPos+partSeqLen), false, massTool);
                PosMassMap posMassMap = new PosMassMap();
                for (int absPos : thisSol.keySet()) {
                    if (thisSol.get(absPos) == -1) continue;
                    posMassMap.put(absPos-trueSeqStartAbsPos, thisSol.get(absPos));
                    tmpPeptide.posVarPtmResMap.put(absPos-trueSeqStartAbsPos, absPos_MassVarPtm_Map.get(absPos).get(thisSol.get(absPos)));
                }
                if (!posMassMap.isEmpty()) tmpPeptide.setVarPTM(posMassMap);
                tmpPeptide.setScore(poolObjVal);
                nModPepsSet.add(tmpPeptide);
                if (java.lang.management.ManagementFactory.getRuntimeMXBean().getInputArguments().toString().indexOf("jdwp") >= 0){
//                    objValThres = -10;
                    double[][] ions = tmpPeptide.getIonMatrixNow();
                    int a = 1;
                }
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

    public void getFeasibleMassPosMapC(int scanNum, Set<Pair<Integer, Map<Double, Integer>>> massTimeResList, TreeMap<Double, Double> plMap, String fullPartSeq,
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

            Set<PosMassMap> triedPtmPattern = new HashSet<>(200);
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

                if (triedPtmPattern.contains(posMassMap)) {
                    continue;
                }else {
                    triedPtmPattern.add(posMassMap);
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
                tmpPeptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, tmpPeptide.getIonMatrix(), ms2Tol));
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

    public void getFeasibleMassPosMapN(int scanNum, Set<Pair<Integer, Map<Double, Integer>>> massTimeResList, TreeMap<Double, Double> plMap, String fullPartSeq,
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

            Set<PosMassMap> triedPtmPattern = new HashSet<>(200);
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

                if (triedPtmPattern.contains(posMassMap)) {
                    continue;
                }else {
                    triedPtmPattern.add(posMassMap);
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
                tmpPeptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, tmpPeptide.getIonMatrix(), ms2Tol));
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
                    peptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, peptide.getIonMatrix(), ms2Tol));
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
                    peptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, peptide.getIonMatrix(), ms2Tol));
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
                    peptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, peptide.getIonMatrix(), ms2Tol));
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
                if (Math.abs(mass - 79.957) < 0.005) {
                    continue;
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
                                if ((leftFlank != 'K' || Math.abs(massTable.get('K') - modEntry.mass) > ms2Tol) && (leftFlank != 'R' || Math.abs(massTable.get('R') - modEntry.mass) > ms2Tol) && (massTable.get(ptmFreePeptide.charAt(1)) + modEntry.mass > ms2Tol)) {  // Fixing missed cleavages caused issue in N-term and the mass of a modified amino acid cannot be 0 or negative.
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
                                if ((rightFlank != 'K' || Math.abs(massTable.get('K') - modEntry.mass) > ms2Tol) && (rightFlank != 'R' || Math.abs(massTable.get('R') - modEntry.mass) > ms2Tol) && (massTable.get(ptmFreePeptide.charAt(ptmFreePeptide.length() - 2)) + modEntry.mass > ms2Tol)) {  // Fixing missed cleavages caused issue in C-term and the mass of a modified amino acid cannot be 0 or negative
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
                        if (massTable.get(partSeq.charAt(i)) + varPtm.mass < ms2Tol) continue;//the mass of a modified amino acid cannot be 0 or negative

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
                        if (massTable.get(partSeq.charAt(i)) + varPtm.mass < ms2Tol) continue;//the mass of a modified amino acid cannot be 0 or negative

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
                        if (massTable.get(partSeq.charAt(i)) + varPtm.mass < ms2Tol) continue;//the mass of a modified amino acid cannot be 0 or negative

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
                        if (massTable.get(partSeq.charAt(i)) + varPtm.mass < ms2Tol) continue;//the mass of a modified amino acid cannot be 0 or negative

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
                        if (massTable.get(partSeq.charAt(i)) + varPtm.mass < ms2Tol) continue;//the mass of a modified amino acid cannot be 0 or negative

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
                        if (massTable.get(partSeq.charAt(i)) + varPtm.mass < ms2Tol) continue;//the mass of a modified amino acid cannot be 0 or negative

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
                                 Map<Integer, TreeMap<Double, VarPtm>> pos_MassVarPtm_Map,
                                 final Map<Integer, Integer> yIdMaxAbsPosMap,
                                 final boolean couldBeProtC,
                                 final int optStartPos,
                                 final int startRefPos,
                                 Map<Integer, List<Byte>> absPos_ptmPositions_Map) {

        int partSeqLen = partSeq.length();
        char aa;
        int absPos;
        for (int relPos = 0; relPos < partSeqLen; ++relPos) {
            aa = partSeq.charAt(relPos);
            absPos = relPos + startRefPos;
            List<Byte> positionsToTry = new ArrayList<>(3);
            absPos_ptmPositions_Map.put(absPos, positionsToTry);
            if (aaWithFixModSet.contains(aa)) continue;

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
                    if ( ! allVarPtmMap.containsKey(position)) continue;
                    for (VarPtm varPtm : allVarPtmMap.get(position)) {
                        double mass = varPtm.mass;
                        if (massTable.get(partSeq.charAt(relPos)) + mass < ms2Tol) continue;
                        String varPtmStr = varPtm.getStr();

                        VarPtm oldVarPtm = dstMap.get(varPtmStr);
                        if (oldVarPtm != null) {
                            if (varPtm.priority > oldVarPtm.priority){
                                dstMap.put(varPtmStr, varPtm);
                            }
                        } else {
                            dstMap.put(varPtmStr, varPtm);
                        }
                    }
                }
                if (!dstMap.isEmpty()) {
                    TreeMap<Double, VarPtm> massVarPtmMap = new TreeMap<>();
                    for (VarPtm varPtm : dstMap.values()){
                        massVarPtmMap.put(varPtm.mass, varPtm);
                    }
                    pos_MassVarPtm_Map.put(absPos, massVarPtmMap);
                }
            }
        } // end for (int i = 0; i < partSeq.length(); ++i) {
    }

    public void prepareInfoMid(int scanNum, String partSeq, Map<Integer, TreeMap<Double, VarPtm>> pos_MassVarPtm_Map, final int startRefPos) {

        int partSeqLen = partSeq.length();
        char aa;
        int absPos;
        for (int relPos = 0; relPos < partSeqLen; ++relPos) {
            aa = partSeq.charAt(relPos);
            absPos = relPos + startRefPos;
            if (aaWithFixModSet.contains(aa)) continue;
            if (aaAllVarPtmMap.containsKey(aa)) {
                Map<Byte, List<VarPtm>> allVarPtmMap = aaAllVarPtmMap.get(aa);
                Map<String, VarPtm> dstMap = new HashMap<>(128);
                if ( ! allVarPtmMap.containsKey(ANYWHERE)) continue;
                for (VarPtm varPtm : allVarPtmMap.get(ANYWHERE)) {
                    double mass = varPtm.mass;
                    if (massTable.get(partSeq.charAt(relPos)) + mass < ms2Tol) continue;
                    String varPtmStr = varPtm.getStr();

                    VarPtm oldVarPtm = dstMap.get(varPtmStr);
                    if (oldVarPtm != null) {
                        if (varPtm.priority > oldVarPtm.priority){
                            dstMap.put(varPtmStr, varPtm);
                        }
                    } else {
                        dstMap.put(varPtmStr, varPtm);
                    }
                }
                if (!dstMap.isEmpty()) {
                    TreeMap<Double, VarPtm> massVarPtmMap = new TreeMap<>();
                    for (VarPtm varPtm : dstMap.values()){
                        massVarPtmMap.put(varPtm.mass, varPtm);
                    }
                    pos_MassVarPtm_Map.put(absPos, massVarPtmMap);
                }
            }
        } // end for (int i = 0; i < partSeq.length(); ++i) {
    }

    public void findPosssible1Ptm(int scanNum, String partSeq, Map<Integer, TreeMap<Double, VarPtm>> pos_MassVarPtm_Map, int startRefPos, TreeSet<Peptide> modPepsSet, double deltaMass) {

        int partSeqLen = partSeq.length();
        int absPos;
        for (int relPos = 0; relPos < partSeqLen; ++relPos) {
            absPos = relPos + startRefPos;
            TreeMap<Double, VarPtm> massVarPtmMap = pos_MassVarPtm_Map.get(absPos);
            if (massVarPtmMap == null) continue;
            massVarPtmMap.subMap(deltaMass-ms2Tol, deltaMass+ms2Tol);
            for (double mass : massVarPtmMap.subMap(deltaMass-ms2Tol, deltaMass+ms2Tol).keySet()) {
                if (massTable.get(partSeq.charAt(relPos)) + mass < ms2Tol) continue;
                VarPtm varPtm = massVarPtmMap.get(mass);

                Peptide tmpPeptide = new Peptide(partSeq, false, massTool);
                PosMassMap posMassMap = new PosMassMap();
                posMassMap.put(relPos, mass);
                tmpPeptide.posVarPtmResMap.put(relPos, varPtm);
                tmpPeptide.setVarPTM(posMassMap);
                tmpPeptide.setScore(0);
                modPepsSet.add(tmpPeptide);
            }
        } // end for (int i = 0; i < partSeq.length(); ++i) {
    }

    public void prepareInfoNTerm(int scanNum, String partSeq,
                                 Map<Integer, TreeMap<Double, VarPtm>> absPos_MassVarPtm_Map,
                                 Map<Integer, Integer> yIdMinAbsPosMap,
                                 boolean couldBeProtN,
                                 int optEndPosP1,
                                 int endRefPos,
                                 String protSeq,
                                 Map<Integer, List<Byte>> absPos_ptmPositions_Map) {

        int partSeqLen = partSeq.length();
        char aa;
        int absPos;
        for (int relPos = 0; relPos < partSeqLen; ++relPos) {
            aa = partSeq.charAt(relPos);
            absPos = relPos + endRefPos - partSeqLen;
            List<Byte> positionsToTry = new ArrayList<>(3);
            absPos_ptmPositions_Map.put(absPos, positionsToTry);
            if (aaWithFixModSet.contains(aa) && absPos != 0) continue;  // if that pos has fix mod but is not N term, dont let it

            positionsToTry.add(ANYWHERE);
            if (yIdMinAbsPosMap.containsValue(absPos)) {
                positionsToTry.add(PEPN);
            }
            if (absPos == optEndPosP1) {
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
                    if ( ! allVarPtmMap.containsKey(position)) continue;
                    for (VarPtm varPtm : allVarPtmMap.get(position)) {
                        double mass = varPtm.mass;
                        if (massTable.get(partSeq.charAt(relPos)) + mass < ms2Tol) continue;
                        String varPtmStr = varPtm.getStr();

                        VarPtm oldVarPtm = dstMap.get(varPtmStr);
                        if (oldVarPtm != null) {
                            if (varPtm.priority > oldVarPtm.priority){
                                dstMap.put(varPtmStr, varPtm);
                            }
                        } else {
                            dstMap.put(varPtmStr, varPtm);
                        }

                    }
                }

                if (!dstMap.isEmpty()) {
                    TreeMap<Double, VarPtm> massVarPtmMap = new TreeMap<>();
                    for (VarPtm varPtm : dstMap.values()){
                        massVarPtmMap.put(varPtm.mass, varPtm);
                    }
                    absPos_MassVarPtm_Map.put(absPos, massVarPtmMap);
                }
            }
        } // end for (int i = 0; i < partSeq.length(); ++i) {

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
                        if (massTable.get(partSeq.charAt(i)) + varPtm.mass < ms2Tol) continue;//the mass of a modified amino acid cannot be 0 or negative

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
                        if (massTable.get(partSeq.charAt(i)) + varPtm.mass < ms2Tol) continue;//the mass of a modified amino acid cannot be 0 or negative

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
                            massesOnAa.add(varPtm.mass);
                        }
                    }
                } else {//aa at middle
                    for (VarPtm varPtm : srcSet) {
                        if (massTable.get(partSeq.charAt(i)) + varPtm.mass < ms2Tol) continue;//the mass of a modified amino acid cannot be 0 or negative

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
                        if (massTable.get(partSeq.charAt(i)) + varPtm.mass < ms2Tol) continue;//the mass of a modified amino acid cannot be 0 or negative

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
                        if (massTable.get(partSeq.charAt(i)) + varPtm.mass < ms2Tol) continue;//the mass of a modified amino acid cannot be 0 or negative

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
                        if (massTable.get(partSeq.charAt(i)) + varPtm.mass < ms2Tol) continue;//the mass of a modified amino acid cannot be 0 or negative

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
                        if (massTable.get(ptmFreePeptide.charAt(i)) + varPtm.mass < ms2Tol) continue;//the mass of a modified amino acid cannot be 0 or negative

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
                        if (massTable.get(ptmFreePeptide.charAt(i)) + varPtm.mass < ms2Tol) continue;//the mass of a modified amino acid cannot be 0 or negative

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
                        if (massTable.get(ptmFreePeptide.charAt(i)) + varPtm.mass < ms2Tol) continue;//the mass of a modified amino acid cannot be 0 or negative

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
