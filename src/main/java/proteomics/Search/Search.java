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

package proteomics.Search;

import proteomics.Index.BuildIndex;
import ProteomicsLibrary.MassTool;
import ProteomicsLibrary.Types.*;
import proteomics.Types.*;
import proteomics.PreSearch.Entry;
import java.util.*;

import static proteomics.PIPI.lszDebugScanNum;

public class Search {

    private static final int rankNum = 5;

    private List<Peptide> ptmOnlyResult = new LinkedList<>();
    private List<Peptide> ptmFreeResult = new LinkedList<>();
    public List<PepWithScore> candidatesList = new LinkedList<>();

    public Search(Entry entry, int scanNum, BuildIndex buildIndex, double precursorMass, SparseVector scanCode, MassTool massTool, double ms1Tolerance
            , double leftInverseMs1Tolerance, double rightInverseMs1Tolerance, int ms1ToleranceUnit, double minPtmMass, double maxPtmMass
            , int localMaxMs2Charge, Set<String> candiSet, String truth) {

        PriorityQueue<ResultEntry> ptmFreeQueue = new PriorityQueue<>(rankNum * 2);
        PriorityQueue<ResultEntry> ptmOnlyQueue = new PriorityQueue<>(rankNum * 2);
        double scanNormSquare = scanCode.norm2square();

        Map<String, PepInfo> peptide0Map = buildIndex.getPepInfoMap();

        boolean containsTruth = false;
        for (String sequence : candiSet) {
            PepInfo pepInfo = peptide0Map.get(sequence);
            double score = 0;
            double temp1 = Math.sqrt(pepInfo.code.norm2square() * scanNormSquare);
            if (temp1 > 1e-6) {
                score = pepInfo.code.dot(scanCode) / temp1;
            }
            double deltaMass = massTool.calResidueMass(sequence) + massTool.H2O - precursorMass; // caution: the order matters under ms1ToleranceUnit == 1 situation
//            if (Math.abs(deltaMass) > 250) continue;
//                    pepSet.add(sequence);


            if (pepInfo.isTarget) {
                if ((Math.abs(deltaMass) <= 0.01)) {

                    // PTM-free
                    if (ptmFreeQueue.size() < rankNum) {
                        ptmFreeQueue.add(new ResultEntry(score, sequence, false));
                    } else {
                        if (score > ptmFreeQueue.peek().score) {
                            ptmFreeQueue.poll();
                            ptmFreeQueue.add(new ResultEntry(score, sequence, false));
                        }
                    }
                }

//                        if ((deltaMass > rightTol) || (deltaMass < -1 * leftTol)) {
                if ((Math.abs(deltaMass) > 0.01)) {

                    // PTM-only
                    if (ptmOnlyQueue.size() < rankNum) {
                        ptmOnlyQueue.add(new ResultEntry(score, sequence, false));
                    } else {
                        if (score > ptmOnlyQueue.peek().score) {
                            ptmOnlyQueue.poll();
                            ptmOnlyQueue.add(new ResultEntry(score, sequence, false));
                        }
                    }
                }
            } else {
//                        if ((deltaMass <= rightTol) && (deltaMass >= -1 * leftTol)) {
                if ((Math.abs(deltaMass) <= 0.01)) {

                    // PTM-free
                    if (ptmFreeQueue.size() < rankNum) {
                        ptmFreeQueue.add(new ResultEntry(score, sequence, true));
                    } else {
                        if (score > ptmFreeQueue.peek().score) {
                            ptmFreeQueue.poll();
                            ptmFreeQueue.add(new ResultEntry(score, sequence, true));
                        }
                    }
                }

//                        if ((deltaMass > rightTol) || (deltaMass < -1 * leftTol)) {
                if ((Math.abs(deltaMass) > 0.01)) {

                    // PTM-only
                    if (ptmOnlyQueue.size() < rankNum) {
                        ptmOnlyQueue.add(new ResultEntry(score, sequence, true));
                    } else {
                        if (score > ptmOnlyQueue.peek().score) {
                            ptmOnlyQueue.poll();
                            ptmOnlyQueue.add(new ResultEntry(score, sequence, true));
                        }
                    }
                }
            }
        }
        if (!(ptmFreeQueue.isEmpty() && ptmOnlyQueue.isEmpty())) {
            mergeResult(ptmFreeQueue, ptmOnlyQueue, peptide0Map);
        }

        if (!(ptmFreeQueue.isEmpty() && ptmOnlyQueue.isEmpty())) {
            entry.ptmOnlyList = convertResult(ptmOnlyQueue, massTool, localMaxMs2Charge);
            entry.ptmFreeList = convertResult(ptmFreeQueue, massTool, localMaxMs2Charge);
        }
    }

    public Search(Entry entry, int scanNum, BuildIndex buildIndex, double precursorMass, SparseVector scanCode, MassTool massTool, double ms1Tolerance
            , double leftInverseMs1Tolerance, double rightInverseMs1Tolerance, int ms1ToleranceUnit, double minPtmMass, double maxPtmMass, int localMaxMs2Charge, List<ThreeExpAA> ncTags) {

        PriorityQueue<ResultEntry> ptmFreeQueue = new PriorityQueue<>(rankNum * 2);
        PriorityQueue<ResultEntry> ptmOnlyQueue = new PriorityQueue<>(rankNum * 2);
        double scanNormSquare = scanCode.norm2square();
        double leftTol = ms1Tolerance;
        double rightTol = ms1Tolerance;
        if (ms1ToleranceUnit == 1) {
            leftTol = precursorMass - (precursorMass * leftInverseMs1Tolerance);
            rightTol = (precursorMass * rightInverseMs1Tolerance) - precursorMass;
        }
        double leftMass = Math.max(precursorMass + minPtmMass - leftTol, buildIndex.getMinPeptideMass());
        double rightMass = Math.min(precursorMass + maxPtmMass + rightTol, buildIndex.getMaxPeptideMass());

        if (leftMass >= rightMass) {
            return;
        }

        Map<String, PepInfo> peptide0Map = buildIndex.getPepInfoMap();
        TreeMap<Double, Set<String>> massPeptideMap = buildIndex.getMassPeptideMap();

        NavigableMap<Double, Set<String>> subMassPeptideMap = massPeptideMap.subMap(leftMass, true, rightMass, true);

        if (!subMassPeptideMap.isEmpty()) {
            for (double mass : subMassPeptideMap.keySet()) {
                for (String sequence : massPeptideMap.get(mass)) {
                    if (lszDebugScanNum.contains(scanNum)) {
                        if(sequence.contentEquals("nGSKLDGTLKGRc") || sequence.contentEquals("nRKLGFADLLERc")) {
                            int a = 1;
                        }
                    }
                    PepInfo pepInfo = peptide0Map.get(sequence);
                    double score = 0;
                    double temp1 = Math.sqrt(pepInfo.code.norm2square() * scanNormSquare);
                    if (temp1 > 1e-6) {
                        score = pepInfo.code.dot(scanCode) / temp1;
                    }
//                    SparseBooleanVector aaa = new SparseBooleanVector(new HashSet<>(Arrays.asList(4305,2435,4339,2421,2635,2556,2541,2429,2430)));
//                    double sss = aaa.dot(scanCode)/temp1;

                    double deltaMass = mass - precursorMass; // caution: the order matters under ms1ToleranceUnit == 1 situation

                    if (pepInfo.isTarget) {
//                        if ((deltaMass <= rightTol) && (deltaMass >= -1 * leftTol)) {
                        if ((Math.abs(deltaMass) <= 0.01)) {

                            // PTM-free
                            if (ptmFreeQueue.size() < rankNum) {
                                ptmFreeQueue.add(new ResultEntry(score, sequence, false));
                            } else {
                                if (score > ptmFreeQueue.peek().score) {
                                    ptmFreeQueue.poll();
                                    ptmFreeQueue.add(new ResultEntry(score, sequence, false));
                                }
                            }
                        }

//                        if ((deltaMass > rightTol) || (deltaMass < -1 * leftTol)) {
                        if ((Math.abs(deltaMass) > 0.01)) {

                            // PTM-only
                            if (ptmOnlyQueue.size() < rankNum) {
                                ptmOnlyQueue.add(new ResultEntry(score, sequence, false));
                            } else {
                                if (score > ptmOnlyQueue.peek().score) {
                                    ptmOnlyQueue.poll();
                                    ptmOnlyQueue.add(new ResultEntry(score, sequence, false));
                                }
                            }
                        }
                    } else {
//                        if ((deltaMass <= rightTol) && (deltaMass >= -1 * leftTol)) {
                        if ((Math.abs(deltaMass) <= 0.01)) {

                            // PTM-free
                            if (ptmFreeQueue.size() < rankNum) {
                                ptmFreeQueue.add(new ResultEntry(score, sequence, true));
                            } else {
                                if (score > ptmFreeQueue.peek().score) {
                                    ptmFreeQueue.poll();
                                    ptmFreeQueue.add(new ResultEntry(score, sequence, true));
                                }
                            }
                        }

//                        if ((deltaMass > rightTol) || (deltaMass < -1 * leftTol)) {
                        if ((Math.abs(deltaMass) > 0.01)) {

                            // PTM-only
                            if (ptmOnlyQueue.size() < rankNum) {
                                ptmOnlyQueue.add(new ResultEntry(score, sequence, true));
                            } else {
                                if (score > ptmOnlyQueue.peek().score) {
                                    ptmOnlyQueue.poll();
                                    ptmOnlyQueue.add(new ResultEntry(score, sequence, true));
                                }
                            }
                        }
                    }
                }
            }
        }

        if (lszDebugScanNum.contains(scanNum)) {
            int a = 1;
        }
        if (!(ptmFreeQueue.isEmpty() && ptmOnlyQueue.isEmpty())) {
            mergeResult(ptmFreeQueue, ptmOnlyQueue, peptide0Map);
        }

        if (!(ptmFreeQueue.isEmpty() && ptmOnlyQueue.isEmpty())) {
            entry.ptmOnlyList = convertResult(ptmOnlyQueue, massTool, localMaxMs2Charge);
            entry.ptmFreeList = convertResult(ptmFreeQueue, massTool, localMaxMs2Charge);
        }
    }
    private void mergeResult(PriorityQueue<ResultEntry> ptmFreeQueue, PriorityQueue<ResultEntry> ptmOnlyQueue, Map<String, PepInfo> peptide0Map) {
        if (!ptmFreeQueue.isEmpty()) {
            for (ResultEntry temp : ptmFreeQueue){
                PepInfo pep0 = peptide0Map.get(temp.peptide);
                candidatesList.add(new PepWithScore(temp.peptide, temp.score, temp.isDecoy(), false, String.join("_", pep0.proteins)));
            }
        }

        if (!ptmOnlyQueue.isEmpty()) {
            for (ResultEntry temp : ptmOnlyQueue){
                PepInfo pep0 = peptide0Map.get(temp.peptide);
                candidatesList.add(new PepWithScore(temp.peptide, temp.score, temp.isDecoy(), true, String.join("_", pep0.proteins)));
            }
        }
    }

    private List<Peptide> convertResult(PriorityQueue<ResultEntry> inputQueue, MassTool massTool, int localMaxMs2Charge) {
        List<Peptide> peptideList = new LinkedList<>();
        int globalRank = inputQueue.size();
        while (!inputQueue.isEmpty()) {
            ResultEntry temp = inputQueue.poll();
            peptideList.add(new Peptide(temp.peptide, temp.isDecoy(), massTool, localMaxMs2Charge, temp.score, globalRank));
            --globalRank;
        }

        return peptideList;
    }
}
