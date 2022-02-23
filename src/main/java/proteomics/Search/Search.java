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

import java.util.*;

public class Search {

    private static final int rankNum = 1;

    private List<Peptide> ptmOnlyResult = new LinkedList<>();
    private List<Peptide> ptmFreeResult = new LinkedList<>();
    public List<PepWithScore> candidatesList = new LinkedList<>();

    public Search(int scanNum, String pepTruth, BuildIndex buildIndex, double precursorMass, SparseVector scanCode, MassTool massTool, double ms1Tolerance, double leftInverseMs1Tolerance, double rightInverseMs1Tolerance, int ms1ToleranceUnit, double minPtmMass, double maxPtmMass, int localMaxMs2Charge) {
        PriorityQueue<ResultEntry> ptmOnlyQueue = new PriorityQueue<>(rankNum);
        double scanNormSquare = scanCode.norm2square();

        Map<String, Peptide0> peptide0Map = buildIndex.getPeptide0Map();

        String sequence = "n"+pepTruth+"c";
        double mass = massTool.calResidueMass(sequence) + massTool.H2O;
        Peptide0 peptide0 = peptide0Map.get(sequence);
        if (peptide0 == null){
            int a = 1;
            System.out.println("contains "+peptide0Map.containsKey(sequence));
        }
        double score = 0;
        double temp1 = Math.sqrt(peptide0.code.norm2square() * scanNormSquare);
        if (temp1 > 1e-6) {
            score = peptide0.code.dot(scanCode) / temp1;
        }
        double deltaMass = mass - precursorMass; // caution: the order matters under ms1ToleranceUnit == 1 situation


        ptmOnlyQueue.add(new ResultEntry(score, sequence, false));

        mergeResult(ptmOnlyQueue, peptide0Map);
        ptmOnlyResult = convertResult(ptmOnlyQueue, massTool, localMaxMs2Charge);
    }

    private void mergeResult(PriorityQueue<ResultEntry> ptmOnlyQueue, Map<String, Peptide0> peptide0Map) {

        for (ResultEntry temp : ptmOnlyQueue){
            Peptide0 pep0 = peptide0Map.get(temp.peptide);
            candidatesList.add(new PepWithScore(temp.peptide, temp.score, temp.isDecoy(), true, String.join("_", pep0.proteins)));
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

    public List<Peptide> getPTMOnlyResult() {
        return ptmOnlyResult;
    }

    public List<Peptide> getPTMFreeResult() {
        return ptmFreeResult;
    }
}
