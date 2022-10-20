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

package proteomics;
import ProteomicsLibrary.*;
import ProteomicsLibrary.Types.SparseBooleanVector;
import proteomics.Index.BuildIndex;

import java.util.*;
import java.util.concurrent.Callable;

public class GeneratePepDb implements Callable<GeneratePepDb.Entry> {
    private Map<String, String> parameterMap;
    private Set<String> reducedProtIdSet;

    private BuildIndex buildIndex;

    private String protId;
    public GeneratePepDb(Map<String, String> parameterMap, BuildIndex buildIndex, String protId) throws Exception {
        this.parameterMap = parameterMap;
        this.buildIndex = buildIndex;
        this.protId = protId;
    }

    @Override
    public Entry call() throws Exception {
        boolean addDecoy = parameterMap.get("add_decoy").contentEquals("1");
        int minPeptideLength = Math.max(5, Integer.valueOf(parameterMap.get("min_peptide_length")));
        int maxPeptideLength = Integer.valueOf(parameterMap.get("max_peptide_length"));

        Entry entry = new Entry();
        entry.protId = protId;
        String protSeq = buildIndex.protSeqMap.get(protId).replace('I', 'L');
        Set<String> peptideSet = buildIndex.massTool.buildPeptideSetPnP(protSeq);

        for (String peptide : peptideSet) {
            if (MassTool.containsNonAAAndNC(peptide)) {
                continue;
            }
            if ((peptide.length() - 2 <= maxPeptideLength) && (peptide.length() - 2 >= minPeptideLength)) { // caution: there are n and c in the sequence
                entry.targetPepMassMap.put(peptide, buildIndex.massTool.calResidueMass(peptide) + buildIndex.massTool.H2O);
                entry.targetPepCodeMap.put(peptide, buildIndex.inferSegment.generateSegmentBooleanVector(peptide));
            }
        }
        if (addDecoy) {
            //1
//                List<Character> proSeqAaList = proSeq.chars().mapToObj( c -> (char)c).collect(Collectors.toList());
//                Collections.shuffle(proSeqAaList);
//                String decoyProSeq = Joiner.on("").join(proSeqAaList);
            //2 swap every two, if same, change
//            String decoyProSeq = DbTool.shuffleSwapEveryTwo(protSeq, parameterMap.get("cleavage_site_1"), parameterMap.get("protection_site_1"), Integer.valueOf(parameterMap.get("is_from_C_term_1")) == 1).replace('I', 'L'); // FixMe: Only consider the first enzyme if the users specify two enzymes.
            //3
//            String decoyProSeq = DbTool.shuffleProtBetweenKR(protSeq, parameterMap.get("cleavage_site_1"), parameterMap.get("protection_site_1"), Integer.valueOf(parameterMap.get("is_from_C_term_1")) == 1).replace('I', 'L'); // FixMe: Only consider the first enzyme if the users specify two enzymes.
            //4
            String decoyProSeq = DbTool.shuffleProtKeepKR(protSeq, parameterMap.get("cleavage_site_1"), parameterMap.get("protection_site_1"), Integer.valueOf(parameterMap.get("is_from_C_term_1")) == 1).replace('I', 'L'); // FixMe: Only consider the first enzyme if the users specify two enzymes.

            entry.decoyProtSeq = decoyProSeq;
            peptideSet = buildIndex. massTool.buildPeptideSetPnP(decoyProSeq);
            for (String peptide : peptideSet) {
                if (MassTool.containsNonAAAndNC(peptide)) {
                    continue;
                }
                if ((peptide.length() - 2 <= maxPeptideLength) && (peptide.length() - 2 >= minPeptideLength)) { // caution: there are n and c in the sequence
                    entry.decoyPepMassMap.put(peptide, buildIndex.massTool.calResidueMass(peptide) + buildIndex.massTool.H2O);
                    entry.decoyPepCodeMap.put(peptide, buildIndex.inferSegment.generateSegmentBooleanVector(peptide));
                }
            }
        }
        return entry;
    }

    public class Entry {
        String protId = "";
        String decoyProtSeq = "";
        Map<String, SparseBooleanVector> targetPepCodeMap = new HashMap<>();
        Map<String, SparseBooleanVector> decoyPepCodeMap = new HashMap<>();
        Map<String, Double> targetPepMassMap = new HashMap<>();
        Map<String, Double> decoyPepMassMap = new HashMap<>();
        Entry() {
        }
    }
}
