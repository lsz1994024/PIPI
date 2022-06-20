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

package proteomics.Index;
//import org.apache.commons.lang.StringUtils;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigInteger;
import java.util.*;

import proteomics.Hash.HashFunc;
//import proteomics.Hash.HashFunc;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import proteomics.PTM.InferPTM;
import proteomics.Segment.InferSegment;
import ProteomicsLibrary.*;
import ProteomicsLibrary.Types.*;
import proteomics.Types.Peptide0;

//import static proteomics.Hash.HashFunc.simhash64;

public class BuildIndex {
    private final int hashbits = 16;
    private final MassTool massTool;
    private Map<Character, Double> fixModMap = new HashMap<>(25, 1);
    private double minPeptideMass = 9999;
    private double maxPeptideMass = 0;
    private InferSegment inferSegment;
    private TreeMap<Double, Set<String>> massPeptideMap = new TreeMap<>();
    private Map<String, Peptide0> peptide0Map;
    public Map<Integer, LinkedList<String>> testMap = new HashMap<>();

    private final String labelling;
    private final DbTool dbTool; // this one doesn't contain contaminant proteins.
    private InferPTM inferPTM;
    public Map<Integer, Integer> truthHashMap = new HashMap<>();
//    public Map<Integer, Long> truthHashMap64 = new HashMap<>();

    public BuildIndex( Map<String, String> parameterMap, String labelling, boolean needCoding, boolean addDecoy, boolean addContaminant) throws Exception {
        // initialize parameters
        int minPeptideLength = Math.max(5, Integer.valueOf(parameterMap.get("min_peptide_length")));
        int maxPeptideLength = Integer.valueOf(parameterMap.get("max_peptide_length"));
        String dbPath = parameterMap.get("db");
        int missedCleavage = Integer.valueOf(parameterMap.get("missed_cleavage"));
        double ms2Tolerance = Double.valueOf(parameterMap.get("ms2_tolerance"));
        double oneMinusBinOffset = 1 - Double.valueOf(parameterMap.get("mz_bin_offset"));
        this.labelling = labelling;

        // Read fix modification
        fixModMap.put('G', Double.valueOf(parameterMap.get("G")));
        fixModMap.put('A', Double.valueOf(parameterMap.get("A")));
        fixModMap.put('S', Double.valueOf(parameterMap.get("S")));
        fixModMap.put('P', Double.valueOf(parameterMap.get("P")));
        fixModMap.put('V', Double.valueOf(parameterMap.get("V")));
        fixModMap.put('T', Double.valueOf(parameterMap.get("T")));
        fixModMap.put('C', Double.valueOf(parameterMap.get("C")));
        fixModMap.put('I', Double.valueOf(parameterMap.get("I")));
        fixModMap.put('L', Double.valueOf(parameterMap.get("L")));
        fixModMap.put('N', Double.valueOf(parameterMap.get("N")));
        fixModMap.put('D', Double.valueOf(parameterMap.get("D")));
        fixModMap.put('Q', Double.valueOf(parameterMap.get("Q")));
        fixModMap.put('K', Double.valueOf(parameterMap.get("K")));
        fixModMap.put('E', Double.valueOf(parameterMap.get("E")));
        fixModMap.put('M', Double.valueOf(parameterMap.get("M")));
        fixModMap.put('H', Double.valueOf(parameterMap.get("H")));
        fixModMap.put('F', Double.valueOf(parameterMap.get("F")));
        fixModMap.put('R', Double.valueOf(parameterMap.get("R")));
        fixModMap.put('Y', Double.valueOf(parameterMap.get("Y")));
        fixModMap.put('W', Double.valueOf(parameterMap.get("W")));
        fixModMap.put('U', Double.valueOf(parameterMap.get("U")));
        fixModMap.put('O', Double.valueOf(parameterMap.get("O")));
        fixModMap.put('n', Double.valueOf(parameterMap.get("n")));
        fixModMap.put('c', Double.valueOf(parameterMap.get("c")));

        // read protein database
        dbTool = new DbTool(dbPath, parameterMap.get("database_type"));
        Map<String, String> proteinPeptideMap;
        DbTool contaminantsDb = null;
        if (addContaminant) {
            contaminantsDb = new DbTool(null, "contaminants");
            proteinPeptideMap = contaminantsDb.getProteinSequenceMap();
            proteinPeptideMap.putAll(dbTool.getProteinSequenceMap()); // using the target sequence to replace contaminant sequence if there is conflict.
        } else {
            proteinPeptideMap = dbTool.getProteinSequenceMap();
        }

        // define a new MassTool object
        massTool = new MassTool(missedCleavage, fixModMap, parameterMap.get("cleavage_site_1").trim(), parameterMap.get("protection_site_1").trim(), parameterMap.get("is_from_C_term_1").trim().contentEquals("1"), parameterMap.getOrDefault("cleavage_site_2", null), parameterMap.getOrDefault("protection_site_2", null), parameterMap.containsKey("is_from_C_term_2") ? parameterMap.get("is_from_C_term_2").trim().contentEquals("1") : null, ms2Tolerance, oneMinusBinOffset, labelling);

        inferPTM = new InferPTM(massTool, fixModMap, parameterMap);

        // build database
        inferSegment = new InferSegment(massTool, parameterMap, fixModMap);

        Set<String> forCheckDuplicate = new HashSet<>(500000);
        Multimap<String, String> peptideProteinMap = HashMultimap.create();
        Map<String, Double> peptideMassMap = new HashMap<>(500000);
        Map<String, String> targetDecoyProteinSequenceMap = new HashMap<>();
        for (String proId : proteinPeptideMap.keySet()) {
            String proSeq = proteinPeptideMap.get(proId);
            if (proId.contentEquals("sp|Q06945|SOX4_HUMAN")){
                int a = 1;
            }
            Set<String> peptideSet = massTool.buildPeptideSetPnP(proSeq);
            for (String peptide : peptideSet) {
                if (MassTool.containsNonAAAndNC(peptide)) {
                    continue;
                }
                if (peptide.contentEquals("nVGGSGGGGHGGGGGGGSSNAGGGGGGASGGGANSKc")) {
                    int a = 1;
                }

                if ((peptide.length() - 2 <= maxPeptideLength) && (peptide.length() - 2 >= minPeptideLength)) { // caution: there are n and c in the sequence
                    if (!forCheckDuplicate.contains(peptide.replace('L', 'I'))) { // don't record duplicate peptide sequences
                        // Add the sequence to the check set for duplicate check
                        forCheckDuplicate.add(peptide.replace('L', 'I'));

                        double mass = massTool.calResidueMass(peptide) + massTool.H2O;
                        // recode min and max peptide mass
                        if (mass < minPeptideMass) {
                            minPeptideMass = mass;
                        }
                        if (mass > maxPeptideMass) {
                            maxPeptideMass = mass;
                        }

                        peptideMassMap.put(peptide, mass);
                        peptideProteinMap.put(peptide, proId);
                    } else if (peptideProteinMap.containsKey(peptide)) {
                        // Considering the case that the sequence has multiple proteins. In the above if block, such a protein ID wasn't recorded. If there are decoy IDs, replace it with the current target ID since the target ID has a higher priority.
                        Set<String> proteinSet = new HashSet<>(peptideProteinMap.get(peptide));
                        peptideProteinMap.get(peptide).clear();
                        for (String protein : proteinSet) {
                            if (!protein.startsWith("DECOY_")) {
                                peptideProteinMap.put(peptide, protein);
                            }
                        }
                        peptideProteinMap.put(peptide, proId);
                    }
                }
            }

            targetDecoyProteinSequenceMap.put(proId, proSeq);

            if (addDecoy) {
                // decoy sequence
                String decoyProSeq = DbTool.shuffleSeq(proSeq, parameterMap.get("cleavage_site_1"), parameterMap.get("protection_site_1"), Integer.valueOf(parameterMap.get("is_from_C_term_1")) == 1); // FixMe: Only consider the first enzyme if the users specify two enzymes.
                peptideSet = massTool.buildPeptideSetPnP(decoyProSeq);

                for (String peptide : peptideSet) {
                    if (MassTool.containsNonAAAndNC(peptide)) {
                        continue;
                    }

                    if ((peptide.length() - 2 <= maxPeptideLength) && (peptide.length() - 2 >= minPeptideLength)) { // caution: there are n and c in the sequence
                        if (!forCheckDuplicate.contains(peptide.replace('L', 'I'))) { // don't record duplicate peptide sequences
                            // Add the sequence to the check set for duplicate check
                            forCheckDuplicate.add(peptide.replace('L', 'I'));

                            double mass = massTool.calResidueMass(peptide) + massTool.H2O;
                            // recode min and max peptide mass
                            if (mass < minPeptideMass) {
                                minPeptideMass = mass;
                            }
                            if (mass > maxPeptideMass) {
                                maxPeptideMass = mass;
                            }

                            peptideMassMap.put(peptide, mass);
                            peptideProteinMap.put(peptide, "DECOY_" + proId);
                        }
                    }
                }
                targetDecoyProteinSequenceMap.put("DECOY_" + proId, decoyProSeq);
            }
        }

        if (addDecoy) {
            // writer concatenated fasta
            Map<String, String> proteinAnnotationMap;
            if (addContaminant) {
                proteinAnnotationMap = contaminantsDb.getProteinAnnotateMap();
                proteinAnnotationMap.putAll(dbTool.getProteinAnnotateMap()); // using the target annotation to replace contaminant sequence if there is conflict.
            } else {
                proteinAnnotationMap = dbTool.getProteinAnnotateMap();
            }

            BufferedWriter writer = new BufferedWriter(new FileWriter(dbPath + ".TD.fasta"));
            for (String proId : targetDecoyProteinSequenceMap.keySet()) {
                writer.write(String.format(Locale.US, ">%s %s\n", proId, proteinAnnotationMap.getOrDefault(proId, "")));
                writer.write(targetDecoyProteinSequenceMap.get(proId) + "\n");
            }
            writer.close();
        }

        Map<String, Integer> pepCountMap = new HashMap<>(); //This is global, do it first
        Map<String, Map<String, Double>> tfMap = new HashMap<>();
//        Map<String, Map<String, Double>> tfidfMap = new HashMap<>();

        Map<String, Peptide0> tempMap = new HashMap<>();
        for (String peptide : peptideMassMap.keySet()) {
            Map<String, Double> tfForOne = new HashMap<>();
            SparseBooleanVector code = null;
            if (needCoding) {
                code = inferSegment.generateSegmentBooleanVector(DbTool.getSequenceOnly(peptide), pepCountMap, tfForOne);
            }
            tfMap.put(peptide, tfForOne);
            Character[] leftRightFlank = DbTool.getLeftRightFlank(peptide, peptideProteinMap, targetDecoyProteinSequenceMap, parameterMap.get("cleavage_site_1"), parameterMap.get("protection_site_1"), parameterMap.get("is_from_C_term_1").contentEquals("1")); // FixMe: Only consider the first enzyme if the users specify two enzymes.
            if (leftRightFlank == null) {
                leftRightFlank = DbTool.getLeftRightFlank(peptide, peptideProteinMap, targetDecoyProteinSequenceMap, parameterMap.get("cleavage_site_2"), parameterMap.get("protection_site_2"), parameterMap.get("is_from_C_term_2").contentEquals("1")); // FixMe: Only consider the first enzyme if the users specify two enzymes.
            }
            if (leftRightFlank != null) {
                tempMap.put(peptide, new Peptide0(code, isTarget(peptideProteinMap.get(peptide)), peptideProteinMap.get(peptide).toArray(new String[0]), leftRightFlank[0], leftRightFlank[1]));

                if (massPeptideMap.containsKey(peptideMassMap.get(peptide))) {
                    massPeptideMap.get(peptideMassMap.get(peptide)).add(peptide);
                } else {
                    Set<String> tempSet = new HashSet<>();
                    tempSet.add(peptide);
                    massPeptideMap.put(peptideMassMap.get(peptide), tempSet);
                }
            }
        }


//        int totalPepNum = tempMap.size();
//        Map<BigInteger, Integer> testMap = new HashMap<>();
//
        peptide0Map = new HashMap<>(tempMap); // Since this map won't be changed any more, using this step to create a HashMap with the capacity exactly equals the actual size.
        int numTargetPep = 0, numDecoyPep = 0;
        int aa = 1;
        int totalPepNum = peptide0Map.size();
        HashFunc hashFunc = new HashFunc();
//        for (String pep : peptide0Map.keySet()){
//            aa++;
//            long pepHash = hashFunc.simhash64(tfMap.get(pep), pepCountMap,totalPepNum);
//            peptide0Map.get(pep).longHash = pepHash;
////            BigInteger bigInt = simHash(tfMap.get(pep), pepCountMap);
////            peptide0Map.get(pep).binHash = bigInt;
//            if (testMap.containsKey(pepHash)) {
////                testMap.put(bigInt, testMap.get(bigInt) + 1);
//                testMap.get(pepHash).add(pep);
//
//            } else {
////                Set<String> testSet = new HashSet<>();
////                testSet.add(pep);
//                LinkedList<String> testList = new LinkedList<>();
//                testList.add(pep);
//                testMap.put(pepHash, testList);
//            }
//        }
//         for(long pepHash : testMap.keySet())  {
////             System.out.println(pepHash + ","+testMap.get(pepHash).size());
//             if (testMap.get(pepHash).size() > 2) {
//                 int a = 1;
//             }
//         }

//        for(long pepHash1 : testMap.keySet())  {
//            for (long pepHash2 : testMap.keySet()) {
//                if (pepHash1== pepHash2 ) continue;
//
//                int hamDist = hashFunc.hammingDistance(pepHash1,pepHash2);
//                if (hamDist <=3 ) {
//                    LinkedList<String> pepSet2 = testMap.get(pepHash2);
//                    LinkedList<String> pepSet1 = testMap.get(pepHash1);
//                    int a = 1;
//                }
//            }
//        }
        int total  = 0;
        for (LinkedList<String> testList  : testMap.values()) total+= testList.size();
//        int total  = 0;
//        for (int num : testMap.values()) total+= num;
        int a = 1;
    }
    public BuildIndex(Map<Integer, String> pepTruth, Map<String, String> parameterMap, String labelling, boolean needCoding, boolean addDecoy, boolean addContaminant) throws Exception {
        // initialize parameters
        int minPeptideLength = Math.max(5, Integer.valueOf(parameterMap.get("min_peptide_length")));
        int maxPeptideLength = Integer.valueOf(parameterMap.get("max_peptide_length"));
        String dbPath = parameterMap.get("db");
        int missedCleavage = Integer.valueOf(parameterMap.get("missed_cleavage"));
        double ms2Tolerance = Double.valueOf(parameterMap.get("ms2_tolerance"));
        double oneMinusBinOffset = 1 - Double.valueOf(parameterMap.get("mz_bin_offset"));
        this.labelling = labelling;

        // Read fix modification
        fixModMap.put('G', Double.valueOf(parameterMap.get("G")));
        fixModMap.put('A', Double.valueOf(parameterMap.get("A")));
        fixModMap.put('S', Double.valueOf(parameterMap.get("S")));
        fixModMap.put('P', Double.valueOf(parameterMap.get("P")));
        fixModMap.put('V', Double.valueOf(parameterMap.get("V")));
        fixModMap.put('T', Double.valueOf(parameterMap.get("T")));
        fixModMap.put('C', Double.valueOf(parameterMap.get("C")));
        fixModMap.put('I', Double.valueOf(parameterMap.get("I")));
        fixModMap.put('L', Double.valueOf(parameterMap.get("L")));
        fixModMap.put('N', Double.valueOf(parameterMap.get("N")));
        fixModMap.put('D', Double.valueOf(parameterMap.get("D")));
        fixModMap.put('Q', Double.valueOf(parameterMap.get("Q")));
        fixModMap.put('K', Double.valueOf(parameterMap.get("K")));
        fixModMap.put('E', Double.valueOf(parameterMap.get("E")));
        fixModMap.put('M', Double.valueOf(parameterMap.get("M")));
        fixModMap.put('H', Double.valueOf(parameterMap.get("H")));
        fixModMap.put('F', Double.valueOf(parameterMap.get("F")));
        fixModMap.put('R', Double.valueOf(parameterMap.get("R")));
        fixModMap.put('Y', Double.valueOf(parameterMap.get("Y")));
        fixModMap.put('W', Double.valueOf(parameterMap.get("W")));
        fixModMap.put('U', Double.valueOf(parameterMap.get("U")));
        fixModMap.put('O', Double.valueOf(parameterMap.get("O")));
        fixModMap.put('n', Double.valueOf(parameterMap.get("n")));
        fixModMap.put('c', Double.valueOf(parameterMap.get("c")));

        // read protein database
        dbTool = new DbTool(dbPath, parameterMap.get("database_type"));
        Map<String, String> proteinPeptideMap;
        DbTool contaminantsDb = null;
        if (addContaminant) {
            contaminantsDb = new DbTool(null, "contaminants");
            proteinPeptideMap = contaminantsDb.getProteinSequenceMap();
            proteinPeptideMap.putAll(dbTool.getProteinSequenceMap()); // using the target sequence to replace contaminant sequence if there is conflict.
        } else {
            proteinPeptideMap = dbTool.getProteinSequenceMap();
        }

        // define a new MassTool object
        massTool = new MassTool(missedCleavage, fixModMap, parameterMap.get("cleavage_site_1").trim(), parameterMap.get("protection_site_1").trim(), parameterMap.get("is_from_C_term_1").trim().contentEquals("1"), parameterMap.getOrDefault("cleavage_site_2", null), parameterMap.getOrDefault("protection_site_2", null), parameterMap.containsKey("is_from_C_term_2") ? parameterMap.get("is_from_C_term_2").trim().contentEquals("1") : null, ms2Tolerance, oneMinusBinOffset, labelling);

        inferPTM = new InferPTM(massTool, fixModMap, parameterMap);

        // build database
        inferSegment = new InferSegment(massTool, parameterMap, fixModMap);

        Set<String> forCheckDuplicate = new HashSet<>(500000);
        Multimap<String, String> peptideProteinMap = HashMultimap.create();
        Map<String, Double> peptideMassMap = new HashMap<>(500000);
        Map<String, String> targetDecoyProteinSequenceMap = new HashMap<>();
        for (String proId : proteinPeptideMap.keySet()) {
            String proSeq = proteinPeptideMap.get(proId);
            if (proId.contentEquals("sp|Q06945|SOX4_HUMAN")){
                int a = 1;
            }
            Set<String> peptideSet = massTool.buildPeptideSetPnP(proSeq);
            for (String peptide : peptideSet) {
                if (MassTool.containsNonAAAndNC(peptide)) {
                    continue;
                }
                if (peptide.contentEquals("nVGGSGGGGHGGGGGGGSSNAGGGGGGASGGGANSKc")) {
                    int a = 1;
                }

                if ((peptide.length() - 2 <= maxPeptideLength) && (peptide.length() - 2 >= minPeptideLength)) { // caution: there are n and c in the sequence
                    if (!forCheckDuplicate.contains(peptide.replace('L', 'I'))) { // don't record duplicate peptide sequences
                        // Add the sequence to the check set for duplicate check
                        forCheckDuplicate.add(peptide.replace('L', 'I'));

                        double mass = massTool.calResidueMass(peptide) + massTool.H2O;
                        // recode min and max peptide mass
                        if (mass < minPeptideMass) {
                            minPeptideMass = mass;
                        }
                        if (mass > maxPeptideMass) {
                            maxPeptideMass = mass;
                        }

                        peptideMassMap.put(peptide, mass);
                        peptideProteinMap.put(peptide, proId);
                    } else if (peptideProteinMap.containsKey(peptide)) {
                        // Considering the case that the sequence has multiple proteins. In the above if block, such a protein ID wasn't recorded. If there are decoy IDs, replace it with the current target ID since the target ID has a higher priority.
                        Set<String> proteinSet = new HashSet<>(peptideProteinMap.get(peptide));
                        peptideProteinMap.get(peptide).clear();
                        for (String protein : proteinSet) {
                            if (!protein.startsWith("DECOY_")) {
                                peptideProteinMap.put(peptide, protein);
                            }
                        }
                        peptideProteinMap.put(peptide, proId);
                    }
                }
            }

            targetDecoyProteinSequenceMap.put(proId, proSeq);

            if (addDecoy) {
                // decoy sequence
                String decoyProSeq = DbTool.shuffleSeq(proSeq, parameterMap.get("cleavage_site_1"), parameterMap.get("protection_site_1"), Integer.valueOf(parameterMap.get("is_from_C_term_1")) == 1); // FixMe: Only consider the first enzyme if the users specify two enzymes.
                peptideSet = massTool.buildPeptideSetPnP(decoyProSeq);

                for (String peptide : peptideSet) {
                    if (MassTool.containsNonAAAndNC(peptide)) {
                        continue;
                    }

                    if ((peptide.length() - 2 <= maxPeptideLength) && (peptide.length() - 2 >= minPeptideLength)) { // caution: there are n and c in the sequence
                        if (!forCheckDuplicate.contains(peptide.replace('L', 'I'))) { // don't record duplicate peptide sequences
                            // Add the sequence to the check set for duplicate check
                            forCheckDuplicate.add(peptide.replace('L', 'I'));

                            double mass = massTool.calResidueMass(peptide) + massTool.H2O;
                            // recode min and max peptide mass
                            if (mass < minPeptideMass) {
                                minPeptideMass = mass;
                            }
                            if (mass > maxPeptideMass) {
                                maxPeptideMass = mass;
                            }

                            peptideMassMap.put(peptide, mass);
                            peptideProteinMap.put(peptide, "DECOY_" + proId);
                        }
                    }
                }
                targetDecoyProteinSequenceMap.put("DECOY_" + proId, decoyProSeq);
            }
        }

        if (addDecoy) {
            // writer concatenated fasta
            Map<String, String> proteinAnnotationMap;
            if (addContaminant) {
                proteinAnnotationMap = contaminantsDb.getProteinAnnotateMap();
                proteinAnnotationMap.putAll(dbTool.getProteinAnnotateMap()); // using the target annotation to replace contaminant sequence if there is conflict.
            } else {
                proteinAnnotationMap = dbTool.getProteinAnnotateMap();
            }

            BufferedWriter writer = new BufferedWriter(new FileWriter(dbPath + ".TD.fasta"));
            for (String proId : targetDecoyProteinSequenceMap.keySet()) {
                writer.write(String.format(Locale.US, ">%s %s\n", proId, proteinAnnotationMap.getOrDefault(proId, "")));
                writer.write(targetDecoyProteinSequenceMap.get(proId) + "\n");
            }
            writer.close();
        }

        Map<String, Integer> pepCountMap = new HashMap<>(); //This is global, do it first
        Map<String, Map<String, Double>> tfMap = new HashMap<>();
//        Map<String, Map<String, Double>> tfidfMap = new HashMap<>();

        Map<String, Peptide0> tempMap = new HashMap<>();
        for (String peptide : peptideMassMap.keySet()) {
            Map<String, Double> tfForOne = new HashMap<>();
            SparseBooleanVector code = null;
            if (needCoding) {
                code = inferSegment.generateSegmentBooleanVector(DbTool.getSequenceOnly(peptide), pepCountMap, tfForOne);
            }
            tfMap.put(peptide, tfForOne);
            Character[] leftRightFlank = DbTool.getLeftRightFlank(peptide, peptideProteinMap, targetDecoyProteinSequenceMap, parameterMap.get("cleavage_site_1"), parameterMap.get("protection_site_1"), parameterMap.get("is_from_C_term_1").contentEquals("1")); // FixMe: Only consider the first enzyme if the users specify two enzymes.
            if (leftRightFlank == null) {
                leftRightFlank = DbTool.getLeftRightFlank(peptide, peptideProteinMap, targetDecoyProteinSequenceMap, parameterMap.get("cleavage_site_2"), parameterMap.get("protection_site_2"), parameterMap.get("is_from_C_term_2").contentEquals("1")); // FixMe: Only consider the first enzyme if the users specify two enzymes.
            }
            if (leftRightFlank != null) {
                tempMap.put(peptide, new Peptide0(code, isTarget(peptideProteinMap.get(peptide)), peptideProteinMap.get(peptide).toArray(new String[0]), leftRightFlank[0], leftRightFlank[1]));

                if (massPeptideMap.containsKey(peptideMassMap.get(peptide))) {
                    massPeptideMap.get(peptideMassMap.get(peptide)).add(peptide);
                } else {
                    Set<String> tempSet = new HashSet<>();
                    tempSet.add(peptide);
                    massPeptideMap.put(peptideMassMap.get(peptide), tempSet);
                }
            }
        }

        int[] bitCount = new int[32];
        HashFunc temp = new HashFunc();
//        Set<Integer> recordedHash = new HashSet<>();
        for (String tag : pepCountMap.keySet()) {

            int v = temp.hash32(tag);
//            recordedHash.add(v);
            for (int i = 32; i >= 1; --i) {
                if (((v >> (32 - i)) & 1) == 1) {

                    bitCount[i - 1] += 1;
                }
//                else
//                    bitCount[i - 1] -= 1;
            }
        }
        for (int i : bitCount) {
//            System.out.println(i);
        }


//        int totalPepNum = tempMap.size();
//        Map<BigInteger, Integer> testMap = new HashMap<>();
//
        peptide0Map = new HashMap<>(tempMap); // Since this map won't be changed any more, using this step to create a HashMap with the capacity exactly equals the actual size.
        int numTargetPep = 0, numDecoyPep = 0;
        int aa = 1;
        int totalPepNum = peptide0Map.size();
        HashFunc hashFunc = new HashFunc();
        for (String pep : peptide0Map.keySet()){
            if (pep.contentEquals("nHAVSEGTKc")) {
                int a = 1;
            }
            aa++;
            int pepHash = hashFunc.simhash32(tfMap.get(pep), pepCountMap,totalPepNum);
            peptide0Map.get(pep).longHash = pepHash;
//            BigInteger bigInt = simHash(tfMap.get(pep), pepCountMap);
//            peptide0Map.get(pep).binHash = bigInt;
            if (testMap.containsKey(pepHash)) {
//                testMap.put(bigInt, testMap.get(bigInt) + 1);
                testMap.get(pepHash).add(pep);

            } else {
//                Set<String> testSet = new HashSet<>();
//                testSet.add(pep);
                LinkedList<String> testList = new LinkedList<>();
                testList.add(pep);
                testMap.put(pepHash, testList);
            }
        }
//         for(long pepHash : testMap.keySet())  {
////             System.out.println(pepHash + ","+testMap.get(pepHash).size());
//             if (testMap.get(pepHash).size() > 2) {
//                 int a = 1;
//             }
//         }

//        for(long pepHash1 : testMap.keySet())  {
//            LinkedList<String> pepSet1 = testMap.get(pepHash1);
//            LinkedList<String> neighbors = new LinkedList<>();
//            neighbors.addAll(pepSet1);
//            for (long pepHash2 : testMap.keySet()) {
//                if (pepHash1== pepHash2 ) continue;
//
//                int hamDist = hashFunc.hammingDistance(pepHash1,pepHash2);
//                if (hamDist < 8 ) {
//                    LinkedList<String> pepSet2 = testMap.get(pepHash2);
//                    neighbors.addAll(pepSet2);
//                    int a = 1;
//                }
//            }
//            if (neighbors.size() > 4 ){
//                int a = 1;
//            }
//
//            if (neighbors.size() > 7 ){
//                int a = 1;
//            }
//
//            int a = 1;
//        }
//        int total  = 0;
//        for (LinkedList<String> testList  : testMap.values()) total+= testList.size();
//        int total  = 0;
//        for (int num : testMap.values()) total+= num;
        int a = 1;
        for (int scan : pepTruth.keySet()) {
            if (scan == 1886){
                int aaa = 1;
            }
            String truthStr = pepTruth.get(scan);
            if (peptide0Map.containsKey(truthStr)){
                truthHashMap.put(scan, peptide0Map.get(truthStr).longHash);
            } else {
//                System.out.println("not containning truth of , "+ scan);
            }
        }
//        for (Peptide0 pep : peptide0Map.values()){
//
//            if (pep.isTarget) {
//                numTargetPep++;
//            }else {
//                numDecoyPep++;
//            }
//        }
    }


    private BigInteger simHash(Map<String, Double> tfMap, Map<String, Integer> pepCountMap) {
        int[] v = new int[hashbits];
        int totalPepNum = pepCountMap.size();
        for (String tag : tfMap.keySet()) {
            double tf = tfMap.get(tag);
            // 2、将每一个分词hash为一组固定长度的数列.比如 64bit 的一个整数.
            Double weight = 1 * tf * Math.log10(totalPepNum/pepCountMap.get(tag)); // 添加权重，权重应改为出现次数，而不是根据词性来指定。
            BigInteger t = this.hash(tag);
            for (int i = 0; i < hashbits; i++) {
                BigInteger bitmask = new BigInteger("1").shiftLeft(i);
                // 3、建立一个长度为64的整数数组(假设要生成64位的数字指纹,也可以是其它数字),
                // 对每一个分词hash后的数列进行判断,如果是1000...1,那么数组的第一位和末尾一位加1,
                // 中间的62位减一,也就是说,逢1加1,逢0减1.一直到把所有的分词hash数列全部判断完毕.
//				if(weight==null) continue;
                // if (wordCount.containsKey(word)) {
                // weight = wordCount.get(word);
                // }
                if (t.and(bitmask).signum() != 0) {
                    // 这里是计算整个文档的所有特征的向量和
                    v[i] += weight;
                } else {
                    v[i] -= weight;
                }
            }
        }

        //binString to BigInteger

        BigInteger fingerprint = new BigInteger("0");
        for (int i = 0; i < this.hashbits; i++) {
            if (v[i] >= 0) {
                fingerprint = fingerprint.add(new BigInteger("1").shiftLeft(i));
            }
        }
        return fingerprint;
    }

    private BigInteger hash(String source) {
        if (source == null || source.length() == 0) {
            return new BigInteger("0");
        } else {
            /**
             * 当sourece 的长度过短，会导致hash算法失效，因此需要对过短的词补偿
             */
            while (source.length() < 3) {
                source = source + source.charAt(0);
            }
            char[] sourceArray = source.toCharArray();
            BigInteger x = BigInteger.valueOf(((long) sourceArray[0]) << 7);
            BigInteger m = new BigInteger("1000003");
            BigInteger mask = new BigInteger("2").pow(this.hashbits).subtract(new BigInteger("1"));
            for (char item : sourceArray) {
                BigInteger temp = BigInteger.valueOf((long) item);
                x = x.multiply(m).xor(temp).and(mask);
            }
            x = x.xor(new BigInteger(String.valueOf(source.length())));
            if (x.equals(new BigInteger("-1"))) {
                x = new BigInteger("-2");
            }
            return x;
        }
    }
    public DbTool getDbTool() {
        return dbTool;
    }

    public MassTool returnMassTool() {
        return massTool;
    }

    public double getMinPeptideMass() {
        return minPeptideMass;
    }

    public double getMaxPeptideMass() {
        return maxPeptideMass;
    }

    public Map<Character, Double> returnFixModMap() {
        return fixModMap;
    }

    public InferSegment getInferSegment() {
        return inferSegment;
    }

    public InferPTM getInferPTM() {
        return inferPTM;
    }

    public TreeMap<Double, Set<String>> getMassPeptideMap() {
        return massPeptideMap;
    }

    public Map<String, Peptide0> getPeptide0Map() {
        return peptide0Map;
    }

    public String getLabelling() {
        return labelling;
    }

    private boolean isTarget(Collection<String> proteinIds) {
        for (String protein : proteinIds) {
            if (protein.startsWith("DECOY_")) { //????? why one decoy means it is not target
                return false;
            }
        }
        return true;
    }
}
