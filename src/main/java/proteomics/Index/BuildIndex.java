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

import java.io.*;
import java.util.*;

import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import proteomics.FM.FMIndex;
import proteomics.FM.SearchInterval;
import proteomics.FM.UnicodeReader;
import proteomics.PTM.InferPTM;
import proteomics.Segment.InferSegment;
import ProteomicsLibrary.*;
import ProteomicsLibrary.Types.*;
import proteomics.Types.Peptide0;

import static java.awt.SystemColor.text;

public class BuildIndex {

    private final MassTool massTool;
    private Map<Character, Double> fixModMap = new HashMap<>(25, 1);
    private double minPeptideMass = 9999;
    private double maxPeptideMass = 0;
    private InferSegment inferSegment;
    private TreeMap<Double, Set<String>> massPeptideMap = new TreeMap<>();
    private Map<String, Peptide0> peptide0Map;
    private final String labelling;
    private final DbTool dbTool; // this one doesn't contain contaminant proteins.
    private InferPTM inferPTM;
    public Map<String, Integer> protLengthMap = new HashMap<>();
    public FMIndex fmIndex;
    public int[] dotPosArr;
    public Map<Integer, String> posProtMap = new HashMap<>();
    public Map<String, String> proteinPeptideMap;

    public BuildIndex(Map<String, String> parameterMap) throws Exception {
        boolean needCoding = true;
        boolean addDecoy = parameterMap.get("add_decoy").contentEquals("1");
        boolean addContaminant = parameterMap.get("add_contaminant").contentEquals("1");
        int minPeptideLength = Math.max(5, Integer.valueOf(parameterMap.get("min_peptide_length")));
        int maxPeptideLength = Integer.valueOf(parameterMap.get("max_peptide_length"));
        String dbPath = parameterMap.get("db");
        int missedCleavage = Integer.valueOf(parameterMap.get("missed_cleavage"));
        double ms2Tolerance = Double.valueOf(parameterMap.get("ms2_tolerance"));
        double oneMinusBinOffset = 1 - Double.valueOf(parameterMap.get("mz_bin_offset"));
        this.labelling = parameterMap.get("15N").trim().contentEquals("1") ? "N15" : "N14";

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
//        Map<String, String> proteinPeptideMap;
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

        BufferedWriter writerProt = new BufferedWriter(new FileWriter("catProt.txt"));
//        String catProtStr = "";
        int dotPos = 0;
        int dotNum = 0;
        dotPosArr = new int[proteinPeptideMap.keySet().size()];
//        ArrayList<Integer>
        for (String proId : proteinPeptideMap.keySet()) {
            dotPosArr[dotNum] = dotPos;
            posProtMap.put(dotNum, proId);
            String proSeq = proteinPeptideMap.get(proId).replace('L', 'I');
            writerProt.write("." + proSeq.replace('L', 'I'));
            dotNum++;
            dotPos += proSeq.length()+1;
//            dotPosArr[dotNum] = dotPos;
            int numOfTags = inferSegment.generateSegmentBooleanVectorForProt(proSeq);
            protLengthMap.put(proId, numOfTags);
            Set<String> peptideSet = massTool.buildPeptideSetPnP(proSeq);

            if (proId.contentEquals("sp|P86452|ZBED6_HUMAN") || proId.contentEquals("sp|Q9UNZ5|L10K_HUMAN")) {
                int a = 1;
            }

            for (String peptide : peptideSet) {
                if (MassTool.containsNonAAAndNC(peptide)) {
                    continue;
                }

                if (peptide.contentEquals("nKIAIIKc") || peptide.contentEquals("nKIIAIKc")) {
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
                String decoyProSeq = DbTool.shuffleSeq(proSeq, parameterMap.get("cleavage_site_1"), parameterMap.get("protection_site_1"), Integer.valueOf(parameterMap.get("is_from_C_term_1")) == 1).replace('L', 'I'); // FixMe: Only consider the first enzyme if the users specify two enzymes.
                peptideSet = massTool.buildPeptideSetPnP(decoyProSeq);

                for (String peptide : peptideSet) {
                    if (MassTool.containsNonAAAndNC(peptide)) {
                        continue;
                    }

                    if (peptide.contentEquals("nKIAIIKc") || peptide.contentEquals("nKIIAIKc")) {
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
                            peptideProteinMap.put(peptide, "DECOY_" + proId);
                        }
                    }
                }
                targetDecoyProteinSequenceMap.put("DECOY_" + proId, decoyProSeq);
            }
        }
        writerProt.close();
        char[] text = loadFile("catProt.txt", true);
        fmIndex = new FMIndex(text);

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

        Map<String, Peptide0> tempMap = new HashMap<>();
        for (String peptide : peptideMassMap.keySet()) {
            if (peptide.contentEquals("nKIAIIKc")) {
                int a = 1;
            }
            SparseBooleanVector code = null;
            if (needCoding) {
                code = inferSegment.generateSegmentBooleanVector(peptide, peptideProteinMap.get(peptide).toArray(new String[0]));
            }

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
        peptide0Map = new HashMap<>(tempMap); // Since this map won't be changed any more, using this step to create a HashMap with the capacity exactly equals the actual size.
        int numTargetPep = 0, numDecoyPep = 0;
        for (Peptide0 pep : peptide0Map.values()){
            if (pep.isTarget) {
                numTargetPep++;
            }else {
                numDecoyPep++;
            }
        }
        System.out.println("Db size "+(numDecoyPep+ numTargetPep)+",targer : decoy = "+numTargetPep+" : "+numDecoyPep);
    }

    public static char[] loadFile(String file, boolean appendTerminalCharacter) throws IOException{
        // read text file, auto recognize bom marker or use
        // system default if markers not found.
        BufferedReader reader = null;
        CharArrayWriter writer = null;
        UnicodeReader r = new UnicodeReader(new FileInputStream(file), null);

        char[] buffer = new char[16 * 1024];   // 16k buffer
        int read;
        try {
            reader = new BufferedReader(r);
            writer = new CharArrayWriter();
            while( (read = reader.read(buffer)) != -1) {
                writer.write(buffer, 0, read);
            }
            if (appendTerminalCharacter) {
                writer.append('\0');
            }
            writer.flush();
            return writer.toCharArray();
        } catch (IOException ex) {
            throw ex;
        } finally {
            try {
                writer.close(); reader.close(); r.close();
            } catch (Exception ex) { }
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
