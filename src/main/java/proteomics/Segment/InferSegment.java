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

package proteomics.Segment;


import ProteomicsLibrary.DbTool;
import ProteomicsLibrary.MassTool;
import ProteomicsLibrary.Types.*;
import org.apache.commons.math3.util.Pair;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.Types.*;

import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class InferSegment {
    private static final Logger logger = LoggerFactory.getLogger(InferSegment.class);
    private static final int minTagNum = 200;
    private static final int regionNum = 10;
    private static final int topNumInEachRegion = 20;
    public static final byte N_TAG = -1;
    public static final byte NON_NC_TAG = 0;
    public static final byte C_TAG = 1;
    private static final Pattern pattern = Pattern.compile("([nc][0-9a-i])?([A-Z#$].?)");
    private final double ms2Tolerance;
    private TreeMap<Segment, Integer> aaVectorTemplate = new TreeMap<>();
    private Map<Double, String> augedMassAaMap = new HashMap<>(35, 1);// augmented amino acid map. Normal aa plus aa with mod
    private double maxAugedMass = 0;
    private final Double[] augedMassArray;
    private Map<String, Double> extraAaMassMap = new HashMap<>(35, 1);
    private double[] nTermPossibleMod = null;
    private double[] cTermPossibleMod = null;
    private MassTool massTool;
    public String aaModLabel = "~!@#$%^&*";
    public String nModLabel = "qwertyuio";
    public String cModLabel = "asdfghjkl";

    public InferSegment(MassTool massTool, Map<String, String> parameterMap, Map<Character, Double> fixModMap) throws Exception {
        this.massTool = massTool;
        this.ms2Tolerance = Double.valueOf(parameterMap.get("ms2_tolerance"));
        Map<Character, Double> massTable = massTool.getMassTable();

        char[] standardAaArray = new char[]{'G', 'A', 'S', 'P', 'V', 'T', 'C', 'L', 'N', 'D', 'Q', 'K', 'E', 'M', 'H', 'F', 'R', 'Y', 'W'};

        Map<Double, Character> oriMassAaMap = new HashMap<>(25, 1);
        for (char aa : standardAaArray) {
            // # = I/L.
            oriMassAaMap.put(massTable.get(aa), aa);
        }

        Character[] aaArray = oriMassAaMap.values().toArray(new Character[0]);

        for (char aa1 : aaArray) {
            for (char aa2 : aaArray) {
                for (char aa3 : aaArray) {
                    aaVectorTemplate.put(new Segment(String.format(Locale.US, "%c%c%c", aa1, aa2, aa3)), 0);
                }
            }
        }

        int idx = 0;
        for (Segment segment : aaVectorTemplate.keySet()) {
            aaVectorTemplate.put(segment, idx);
            ++idx;
        }

        // generate a mass aa map containing modified amino acid
        for (double k : oriMassAaMap.keySet()) {
            augedMassAaMap.put(k, oriMassAaMap.get(k).toString());
        }
        int iLabelAa = 0;
        int iLabelN = 0;
        int iLabelC = 0;
        for (String k : parameterMap.keySet()) {
            if (!k.startsWith("mod")) continue;

            String v = parameterMap.get(k);
            if (v.startsWith("0.0")) break;

            //15.994915,M,0,Oxidation
            String[] modStr = v.split(",");
            double modMass = Double.valueOf(modStr[0]);
            char modSite = modStr[1].charAt(0);
            char label;
            int modPosition = Integer.valueOf(modStr[2]);
            VarPtm varPtmByUser = new VarPtm(modMass, modSite, modPosition, modStr[3], "ByUser", 1);
            if (modPosition == 4) {//  position anywhere, highest prority
                //must happen in one aa, X (any aa) not allowed
                boolean added = false;
                label = aaModLabel.charAt(iLabelAa); // todo limit the max ilabel
                double augedMass = massTable.get(modSite) + modMass;
                boolean shouldAdd = true;
                for (double tempMass : augedMassAaMap.keySet()) {
                    if (Math.abs(tempMass - augedMass) <= ms2Tolerance) {
                        logger.warn(String.format(Locale.US, "%s is bad. %s and %s have conflict mass values(%f vs %f).", v, modSite, augedMassAaMap.get(tempMass), augedMass, tempMass));
                        shouldAdd = false;
                    }
                }
                if (Math.abs(fixModMap.get(modSite)) < 0.1 && shouldAdd) {
                    String augedAa = ""+modSite+label;
                    augedMassAaMap.put(augedMass, augedAa);
                    extraAaMassMap.put(augedAa, modMass);
                    added = true;
//                    massTool.labelMassMap.put(label, modMass);
                }

                if (added) {
                    iLabelAa++;
                    massTool.labelVarPtmMap.put(label, varPtmByUser);
                }

            } else if (modPosition == 0 || modPosition == 2) {// position N term, middle  prority // when find tags we can't differ pepN from protN
                label = nModLabel.charAt(iLabelN);
                boolean added = false;
                if (modSite == 'X') { //on any aa
                    for (char oriAa : aaArray){
                        double augedMass = massTable.get(oriAa) + modMass;
                        boolean shouldAdd = true;
                        for (double tempMass : augedMassAaMap.keySet()) {
                            if (Math.abs(tempMass - augedMass) <= ms2Tolerance) {
//                                        throw new Exception(String.format(Locale.US, "%s is bad. %s and %s have conflict mass values(%f vs %f).", v, oriAa, augedMassAaMap.get(tempMass), augedMass, tempMass));
                                logger.warn(String.format(Locale.US, "%s is bad. %s and %s have conflict mass values(%f vs %f).", v, oriAa, augedMassAaMap.get(tempMass), augedMass, tempMass));
                                shouldAdd = false;
                            }
                        }
                        if (Math.abs(fixModMap.get(oriAa)) < 0.1 && shouldAdd) {
                            // fix modification and var modification cannot be coexist
                            String augedAa = ""+oriAa+label;
                            augedMassAaMap.put(augedMass, augedAa);
                            extraAaMassMap.put(augedAa, modMass);
                            added = true;
                        }
                    }
                } else { //on single aa
                    double augedMass = massTable.get(modSite) + modMass;
                    boolean shouldAdd = true;
                    for (double tempMass : augedMassAaMap.keySet()) {
                        if (Math.abs(tempMass - augedMass) <= ms2Tolerance) {
                            logger.warn(String.format(Locale.US, "%s and %s have conflict mass values(%f vs %f).", v, augedMassAaMap.get(tempMass), augedMass, tempMass));
                            shouldAdd = false;
                        }
                    }
                    if (Math.abs(fixModMap.get(modSite)) < 0.1 && shouldAdd) {
                        // fix modification and var modification cannot be coexist
                        String augedAa = ""+modSite+label;
                        augedMassAaMap.put(augedMass, augedAa);
                        extraAaMassMap.put(augedAa, modMass);
                        added = true;
                    }
                }

                if (added) {
                    iLabelN++;
                    massTool.labelVarPtmMap.put(label, varPtmByUser);
                }
            } else {// position C term, lowest  prority
                label = cModLabel.charAt(iLabelC);
                boolean added = false;
                if (modSite == 'X') { //on any aa
                    for (char oriAa : aaArray){
                        double augedMass = massTable.get(oriAa) + modMass;
                        boolean shouldAdd = true;
                        for (double tempMass : augedMassAaMap.keySet()) {
                            if (Math.abs(tempMass - augedMass) <= ms2Tolerance) {
                                logger.warn(String.format(Locale.US, "%s is bad. %s and %s have conflict mass values(%f vs %f).", v, oriAa, augedMassAaMap.get(tempMass), augedMass, tempMass));
                                shouldAdd = false;
                            }
                        }
                        if (Math.abs(fixModMap.get(oriAa)) < 0.1 && shouldAdd) {
                            // fix modification and var modification cannot be coexist
                            String augedAa = ""+oriAa+label;
                            augedMassAaMap.put(augedMass, augedAa);
                            extraAaMassMap.put(augedAa, modMass);
                            added = true;
                        }
                    }
                } else {
                    double augedMass = massTable.get(modSite) + modMass;
                    boolean shouldAdd = true;
                    for (double tempMass : augedMassAaMap.keySet()) {
                        if (Math.abs(tempMass - augedMass) <= ms2Tolerance) {
                            logger.warn(String.format(Locale.US, "%s and %s have conflict mass values(%f vs %f).", v, augedMassAaMap.get(tempMass), augedMass, tempMass));
                            shouldAdd = false;
                        }
                    }
                    if (Math.abs(fixModMap.get(modSite)) < 0.1 && shouldAdd) {
                        // fix modification and var modification cannot be coexist
                        String augedAa = ""+modSite+label;
                        augedMassAaMap.put(augedMass, augedAa);
                        extraAaMassMap.put(augedAa, modMass);
                        added = true;
                    }
                }
                if (added) {
                    iLabelC++;
                    massTool.labelVarPtmMap.put(label, varPtmByUser);
                }
            }
        }
        augedMassArray = augedMassAaMap.keySet().toArray(new Double[0]);
        maxAugedMass = Collections.max(augedMassAaMap.keySet()) + ms2Tolerance;
    }

    public List<ExpTag> inferSegmentLocationFromSpectrum(double precursorMass, TreeMap<Double, Double> finalPlMap, int scanNum) throws Exception {
        return inferThreeAAFromSpectrum(finalPlMap, precursorMass - massTool.H2O + MassTool.PROTON);
    }

    public List<ExpTag> getAllTag4(double precursorMass, TreeMap<Double, Double> finalPlMap, int scanNum) throws Exception {
        return getTag4(finalPlMap, precursorMass - massTool.H2O + MassTool.PROTON, scanNum);
    }

    public List<ExpTag> getAllTag5(double precursorMass, TreeMap<Double, Double> finalPlMap, int scanNum) throws Exception {
        return getTag5(finalPlMap, precursorMass - massTool.H2O + MassTool.PROTON, scanNum);
    }

    public SparseVector generateSegmentIntensityVector(List<ExpTag> inputList) {
        SparseVector finalVector = new SparseVector();
        if (inputList.isEmpty()) {
            return finalVector;
        } else {
            for (ExpTag expAaList : inputList) {
                double totalIntensity = expAaList.getTotalIntensity();
                int idx = aaVectorTemplate.get(new Segment(expAaList.getFreeAaString()));
                double value = Math.max(totalIntensity, finalVector.get(idx));
                finalVector.put(idx, value);
            }
            return finalVector;
        }
    }

    private int minTagLen = 3;
    private int maxTagLen = 6;
    public Map<String, Double> getTagStrMap(List<ExpTag> longTagList) {
        Map<String, Double> finalMap = new HashMap<>();
        if (!longTagList.isEmpty()) {
            Segment tempSeg;
            for (ExpTag tag : longTagList) {
                String seq = tag.getFreeAaString();
                for (int i = minTagLen; i <= Math.min(tag.size(), maxTagLen); i++){
                    for (int j = 0; j <= tag.size()-i; j++){
                        tempSeg = new Segment(seq.substring(j, j+i));
                        finalMap.put(tempSeg.toString(), tag.getIntesFrom(j, j+i-1));
                    }
                }
            }
        }
        return finalMap;
    }

    public Set<String> getTagStrSet(String pepSeq) {
        Set<String> finalSet = new HashSet<>();
        Segment tempSeg;
        for (int i = minTagLen; i <= Math.min(pepSeq.length(), maxTagLen); i++){
            for (int j = 0; j <= pepSeq.length()-i; j++){
                tempSeg = new Segment(pepSeq.substring(j, j+i));
                finalSet.add(tempSeg.toString());
            }
        }
        return finalSet;
    }

    public double norm2square(Map<String, Double> scanTagStrMap){
        double res = 0;
        for (double value : scanTagStrMap.values()) {
            res += value*value;
        }
        return res;
    }
//    public double norm2square(Set<Double> candSeqSet){
//        double res = 0;
//        for (double value : scanTagStrMap.values()) {
//            res += value*value;
//        }
//        return res;
//    }


    public SparseBooleanVector generateSegmentBooleanVector(String peptide) {

        String normalizedPeptide = normalizeSequence(DbTool.getSequenceOnly(peptide));
//        for (int i = 0; i <= normalizedPeptide.length() - 4; ++i) {
//
//            Segment seg = new Segment(normalizedPeptide.substring(i, i + 4));
//            String tag = seg.toString();
//
//            if (tagPepMap.containsKey(tag)) {
//                tagPepMap.get(tag).add(peptide);
//            } else {
//                Set<String> pepSet = new HashSet<>();
//                pepSet.add(peptide);
//                tagPepMap.put(tag, pepSet);
//            }
//        }

        Set<Integer> tempSet = new HashSet<>(DbTool.getSequenceOnly(peptide).length() + 1, 1);
        for (int i = 0; i <= normalizedPeptide.length() - 3; ++i) {
            tempSet.add(aaVectorTemplate.get(new Segment(normalizedPeptide.substring(i, i + 3))));
        }
        return new SparseBooleanVector(tempSet);
    }
    public static String normalizeSequence(String seq) {
        return seq.replaceAll("I", "L");
    }

    private List<ExpTag> inferThreeAAFromSpectrum(TreeMap<Double, Double> plMap, double cTermMz) throws Exception {
        Double[] mzArray = plMap.keySet().toArray(new Double[0]);
        Double[] intensityArray = plMap.values().toArray(new Double[0]);
        Set<ExpTag> tempSet = new HashSet<>();
        List<ExpTag> outputList = new LinkedList<>();
        for (int i = 0; i < mzArray.length - 3; ++i) {
            double mz1 = mzArray[i];
            double intensity1 = intensityArray[i];
            for (int j = i + 1; j < mzArray.length - 2; ++j) {
                double mz2 = mzArray[j];
                double intensity2 = intensityArray[j];
                String aa1 = inferAA(mz1, mz2, 0);
                if (aa1 != null) {
                    Matcher matcher = pattern.matcher(aa1);
                    char ptmFreeAA = '\0';
                    double mod = 0;
                    double nTermMod = 0;
                    if (matcher.matches()) {
                        if (extraAaMassMap.containsKey(matcher.group(2))) {
                            mod = extraAaMassMap.get(matcher.group(2));
                        }
                        ptmFreeAA = matcher.group(2).charAt(0);
                        if (matcher.group(1) != null) {
                            if ((matcher.group(1).charAt(1) - '0' >= 0) && (matcher.group(1).charAt(1) - '0' < 10)) {
                                nTermMod = nTermPossibleMod[matcher.group(1).charAt(1) - '0'];
                            } else {
                                throw new Exception("Something is wrong in inferring tags.");
                            }
                        }
                    } else {
                        throw new NullPointerException(String.format(Locale.US, "Cannot find the PTM free amino acid for %s.", aa1));
                    }
                    ExpAa expAa1 = new ExpAa(aa1, ptmFreeAA, mz1, mz2, intensity1, intensity2, 0);
                    List<List<ExpAa>> tempAasList2 = new LinkedList<>();
                    for (int k = j + 1; k < mzArray.length - 1; ++k) {
                        double mz3 = mzArray[k];
                        double intensity3 = intensityArray[k];
                        String aa2 = inferAA(mz2, mz3, 0);
                        if (aa2 != null) {
                            mod = 0;
                            if (extraAaMassMap.containsKey(aa2)) {
                                mod = extraAaMassMap.get(aa2);
                            }
                            ExpAa expAa2 = new ExpAa(aa2, aa2.charAt(0), mz2, mz3, intensity2, intensity3, 0);
                            List<ExpAa> tempAasList3 = new LinkedList<>();
                            for (int l = k + 1; l < mzArray.length; ++l) {
                                double mz4 = mzArray[l];
                                double intensity4 = intensityArray[l];
                                String aa3 = inferAA(mz3, mz4, 0);
                                if (aa3 != null) {
                                    matcher = pattern.matcher(aa3);
                                    ptmFreeAA = '\0';
                                    mod = 0;
                                    double cTermMod = 0;
                                    if (matcher.matches()) {
                                        if (extraAaMassMap.containsKey(matcher.group(2))) {
                                            mod = extraAaMassMap.get(matcher.group(2));
                                        }
                                        ptmFreeAA = matcher.group(2).charAt(0);
                                        if (matcher.group(1) != null) {
                                            if ((matcher.group(1).charAt(1) - '0' >= 0) && (matcher.group(1).charAt(1) - '0' < 10)) {
                                                cTermMod = cTermPossibleMod[matcher.group(1).charAt(1) - '0'];
                                            } else {
                                                throw new Exception("Something is wrong in inferring tags.");
                                            }
                                        }
                                    } else {
                                        throw new NullPointerException(String.format(Locale.US, "Cannot find the PTM free amino acid for %s.", aa3));
                                    }
                                    ExpAa expAa3 = new ExpAa(aa3, ptmFreeAA, mz3, mz4, intensity3, intensity4, 0);
                                    tempAasList3.add(expAa3);
                                }
                            }
                            for (ExpAa expAas3 : tempAasList3) {
                                List<ExpAa> tempList2 = new LinkedList<>();
                                tempList2.add(expAa2);
                                tempList2.add(expAas3);
                                tempAasList2.add(tempList2);
                            }
                        }
                    }
                    for (List<ExpAa> expAas2 : tempAasList2) {
                        ExpTag expTag = new ExpTag(expAa1, expAas2.get(0), expAas2.get(1));
                        tempSet.add(expTag);
                    }
                }
            }
        }

        // eliminate "overlapped" tags
        ExpTag[] tempArray = tempSet.toArray(new ExpTag[0]);
        List<ExpTag> tempList = new LinkedList<>();
        for (int i = 0; i < tempArray.length; ++i) {
            boolean keep = true;
            for (int j = 0; j < tempArray.length; ++j) {
                if (i != j) {
                    if (tempArray[i].approximateEquals(tempArray[j], 2 * ms2Tolerance)) {
                        if (tempArray[i].getTotalIntensity() < tempArray[j].getTotalIntensity()) {
                            keep = false;
                            break;
                        }
                    }
                }
            }
            if (keep) {
                tempList.add(tempArray[i]);
            }
        }

        if (tempList.size() > minTagNum) {
            double minMz = plMap.firstKey();
            double regionWindow = Math.ceil((plMap.lastKey() - minMz) / regionNum);
            for (ExpTag expAa : tempList) {
                expAa.setRegionIdx((int) Math.floor((expAa.getHeadLocation() - minMz) / regionWindow));
            }
            List<List<Double>> regionIntensityList = new ArrayList<>(20);
            for (int i = 0; i < regionNum; ++i) {
                regionIntensityList.add(new ArrayList<>(100));
            }
            for (ExpTag expAa : tempList) {
                regionIntensityList.get(expAa.getRegionIdx()).add(expAa.getTotalIntensity());
            }
            double[] intensityTArray = new double[regionNum];
            for (int i = 0; i < regionNum; ++i) {
                List<Double> intensityList = regionIntensityList.get(i);
                Collections.sort(intensityList, Collections.reverseOrder());
                if (intensityList.size() > topNumInEachRegion) {
                    intensityTArray[i] = intensityList.get(topNumInEachRegion);
                }
            }
            for (ExpTag expAa : tempList) {
                if (expAa.getTotalIntensity() > intensityTArray[expAa.getRegionIdx()]) {
                    outputList.add(expAa);
                }
            }
            return outputList;
        } else {
            return tempList;
        }
    }
    public List<ExpTag> cleanAbundantTagsPrefix(List<ExpTag> allLongTagsList) {
        List<ExpTag> cleanTagList = new LinkedList<>();
        Set<String> prefixSet = new HashSet<>();
//        Set<String> suffixSet = new HashSet<>();
        for (ExpTag tag : allLongTagsList) {
            String prefix = tag.getFreeAaString().substring(0,4);
//            String suffix = tag.getFreeAaString().substring(tag.size()-4,tag.size());
            boolean shouldAdd = true;
            if (!prefixSet.contains(prefix)) {
                prefixSet.add(prefix);
            } else {
                shouldAdd = false;
            }
//            if (!suffixSet.contains(suffix)) {
//                suffixSet.add(suffix);
//            } else {
//                shouldAdd = false;
//            }
            if (shouldAdd) {
                cleanTagList.add(tag);
            }
        }
        return cleanTagList;
    }

    public List<ExpTag> cleanAbundantTags(List<ExpTag> allLongTagsList) {
        Map<String, List<ExpTag>> tagStrClusterMap = new HashMap<>();
        for (ExpTag tag : allLongTagsList.subList(0, allLongTagsList.size())){
            boolean shouldBeNewCluster = true;
            for (String clusterStr : tagStrClusterMap.keySet()) {
                String lcs = LCS(clusterStr, tag.getFreeAaString());
                double leveDist = minDistance(clusterStr, tag.getFreeAaString());
                if (lcs.length() >= tag.size()*0.5) {
//                if (leveDist < tag.size()*0.5) {

                    tagStrClusterMap.get(clusterStr).add(tag);
                    shouldBeNewCluster = false;
                    break; // only add into the first possible cluster
                } else {
                    continue;
                }
            }
            if (shouldBeNewCluster) {
                List<ExpTag> tmpList = new ArrayList<>();
                tmpList.add(tag);
                tagStrClusterMap.put(tag.getFreeAaString(), tmpList);
            }
        }
        List<ExpTag> cleanTagList = new LinkedList<>(); //todo what should I put into the list, cluster Name tags or the top tags in each cluster
        for (String clusterStr : tagStrClusterMap.keySet()) {
            List<ExpTag> clusterTagList = tagStrClusterMap.get(clusterStr);
            Collections.sort(clusterTagList, Comparator.comparing(o->o.getTotalIntensity(), Comparator.reverseOrder()));
            cleanTagList.add(clusterTagList.get(0));
        }
        Collections.sort(cleanTagList, Comparator.comparing(o->o.getTotalIntensity(), Comparator.reverseOrder()));
        return cleanTagList;
    }

    public List<ExpTag> cleanAbundantTagsNew(List<ExpTag> allLongTagsList) {
        Map<String, List<ExpTag>> tagStrClusterMap = new HashMap<>();
        for (ExpTag tag : allLongTagsList.subList(0, allLongTagsList.size())){
            boolean shouldBeNewCluster = true;
            for (String clusterStr : tagStrClusterMap.keySet()) {
                String lcs = LCS(clusterStr, tag.getFreeAaString());
                double leveDist = minDistance(clusterStr, tag.getFreeAaString());
                if (lcs.length() >= tag.size()*0.5) {
//                if (leveDist < tag.size()*0.5) {
                    tagStrClusterMap.get(clusterStr).add(tag);
//                    if (tagStrClusterMap.get(clusterStr).size() == 2) { // before add this tag it only has one tag, which means need to modify the tagstr key
//                        tagStrClusterMap.put(lcs,tagStrClusterMap.get(clusterStr));
//                        tagStrClusterMap.remove(clusterStr);
//                    }
                    shouldBeNewCluster = false;
                    break; // only add into the first possible cluster
                }
            }
            if (shouldBeNewCluster) {
                List<ExpTag> tmpList = new ArrayList<>();
                tmpList.add(tag);
                tagStrClusterMap.put(tag.getFreeAaString(), tmpList);
            }
        }
        List<ExpTag> cleanTagList = new LinkedList<>(); //todo what should I put into the list, cluster Name tags or the top tags in each cluster

        for (String clusterStr : tagStrClusterMap.keySet()) {
            List<ExpTag> clusterTagList = tagStrClusterMap.get(clusterStr);
            Collections.sort(clusterTagList, Comparator.comparing(o->o.getTotalIntensity(), Comparator.reverseOrder()));
            cleanTagList.add(clusterTagList.get(0));
            if (clusterTagList.size() > 1 && (clusterTagList.get(1).getTotalIntensity()>0.95*clusterTagList.get(0).getTotalIntensity())) {
                cleanTagList.add(clusterTagList.get(1));
            }
        }
        Collections.sort(cleanTagList, Comparator.comparing(o->o.getTotalIntensity(), Comparator.reverseOrder()));
        return cleanTagList;
    }
    private static int minDistance(String word1, String word2) {

        int dp[][] = new int[word1.length() + 1][word2.length() + 1];

        for (int i = 0; i < word1.length() + 1; i++) {
            // 从i个字符变成0个字符，需要i步（删除）
            dp[i][0] = i;
        }
        for (int i = 0; i < word2.length() + 1; i++) {
            // 当从0个字符变成i个字符，需要i步(增加)
            dp[0][i] = i;
        }

        for (int i = 1; i < word1.length() + 1; i++) {
            for (int j = 1; j < word2.length() + 1; j++) {
                //当相同的时，dp[i][j] = dp[i - 1][j - 1]
                if (word1.charAt(i - 1) == word2.charAt(j - 1)) {
                    dp[i][j] = dp[i - 1][j - 1];
                } else {
                    //当不同的时候，我们需要求三种操作的最小值
                    //其中dp[i - 1][j - 1]表示的是替换，dp[i - 1][j]表示删除字符，do[i][j - 1]表示的是增加字符
                    dp[i][j] = 1 + Math.min(dp[i - 1][j - 1], Math.min(dp[i - 1][j], dp[i][j - 1]));
                }
            }
        }
        return dp[word1.length()][word2.length()];
    }
    private static String LCS(String a, String b) {
        char[] arrayA = a.toCharArray();
        char[] arrayB = b.toCharArray();
        String maxString = "";
        for (int i = 0; i < arrayA.length; i++) {
            for (int j = 0; j < arrayB.length; j++) {
                if (arrayA[i]==arrayB[j]){
                    // i为起点
                    Integer i1 = i;
                    Integer j1 = j;
                    // 找到相同的字母后，以当前j为起点，持续比较下去，看相同部分到哪里结束
                    for (; i1 <arrayA.length && j1 < arrayB.length && arrayA[i1]==arrayB[j1]; i1++,j1++) {}
                    // i1为终点
                    if (i!=i1-1){
                        // 不断用common共同子串与max最长子串比较长度，只保留最长那个子串
                        String common = a.substring(i,i1);
                        if (maxString.length()<common.length()){
                            maxString = common;
                        }
                    }
                }
            }
        }
        return maxString;
    }
    public List<ExpTag> getLongTag(TreeMap<Double, Double> plMap, double cTermMz, int scanNum, int minTagLenToExtractLocal, int maxTagLenToExtractLocal) throws Exception {
        Double[] mzArray = plMap.keySet().toArray(new Double[0]);
        Double[] intensityArray = plMap.values().toArray(new Double[0]);
        List<ExpTag> outputList = new LinkedList<>();

        Set<Pair<Integer, Integer>> edgeSet = new HashSet<>();
        Set<Integer> nodeSet = new HashSet<>();
        Set<Integer> startNodeSet = IntStream.range(0, mzArray.length).boxed().collect(Collectors.toSet());
        Set<Integer> endNodeSet = IntStream.range(0, mzArray.length).boxed().collect(Collectors.toSet());
        Map<Pair<Integer, Integer>, ExpAa> edgeInfoMap = new HashMap<>();
        int numEdges = 0;
        double minEdgeWeight = 0.0;
        for (int i = 0; i < mzArray.length - 1; ++i) {
            double mz1 = mzArray[i];
            int isNorC = NON_NC_TAG;
            if (Math.abs(mz1 - MassTool.PROTON) <= ms2Tolerance) {
                isNorC = N_TAG;
            } else if (Math.abs(mz1 - MassTool.PROTON - massTool.H2O) <= ms2Tolerance) {
                isNorC = C_TAG;
            }
            double intensity1 = intensityArray[i];
            if (intensity1 < 0.1) continue;//todo

            for (int j = i + 1; j < mzArray.length; ++j) {
                double mz2 = mzArray[j];
                if (mz2 > mz1 + maxAugedMass + 5) break; // no need to try further because the can not match to any aa
                double intensity2 = intensityArray[j];
                if (intensity2 < 0.1) continue;//todo

                if ( (intensity1 + intensity2) < 0.4) continue;  //todo

                String aa = inferAA(mz1, mz2, isNorC);

                if (aa != null ) {
                    Pair<Integer, Integer> edge = new Pair<>(i,j);
                    edgeSet.add(edge);
                    nodeSet.add(i);
                    nodeSet.add(j);
                    startNodeSet.remove(j);
                    endNodeSet.remove(i);
                    edgeInfoMap.put(edge, new ExpAa(aa, aa.charAt(0), mz1, mz2, intensity1, intensity2, isNorC));
                }
            }
        }
        if (edgeInfoMap.size() > 120){ //only use the top 100 edges
            List<Map.Entry<Pair<Integer, Integer>, ExpAa>> edgeInfoList = new ArrayList<>(edgeInfoMap.entrySet());
            Collections.sort(edgeInfoList, Comparator.comparing(o -> o.getValue().getTotalIntensity(), Comparator.reverseOrder()));
            nodeSet.clear();
            edgeSet.clear();
            startNodeSet = IntStream.range(0, mzArray.length).boxed().collect(Collectors.toSet());
            endNodeSet = IntStream.range(0, mzArray.length).boxed().collect(Collectors.toSet());
            int i = 0;
            for (Map.Entry<Pair<Integer, Integer>, ExpAa> entry : edgeInfoList){
                if (i > 120) {
                    edgeInfoMap.remove(entry.getKey());
                } else {
                    edgeSet.add(entry.getKey());
                    nodeSet.add(entry.getKey().getFirst());
                    nodeSet.add(entry.getKey().getSecond());
                    startNodeSet.remove(entry.getKey().getSecond());
                    endNodeSet.remove(entry.getKey().getFirst());
                }
                i++;
            }

        }

        startNodeSet.retainAll(nodeSet);
        endNodeSet.retainAll(nodeSet);
        Graph g = new Graph(edgeSet, nodeSet);
        ArrayList<ArrayList<Integer>> allPath = g.getAllPaths(startNodeSet, endNodeSet);
        for (ArrayList<Integer> path : allPath) {
            if (path.size() < minTagLenToExtractLocal+1) continue; //aa length is 6 . peaks number is 7.
            if (path.size() > maxTagLenToExtractLocal + 1) {
                for (int i = 0; i < path.size()-maxTagLenToExtractLocal; i++) {
                    List<ExpAa> expAaList = new ArrayList<>();
                    for (int ii = i; ii < i+maxTagLenToExtractLocal; ii++){
                        int j = ii + 1;
                        expAaList.add(edgeInfoMap.get(new Pair<>(path.get(ii),path.get(j))));
                    }
                    outputList.add(new ExpTag(expAaList));
                }
            } else {
                List<ExpAa> expAaList = new ArrayList<>();
                for (int i = 0; i < path.size()-1; i++){
                    int j = i + 1;
                    expAaList.add(edgeInfoMap.get(new Pair<>(path.get(i),path.get(j))));
                }
                outputList.add(new ExpTag(expAaList));
            }
        }
        outputList.sort(Comparator.comparingDouble(ExpTag::getTotalIntensity).reversed());
        boolean[] shouldKeep = new boolean[outputList.size()];
        Set<String> addedTags = new HashSet<>();
        int i = 0;
        for(ExpTag tag : outputList) {
            if (addedTags.contains(tag.getFreeAaString())) {
                shouldKeep[i] = false;
            } else {
                shouldKeep[i] = true;
                addedTags.add(tag.getFreeAaString());
            }
            i++;
        }
        List<ExpTag> finalList = new LinkedList<>();
        for (int j = 0; j < outputList.size(); j++){
            if (shouldKeep[j]) {
                finalList.add(outputList.get(j));
            }
        }
        return finalList;
    }

    public List<ExpTag> getLongTagCutEnds(TreeMap<Double, Double> plMap, double cTermMz, int scanNum, int minTagLenToExtractLocal, int maxTagLenToExtractLocal) throws Exception {
        Double[] mzArray = plMap.keySet().toArray(new Double[0]);
        Double[] intensityArray = plMap.values().toArray(new Double[0]);
        List<ExpTag> outputList = new LinkedList<>();

        Set<Integer> nodeSet = new HashSet<>();
        Set<Integer> startNodeSet = new HashSet<>();
        Set<Integer> endNodeSet = new HashSet<>();
        Map<Pair<Integer, Integer>, ExpAa> edgeInfoMap = new HashMap<>();
        TreeMap<Integer, Double> nodeMap = new TreeMap<>();
        Map<Integer, Map<Integer, Double>> inEdgesMap = new HashMap<>();
        Map<Integer, Map<Integer, Double>> outEdgesMap = new HashMap<>();

        int maxEdgeNum = 120;
        Map<Double, Integer> mzIdMap = new HashMap<>();
        int ii = 0;
        for (double mz : plMap.keySet()) {
            mzIdMap.put(mz, ii);
            ii++;
        }

        for (int i = 0; i < mzArray.length - 1; ++i) {
            double mz1 = mzArray[i];
            int isNorC = NON_NC_TAG;
            if (Math.abs(mz1 - MassTool.PROTON) <= ms2Tolerance) {
                isNorC = N_TAG;
            } else if (Math.abs(mz1 - MassTool.PROTON - massTool.H2O) <= ms2Tolerance) {
                isNorC = C_TAG;
            }
            double intensity1 = intensityArray[i];
            if (intensity1 < 0.1) continue;//todo

            for (int j = i + 1; j < mzArray.length; ++j) {
                double mz2 = mzArray[j];
                if (mz2 > mz1 + maxAugedMass + 5) break; // no need to try further because the can not match to any aa
                double intensity2 = intensityArray[j];
                if (intensity2 < 0.1) continue;//todo

//                if ( (intensity1 + intensity2) < 0.4) continue;  //todo

                String aaStr = inferAA(mz1, mz2, isNorC);

                if (aaStr != null ) {
                    if ((mz1 == MassTool.PROTON + massTool.H2O) && (!"KR".contains(aaStr.substring(0,1)))) {
                        continue;
                    }
                    Pair<Integer, Integer> edge = new Pair<>(i,j);
                    ExpAa aa = new ExpAa(aaStr, aaStr.charAt(0), mz1, mz2, intensity1, intensity2, isNorC);
                    edgeInfoMap.put(edge, aa);
                    nodeMap.put(i, intensity1);
                    nodeMap.put(j, intensity2);
                    if (outEdgesMap.containsKey(i)) { // from i out to many j
                        outEdgesMap.get(i).put(j, 0.5 * (intensity1 + intensity2));
                    } else {
                        Map<Integer, Double> temp = new HashMap<>();
                        temp.put(j, 0.5 * (intensity1 + intensity2));
                        outEdgesMap.put(i, temp);
                    }
                    if (inEdgesMap.containsKey(j)) { // from many i in to j
                        inEdgesMap.get(j).put(i, 0.5 * (intensity1 + intensity2));
                    } else {
                        Map<Integer, Double> temp = new HashMap<>();
                        temp.put(i, 0.5 * (intensity1 + intensity2));
                        inEdgesMap.put(j, temp);
                    }
                }
            }
        }

        // re-calculate the weight of edge , because when edge is at the beginning , w should be intensity1 + 0.5*intensity2, not 0.5 * (intensity1 + intensity2)
        for (int node : nodeMap.keySet()) {
            if (!inEdgesMap.containsKey(node)) {
                for (int n2 : outEdgesMap.get(node).keySet()) {
                    outEdgesMap.get(node).put(n2, outEdgesMap.get(node).get(n2) + nodeMap.get(node) / 2);
                    inEdgesMap.get(n2).put(node, inEdgesMap.get(n2).get(node) + nodeMap.get(node) / 2);
                }
            }
            if (!outEdgesMap.containsKey(node)) {
                for (int n1 : inEdgesMap.get(node).keySet()) {
                    outEdgesMap.get(n1).put(node, outEdgesMap.get(n1).get(node) + nodeMap.get(node) / 2);
                    inEdgesMap.get(node).put(n1, inEdgesMap.get(node).get(n1) + nodeMap.get(node) / 2);
                }
            }
        }

        Set<Pair<Integer, Integer>> edgeToDel = new HashSet<>();

        for (int n1 : outEdgesMap.keySet()) {
            Map<Integer, Double> outNodesMap = outEdgesMap.get(n1);
            if (outNodesMap.size() == 1) continue; //no need to cut when there are no branches

            List<Integer> nOutToCut_List = new ArrayList<>();
            boolean hasLongOutEdge = false;
            for (int n2 : outNodesMap.keySet()) {
                if (outEdgesMap.containsKey(n2)) { //this out node has outnode, it is long
                    hasLongOutEdge = true;
                } else {
                    nOutToCut_List.add(n2);
                }
            }
            if (nOutToCut_List.isEmpty()) continue;
            Collections.sort(nOutToCut_List, Comparator.comparing(o->nodeMap.get(o)));
            for (int nOutToCut : nOutToCut_List.subList(0, nOutToCut_List.size()-1)) {
                edgeToDel.add(new Pair(n1, nOutToCut));
            }
            if (hasLongOutEdge) {
                edgeToDel.add(new Pair(n1, nOutToCut_List.get(nOutToCut_List.size()-1)));
            }
        }
        for (Pair<Integer, Integer> edge : edgeToDel){
            inEdgesMap.get(edge.getSecond()).remove(edge.getFirst());
            outEdgesMap.get(edge.getFirst()).remove(edge.getSecond());
            edgeInfoMap.remove(new Pair(edge.getFirst(), edge.getSecond()));
        }

        edgeToDel.clear();
        for (int n2 : inEdgesMap.keySet()) {
            Map<Integer, Double> inNodesMap = inEdgesMap.get(n2);
            if (inNodesMap.size() == 1) continue; //no need to cut when there are no branches

            List<Integer> nInToCut_List = new ArrayList<>();
            boolean hasLongInEdge = false;
            for (int n1 : inNodesMap.keySet()) {
                if (inEdgesMap.containsKey(n1)) { //this out node has outnode, it is long
                    hasLongInEdge = true;
                } else {
                    nInToCut_List.add(n1);
                }
            }
            if (nInToCut_List.isEmpty()) continue;
            Collections.sort(nInToCut_List, Comparator.comparing(o->nodeMap.get(o)));
            for (int nOutToCut : nInToCut_List.subList(0, nInToCut_List.size()-1)) {
                edgeToDel.add(new Pair(nOutToCut, n2));
            }
            if (hasLongInEdge) {
                edgeToDel.add(new Pair(nInToCut_List.get(nInToCut_List.size()-1), n2));
            }
        }
        for (Pair<Integer, Integer> edge : edgeToDel){
            inEdgesMap.get(edge.getSecond()).remove(edge.getFirst());
            outEdgesMap.get(edge.getFirst()).remove(edge.getSecond());
            edgeInfoMap.remove(new Pair(edge.getFirst(), edge.getSecond()));
        }

        if (edgeInfoMap.size() > maxEdgeNum){ //only use the top 100 edges
            List<Map.Entry<Pair<Integer, Integer>, ExpAa>> edgeInfoList = new ArrayList<>(edgeInfoMap.entrySet());
            Collections.sort(edgeInfoList, Comparator.comparing(o -> o.getValue().getTotalIntensity(), Comparator.reverseOrder()));
            int i = 0;
            for (Map.Entry<Pair<Integer, Integer>, ExpAa> entry : edgeInfoList){
                if (i > maxEdgeNum) {
                    Pair<Integer, Integer> edge = entry.getKey();
                    edgeInfoMap.remove(entry.getKey());
                    inEdgesMap.get(edge.getSecond()).remove(edge.getFirst());
                    outEdgesMap.get(edge.getFirst()).remove(edge.getSecond());
                }
                i++;
            }
        }
        startNodeSet.clear();
        endNodeSet.clear();
        for (int n2 : inEdgesMap.keySet()){
            if (!inEdgesMap.get(n2).isEmpty()) {
                if (!outEdgesMap.containsKey(n2)) {
                    endNodeSet.add(n2);
                } else {
                    if (outEdgesMap.get(n2).isEmpty()) {
                        endNodeSet.add(n2);
                    }
                }
            }
        }
        for (int n1 : outEdgesMap.keySet()){
            if (!outEdgesMap.get(n1).isEmpty()) {
                if (!inEdgesMap.containsKey(n1)) {
                    startNodeSet.add(n1);
                } else {
                    if (inEdgesMap.get(n1).isEmpty()) {
                        startNodeSet.add(n1);
                    }
                }
            }
        }
        for (Pair<Integer, Integer> edge : edgeInfoMap.keySet()){
            nodeSet.add(edge.getFirst());
            nodeSet.add(edge.getSecond());
        }
        Graph g = new Graph(edgeInfoMap.keySet(), nodeSet);

        ArrayList<ArrayList<Integer>> allPath = g.getAllPaths(startNodeSet, endNodeSet);
        for (ArrayList<Integer> path : allPath) {
            if (path.size() < minTagLenToExtractLocal+1) continue; //aa length is 6 . peaks number is 7.
            if (path.size() > maxTagLenToExtractLocal + 1) {
                for (int i = 0; i < path.size()-maxTagLenToExtractLocal; i++) {
                    List<ExpAa> expAaList = new ArrayList<>();
                    for (int k = i; k < i+maxTagLenToExtractLocal; k++){
                        int j = k + 1;
                        expAaList.add(edgeInfoMap.get(new Pair<>(path.get(k),path.get(j))));
                    }
                    outputList.add(new ExpTag(expAaList));
                }
            } else {
                List<ExpAa> expAaList = new ArrayList<>();
                for (int i = 0; i < path.size()-1; i++){
                    int j = i + 1;
                    expAaList.add(edgeInfoMap.get(new Pair<>(path.get(i),path.get(j))));
                }
                outputList.add(new ExpTag(expAaList));
            }
        }
        outputList.sort(Comparator.comparingDouble(ExpTag::getTotalIntensity).reversed());
        boolean[] shouldKeep = new boolean[outputList.size()];
        Set<String> addedTags = new HashSet<>();
        int i = 0;
        for(ExpTag tag : outputList) {
            if (addedTags.contains(tag.getFreeAaString())) {
                shouldKeep[i] = false;
            } else {
                shouldKeep[i] = true;
                addedTags.add(tag.getFreeAaString());
            }
            i++;
        }
        List<ExpTag> finalList = new LinkedList<>();
        for (int j = 0; j < outputList.size(); j++){
            if (shouldKeep[j]) {
                finalList.add(outputList.get(j));
            }
        }
        return finalList;
    }
    public List<ExpTag> getLongDenoisedTag(TreeMap<Double, Double> plMap, double cTermMz, int scanNum, Set<Pair<Integer, Integer>> edge1ToDel) throws Exception {
        Double[] mzArray = plMap.keySet().toArray(new Double[0]);
        Double[] intensityArray = plMap.values().toArray(new Double[0]);
        List<ExpTag> outputList = new LinkedList<>();

        Set<Pair<Integer, Integer>> edgeSet = new HashSet<>();
        Set<Integer> nodeSet = new HashSet<>();
        Set<Integer> startNodeSet = IntStream.range(0, mzArray.length).boxed().collect(Collectors.toSet());
        Set<Integer> endNodeSet = IntStream.range(0, mzArray.length).boxed().collect(Collectors.toSet());
        Map<Pair<Integer, Integer>, ExpAa> edgeInfoMap = new HashMap<>();

        Map<Double, Integer> mzIdMap = new HashMap<>();
        int ii = 0;
        for (double mz : plMap.keySet()) {
            mzIdMap.put(mz, ii);
            ii++;
        }
        TreeMap<Integer, Double> nodeMap = new TreeMap<>();
        Map<Integer, Map<Integer, Double>> inEdgeMap = new HashMap<>();
        Map<Integer, Map<Integer, Double>> outEdgeMap = new HashMap<>();

        for (int i = 0; i < mzArray.length - 1; ++i) {
            double mz1 = mzArray[i];
            double intensity1 = intensityArray[i];
            for (int j = i + 1; j < mzArray.length; ++j) {
                double mz2 = mzArray[j];
                double intensity2 = intensityArray[j];
                int isNorC = NON_NC_TAG;
                if (Math.abs(mz1 - MassTool.PROTON) <= ms2Tolerance) {
                    isNorC = N_TAG;
                } else if (Math.abs(mz1 - MassTool.PROTON - massTool.H2O) <= ms2Tolerance) {
                    isNorC = C_TAG;
                }
                String aaStr = inferAA(mz1, mz2, isNorC);
                if (aaStr != null) {
                    Pair<Integer, Integer> edge = new Pair<>(i,j);
                    nodeSet.add(i);
                    nodeSet.add(j);
                    startNodeSet.remove(j);
                    endNodeSet.remove(i);
                    ExpAa aa = new ExpAa(aaStr, aaStr.charAt(0), mz1, mz2, intensity1, intensity2, 0);
                    edgeInfoMap.put(edge, aa);

                    if ((aa.getHeadLocation() == MassTool.PROTON + massTool.H2O) && (!"KR".contains(aa.getAA()))) {
                        continue;
                    }
                    nodeMap.put(i, intensity1);
                    nodeMap.put(j, intensity2);
                    if (outEdgeMap.containsKey(i)) {
                        outEdgeMap.get(i).put(j, 0.5 * (intensity1 + intensity2));
                    } else {
                        Map<Integer, Double> temp = new HashMap<>();
                        temp.put(j, 0.5 * (intensity1 + intensity2));
                        outEdgeMap.put(i, temp);
                    }

                    if (inEdgeMap.containsKey(j)) {
                        inEdgeMap.get(j).put(i, 0.5 * (intensity1 + intensity2));
                    } else {
                        Map<Integer, Double> temp = new HashMap<>();
                        temp.put(i, 0.5 * (intensity1 + intensity2));
                        inEdgeMap.put(j, temp);
                    }

                    if ((aa.getHeadLocation() == MassTool.PROTON + massTool.H2O) && (!"KR".contains(aa.getAA()))) {
                        continue;
                    }
                    nodeMap.put(i, intensity1);
                    nodeMap.put(j, intensity2);
                    if (outEdgeMap.containsKey(i)) {
                        outEdgeMap.get(i).put(j, 0.5 * (intensity1 + intensity2));
                    } else {
                        Map<Integer, Double> temp = new HashMap<>();
                        temp.put(j, 0.5 * (intensity1 + intensity2));
                        outEdgeMap.put(i, temp);
                    }

                    if (inEdgeMap.containsKey(j)) {
                        inEdgeMap.get(j).put(i, 0.5 * (intensity1 + intensity2));
                    } else {
                        Map<Integer, Double> temp = new HashMap<>();
                        temp.put(i, 0.5 * (intensity1 + intensity2));
                        inEdgeMap.put(j, temp);
                    }
                }
            }
        }

        for (int node : nodeMap.keySet()) {
            if (!inEdgeMap.containsKey(node)) {
                for (int n2 : outEdgeMap.get(node).keySet()) {
                    outEdgeMap.get(node).put(n2, outEdgeMap.get(node).get(n2) + nodeMap.get(node) / 2);
                    inEdgeMap.get(n2).put(node, inEdgeMap.get(n2).get(node) + nodeMap.get(node) / 2);
                }
            }
            if (!outEdgeMap.containsKey(node)) {
                for (int n1 : inEdgeMap.get(node).keySet()) {
                    outEdgeMap.get(n1).put(node, outEdgeMap.get(n1).get(node) + nodeMap.get(node) / 2);
                    inEdgeMap.get(node).put(n1, inEdgeMap.get(node).get(n1) + nodeMap.get(node) / 2);
                }
            }
        }
        // finish prepare graph

        Set<Pair<Integer, Integer>> edgeToDel = new HashSet<>();
        Map<Integer, Integer> tempToDel = new HashMap<>();

        for (int node : nodeMap.descendingKeySet()) {
            if (outEdgeMap.containsKey(node)) {
                Map<Integer, Double> tempOutMap = outEdgeMap.get(node);
                double maxOut = 0d;
                int maxOutPeak = -1;
                for (int n2 : tempOutMap.keySet()) {
                    if (tempOutMap.get(n2) > maxOut) {
                        maxOut = tempOutMap.get(n2);
                        maxOutPeak = n2;
                    }
                }
                for (int n2 : tempOutMap.keySet()) {
                    if (n2 != maxOutPeak) {
                        edgeToDel.add(new Pair(node, n2));
                        tempToDel.put(node, n2);

                    }
                }
                if (inEdgeMap.containsKey(node)) {

                    for (int n1 : inEdgeMap.get(node).keySet()) {
                        outEdgeMap.get(n1).put(node, outEdgeMap.get(n1).get(node) + maxOut);
                        inEdgeMap.get(node).put(n1, inEdgeMap.get(node).get(n1) + maxOut);

                    }
                }
            }
        }

        for (int n1 : tempToDel.keySet()) {
            int n2 = tempToDel.get(n1);
            inEdgeMap.get(n2).remove(n1);
            outEdgeMap.get(n1).remove(n2);
        }
        for (int node : nodeMap.keySet()) {
            if (inEdgeMap.containsKey(node)) {
                Map<Integer, Double> tempInMap = inEdgeMap.get(node);
                double maxIn = 0d;
                int maxInPeak = -1;
                for (int n1 : tempInMap.keySet()) {
                    if (tempInMap.get(n1) > maxIn) {
                        maxIn = tempInMap.get(n1);
                        maxInPeak = n1;
                    }
                }
                for (int n1 : tempInMap.keySet()) {
                    if (n1 != maxInPeak) {
                        edgeToDel.add(new Pair(n1, node));

                    }
                }
                if (outEdgeMap.containsKey(node)) {

                    for (int n2 : outEdgeMap.get(node).keySet()) {
                        outEdgeMap.get(node).put(n2, outEdgeMap.get(node).get(n2) + maxIn);
                        inEdgeMap.get(n2).put(node, inEdgeMap.get(n2).get(node) + maxIn);

                    }
                }
            }
        }
        Set<Pair<Integer, Integer>> pairToDel = new HashSet<>();
        for (Pair<Integer, Integer> e : edgeToDel) {
            pairToDel.add(new Pair(e.getFirst(), e.getSecond()));
        }

//        if (!edgeToDel.contains(edge)) edgeSet.add(edge);
        startNodeSet.retainAll(nodeSet);
        endNodeSet.retainAll(nodeSet);

        for (Pair<Integer, Integer> edge : edgeToDel) {

        }



        Graph g = new Graph(edgeSet, nodeSet);
        ArrayList<ArrayList<Integer>> allPath = g.getAllPaths(startNodeSet, endNodeSet);
//        Map<String, Double>
        for (ArrayList<Integer> path : allPath) {
            if (path.size() < 7) continue; //aa length is 6 . peaks number is 7.
            List<ExpAa> expAaList = new ArrayList<>();
            for (int i = 0; i < path.size()-1; i++){
                int j = i + 1;
                expAaList.add(edgeInfoMap.get(new Pair<>(path.get(i),path.get(j))));
            }
            outputList.add(new ExpTag(expAaList));
        }
        outputList.sort(Comparator.comparingDouble(ExpTag::getTotalIntensity).reversed());
        boolean[] shouldKeep = new boolean[outputList.size()];
        Set<String> addedTags = new HashSet<>();
        int i = 0;
        for(ExpTag tag : outputList) {
            if (addedTags.contains(tag.getFreeAaString())) {
                shouldKeep[i] = false;
            } else {
                shouldKeep[i] = true;
                addedTags.add(tag.getFreeAaString());
            }
            i++;
        }
        List<ExpTag> finalList = new LinkedList<>();
        for (int j = 0; j < outputList.size(); j++){
            if (shouldKeep[j]) {
                finalList.add(outputList.get(j));
            }
        }
        return finalList;
    }
    public int getLongTagNumForProt(String prot) {
        String normalizedProt = normalizeSequence(DbTool.getSequenceOnly(prot));
        Set<String> tagSet = new HashSet<>();
        for (int i = 0; i <= normalizedProt.length() - 7; ++i) {
            Segment seg = new Segment(normalizedProt.substring(i, i + 7));
            String tag = seg.toString();
            tagSet.add(tag);
//            Segment segL = new Segment(normalizedProt.substring(i, i + 4).replace('#','L'));
//            String tagL = segL.toString();
        }

        return tagSet.size();
    }
    private List<ExpTag> getTag4(TreeMap<Double, Double> plMap, double cTermMz, int scanNum) throws Exception {
        Double[] mzArray = plMap.keySet().toArray(new Double[0]);
        Double[] intensityArray = plMap.values().toArray(new Double[0]);
        Set<ExpTag> tempSet = new HashSet<>();
        List<ExpTag> outputList = new LinkedList<>();
        for (int i = 0; i < mzArray.length - 4; ++i) {
            double mz1 = mzArray[i];
            double intensity1 = intensityArray[i];
            for (int j = i + 1; j < mzArray.length - 3; ++j) {
                double mz2 = mzArray[j];
                double intensity2 = intensityArray[j];
                if ((intensity1+intensity2)/2 < 0.2) continue;
                String aa1 = inferAA(mz1, mz2, 0);
                if (aa1 != null) {
                    Matcher matcher = pattern.matcher(aa1);
                    char ptmFreeAA = '\0';
                    double mod = 0;
                    double nTermMod = 0;
                    if (matcher.matches()) {
                        if (extraAaMassMap.containsKey(matcher.group(2))) {
                            mod = extraAaMassMap.get(matcher.group(2));
                        }
                        ptmFreeAA = matcher.group(2).charAt(0);
                        if (matcher.group(1) != null) {
                            if ((matcher.group(1).charAt(1) - '0' >= 0) && (matcher.group(1).charAt(1) - '0' < 10)) {
                                nTermMod = nTermPossibleMod[matcher.group(1).charAt(1) - '0'];
                            } else {
                                throw new Exception("Something is wrong in inferring tags.");
                            }
                        }
                    } else {
                        throw new NullPointerException(String.format(Locale.US, "Cannot find the PTM free amino acid for %s.", aa1));
                    }
                    ExpAa expAa1 = new ExpAa(aa1, ptmFreeAA, mz1, mz2, intensity1, intensity2, 0);
                    List<List<ExpAa>> tempAasList2 = new LinkedList<>();
                    for (int k = j + 1; k < mzArray.length - 2; ++k) {
                        double mz3 = mzArray[k];
                        double intensity3 = intensityArray[k];
                        if ((intensity2+intensity3)/2 < 0.2) continue;
                        String aa2 = inferAA(mz2, mz3, 0);
                        if (aa2 != null) {
                            mod = 0;
                            if (extraAaMassMap.containsKey(aa2)) {
                                mod = extraAaMassMap.get(aa2);
                            }
                            ExpAa expAa2 = new ExpAa(aa2, aa2.charAt(0), mz2, mz3, intensity2, intensity3,  0);
                            List<ExpAa> tempAasList4 = new LinkedList<>();
                            for (int l = k + 1; l < mzArray.length - 1; ++l) {
                                double mz4 = mzArray[l];
                                double intensity4 = intensityArray[l];
                                if ((intensity3+intensity4)/2 < 0.2) continue;
                                String aa3 = inferAA(mz3, mz4, 0);
                                if (aa3 != null) {
                                    mod = 0;
                                    if (extraAaMassMap.containsKey(matcher.group(2))) {
                                        mod = extraAaMassMap.get(matcher.group(2));
                                    }
                                    ExpAa expAa3 = new ExpAa(aa3, aa3.charAt(0), mz3, mz4, intensity3, intensity4,  0);
                                    tempAasList4.add(expAa3);
                                    for (int m = l + 1; m < mzArray.length; ++m) {
                                        double mz5 = mzArray[m];
                                        double intensity5 = intensityArray[m];
                                        if ((intensity4+intensity5)/2 < 0.2) continue;
                                        String aa4 = inferAA(mz4, mz5, 0);
                                        if (aa4 != null) {
                                            matcher = pattern.matcher(aa4);
                                            ptmFreeAA = '\0';
                                            mod = 0;
                                            double cTermMod = 0;
                                            if (matcher.matches()) {
                                                if (extraAaMassMap.containsKey(matcher.group(2))) {
                                                    mod = extraAaMassMap.get(matcher.group(2));
                                                }
                                                ptmFreeAA = matcher.group(2).charAt(0);
                                                if (matcher.group(1) != null) {
                                                    if ((matcher.group(1).charAt(1) - '0' >= 0) && (matcher.group(1).charAt(1) - '0' < 10)) {
                                                        cTermMod = cTermPossibleMod[matcher.group(1).charAt(1) - '0'];
                                                    } else {
                                                        throw new Exception("Something is wrong in inferring tags.");
                                                    }
                                                }
                                            } else {
                                                throw new NullPointerException(String.format(Locale.US, "Cannot find the PTM free amino acid for %s.", aa4));
                                            }
                                            ExpAa expAa4 = new ExpAa(aa4, ptmFreeAA, mz4, mz5, intensity4, intensity5, 0);
                                            tempAasList4.add(expAa4);
                                            tempSet.add(new ExpTag(expAa1, expAa2, expAa3, expAa4));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        // eliminate "overlapped" tags
        ExpTag[] tempArray = tempSet.toArray(new ExpTag[0]);
        List<ExpTag> tempList = new LinkedList<>();
        for (int i = 0; i < tempArray.length; ++i) {
            boolean keep = true;
            for (int j = 0; j < tempArray.length; ++j) {
                if (i != j) {
                    if (tempArray[i].approximateEquals(tempArray[j], 2 * ms2Tolerance)) {
                        if (tempArray[i].getTotalIntensity() < tempArray[j].getTotalIntensity()) {
                            keep = false;
                            break;
                        }
                    }
                }
            }
            if (keep) {
                tempList.add(tempArray[i]);
            }
        }

        if (tempList.size() > minTagNum) {
            double minMz = plMap.firstKey();
            double regionWindow = Math.ceil((plMap.lastKey() - minMz) / regionNum);
            for (ExpTag expAa : tempList) {
                expAa.setRegionIdx((int) Math.floor((expAa.getHeadLocation() - minMz) / regionWindow));
            }
            List<List<Double>> regionIntensityList = new ArrayList<>(20);
            for (int i = 0; i < regionNum; ++i) {
                regionIntensityList.add(new ArrayList<>(100));
            }
            for (ExpTag expAa : tempList) {
                regionIntensityList.get(expAa.getRegionIdx()).add(expAa.getTotalIntensity());
            }
            double[] intensityTArray = new double[regionNum];
            for (int i = 0; i < regionNum; ++i) {
                List<Double> intensityList = regionIntensityList.get(i);
                Collections.sort(intensityList, Collections.reverseOrder());
                if (intensityList.size() > topNumInEachRegion) {
                    intensityTArray[i] = intensityList.get(topNumInEachRegion);
                }
            }
            for (ExpTag expAa : tempList) {
                if (expAa.getTotalIntensity() > intensityTArray[expAa.getRegionIdx()]) {
                    outputList.add(expAa);
                }
            }
            return outputList;
        } else {
            return tempList;
        }
    }

    private List<ExpTag> getTag5(TreeMap<Double, Double> plMap, double cTermMz, int scanNum) throws Exception {
        Double[] mzArray = plMap.keySet().toArray(new Double[0]);
        Double[] intensityArray = plMap.values().toArray(new Double[0]);
        Set<ExpTag> tempSet = new HashSet<>();
        List<ExpTag> outputList = new LinkedList<>();
        for (int i = 0; i < mzArray.length - 5; ++i) {
            double mz1 = mzArray[i];
            double intensity1 = intensityArray[i];
            for (int j = i + 1; j < mzArray.length - 4; ++j) {
                double mz2 = mzArray[j];
                double intensity2 = intensityArray[j];
                if ((intensity1+intensity2)/2 < 0.5) continue;
                String aa1 = inferAA(mz1, mz2, 0);
                if (aa1 != null) {
//                    Matcher matcher = pattern.matcher(aa1);
//                    char ptmFreeAA = '\0';
//                    double mod = 0;
//                    double nTermMod = 0;
//                    if (matcher.matches()) {
//                        if (modifiedAAMassMap.containsKey(matcher.group(2))) {
//                            mod = modifiedAAMassMap.get(matcher.group(2));
//                        }
//                        ptmFreeAA = matcher.group(2).charAt(0);
//                        if (matcher.group(1) != null) {
//                            if ((matcher.group(1).charAt(1) - '0' >= 0) && (matcher.group(1).charAt(1) - '0' < 10)) {
//                                nTermMod = nTermPossibleMod[matcher.group(1).charAt(1) - '0'];
//                            } else {
//                                throw new Exception("Something is wrong in inferring tags.");
//                            }
//                        }
//                    } else {
//                        throw new NullPointerException(String.format(Locale.US, "Cannot find the PTM free amino acid for %s.", aa1));
//                    }
                    ExpAa expAa1 = new ExpAa(aa1, aa1.charAt(0), mz1, mz2, intensity1, intensity2, 0);
                    for (int k = j + 1; k < mzArray.length - 3; ++k) {
                        double mz3 = mzArray[k];
                        double intensity3 = intensityArray[k];
                        if ((intensity2+intensity3)/2 < 0.5) continue;
                        String aa2 = inferAA(mz2, mz3, 0);
                        if (aa2 != null) {
//                            mod = 0;
//                            if (modifiedAAMassMap.containsKey(aa2)) {
//                                mod = modifiedAAMassMap.get(aa2);
//                            }
                            ExpAa expAa2 = new ExpAa(aa2, aa2.charAt(0), mz2, mz3, intensity2, intensity3, 0);
                            for (int l = k + 1; l < mzArray.length - 2; ++l) {
                                double mz4 = mzArray[l];
                                double intensity4 = intensityArray[l];
                                if ((intensity3+intensity4)/2 < 0.5) continue;
                                String aa3 = inferAA(mz3, mz4, 0);
                                if (aa3 != null) {
//                                    mod = 0;
//                                    if (modifiedAAMassMap.containsKey(matcher.group(2))) {
//                                        mod = modifiedAAMassMap.get(matcher.group(2));
//                                    }
                                    ExpAa expAa3 = new ExpAa(aa3, aa3.charAt(0), mz3, mz4, intensity3, intensity4, 0);
                                    for (int m = l + 1; m < mzArray.length; ++m) {
                                        double mz5 = mzArray[m];
                                        double intensity5 = intensityArray[m];
                                        if ((intensity4+intensity5)/2 < 0.5) continue;
                                        String aa4 = inferAA(mz4, mz5, 0);
                                        if (aa4 != null) {
//                                            mod = 0;
//                                            if (modifiedAAMassMap.containsKey(matcher.group(2))) {
//                                                mod = modifiedAAMassMap.get(matcher.group(2));
//                                            }
                                            ExpAa expAa4 = new ExpAa(aa4, aa4.charAt(0), mz4, mz5, intensity4, intensity5, 0);
                                            for (int n = m + 1; n < mzArray.length; ++n) {
                                                double mz6 = mzArray[n];
                                                double intensity6 = intensityArray[n];
                                                if ((intensity5+intensity6)/2 < 0.5) continue;
                                                String aa5 = inferAA(mz5, mz6, 0);
                                                if (aa5 != null) {
//                                                    matcher = pattern.matcher(aa5);
//                                                    ptmFreeAA = '\0';
//                                                    mod = 0;
//                                                    double cTermMod = 0;
//                                                    if (matcher.matches()) {
//                                                        if (modifiedAAMassMap.containsKey(matcher.group(2))) {
//                                                            mod = modifiedAAMassMap.get(matcher.group(2));
//                                                        }
//                                                        ptmFreeAA = matcher.group(2).charAt(0);
//                                                        if (matcher.group(1) != null) {
//                                                            if ((matcher.group(1).charAt(1) - '0' >= 0) && (matcher.group(1).charAt(1) - '0' < 10)) {
//                                                                cTermMod = cTermPossibleMod[matcher.group(1).charAt(1) - '0'];
//                                                            } else {
//                                                                throw new Exception("Something is wrong in inferring tags.");
//                                                            }
//                                                        }
//                                                    } else {
//                                                        throw new NullPointerException(String.format(Locale.US, "Cannot find the PTM free amino acid for %s.", aa4));
//                                                    }
                                                    ExpAa expAa5 = new ExpAa(aa5, aa5.charAt(0), mz5, mz6, intensity5, intensity6, 0);
                                                    tempSet.add(new ExpTag(expAa1, expAa2, expAa3, expAa4, expAa5));
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

        // eliminate "overlapped" tags
        ExpTag[] tempArray = tempSet.toArray(new ExpTag[0]);
        List<ExpTag> tempList = new LinkedList<>();
        for (int i = 0; i < tempArray.length; ++i) {
            boolean keep = true;
            for (int j = 0; j < tempArray.length; ++j) {
                if (i != j) {
                    if (tempArray[i].approximateEquals(tempArray[j], 2 * ms2Tolerance)) {
                        if (tempArray[i].getTotalIntensity() < tempArray[j].getTotalIntensity()) {
                            keep = false;
                            break;
                        }
                    }
                }
            }
            if (keep) {
                tempList.add(tempArray[i]);
            }
        }

        if (tempList.size() > minTagNum) {
            double minMz = plMap.firstKey();
            double regionWindow = Math.ceil((plMap.lastKey() - minMz) / regionNum);
            for (ExpTag expAa : tempList) {
                expAa.setRegionIdx((int) Math.floor((expAa.getHeadLocation() - minMz) / regionWindow));
            }
            List<List<Double>> regionIntensityList = new ArrayList<>(20);
            for (int i = 0; i < regionNum; ++i) {
                regionIntensityList.add(new ArrayList<>(100));
            }
            for (ExpTag expAa : tempList) {
                regionIntensityList.get(expAa.getRegionIdx()).add(expAa.getTotalIntensity());
            }
            double[] intensityTArray = new double[regionNum];
            for (int i = 0; i < regionNum; ++i) {
                List<Double> intensityList = regionIntensityList.get(i);
                Collections.sort(intensityList, Collections.reverseOrder());
                if (intensityList.size() > topNumInEachRegion) {
                    intensityTArray[i] = intensityList.get(topNumInEachRegion);
                }
            }
            for (ExpTag expAa : tempList) {
                if (expAa.getTotalIntensity() > intensityTArray[expAa.getRegionIdx()]) {
                    outputList.add(expAa);
                }
            }
            return outputList;
        } else {
            return tempList;
        }
    }
    private String inferAA(double mz1, double mz2, int isNorC) {
        double mzDiff = mz2 - mz1;
        String aa = null;
        for (Map.Entry<Double, String> massAa : augedMassAaMap.entrySet()) {
            String augedAa = massAa.getValue();
            if ((isNorC == N_TAG && cModLabel.contains(augedAa.substring(augedAa.length()-1)))  // use n peak then should no C mod
                    || (isNorC == C_TAG && nModLabel.contains(augedAa.substring(augedAa.length()-1))) // use c peak then should no N mod
                    || (isNorC == NON_NC_TAG && cModLabel.contains(augedAa.substring(augedAa.length()-1))) // use non-nc peak then should no C mod
                    || (isNorC == NON_NC_TAG && nModLabel.contains(augedAa.substring(augedAa.length()-1))) // use non-nc peak then should no N mod
            ) {
                continue;
            }

            if (Math.abs(mzDiff - massAa.getKey()) <= 2 * ms2Tolerance) {
                aa = massAa.getValue();// including M~
            }
        }

        if (aa != null && isNorC == C_TAG && (!"KR".contains(aa))) {
            aa = null;   // if it is C tag, it must start with KR
        }

        if (aa != null && aa.length() > 1) {
            int a = 1;
        }
        return aa;
    }

    public TreeMap<Double, Double> addVirtualPeaks(double precursorMass, TreeMap<Double, Double> plMap) {
        double totalMass = precursorMass + 2 * MassTool.PROTON;
        TreeMap<Double, Double> finalPlMap = new TreeMap<>();
        for (double mz : plMap.keySet()) {
            finalPlMap.put(mz, plMap.get(mz));
        }
        for (double mz : plMap.keySet()) {
            double anotherMz = totalMass - mz;
            double leftMz = anotherMz - ms2Tolerance;
            double rightMz = anotherMz + ms2Tolerance;
            NavigableMap<Double, Double> temp = null;
            try {
                temp = plMap.subMap(leftMz, true, rightMz, true);
            } catch (IllegalArgumentException ex) {}

//            if ((temp == null) || (temp.isEmpty())) {
//                finalPlMap.put(anotherMz, plMap.get(mz));
//            }
        }

        // Add two virtual peak. Because we have convert all y-ions to b-ions.
        finalPlMap.put(MassTool.PROTON, 1d);
        double cTermMz = precursorMass - massTool.H2O + MassTool.PROTON;
        double leftMz = cTermMz - ms2Tolerance;
        double rightMz = cTermMz + ms2Tolerance;
        NavigableMap<Double, Double> temp = null;
        try {
            temp = plMap.subMap(leftMz, true, rightMz, true);
        } catch (IllegalArgumentException ex) {}
        if ((temp == null) || (temp.isEmpty())) {
            finalPlMap.put(cTermMz, 1d);
        }
        finalPlMap.put(MassTool.PROTON + massTool.H2O, 1d);


        //shouldDoMyDeisotope
        return finalPlMap;
    }
}
