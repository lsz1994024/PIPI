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

import ProteomicsLibrary.MassTool;
import ProteomicsLibrary.Score;
import ProteomicsLibrary.SpecProcessor;
import ProteomicsLibrary.Types.SparseBooleanVector;
import ProteomicsLibrary.Types.SparseVector;
import com.google.common.collect.Sets;
import gurobi.GRB;
import gurobi.GRBEnv;
import org.apache.commons.math3.util.Pair;
import proteomics.FM.FMRes;
import proteomics.Index.BuildIndex;
import proteomics.PTM.InferPTM;
import proteomics.Segment.InferSegment;
import proteomics.Spectrum.DatasetReader;
import proteomics.Types.*;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;

import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.locks.ReentrantLock;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

import static proteomics.PIPI.*;
import static proteomics.PTM.InferPTM.*;
import static proteomics.Segment.InferSegment.*;

public final class MainSearch implements Callable<MainSearch.Entry> {
    private static final int candisNum = 20;
    private final BuildIndex buildIndex;
    private final MassTool massTool;
    private final double ms1Tolerance;
    private final double minPtmMass;
    private final double maxPtmMass;
    private final JMzReader spectraParser;
    private final ReentrantLock lock;
    private final String scanName;
    private final int precursorCharge;
    private final double precursorMass;
    private final SpecProcessor specProcessor;
    private final int scanNum;
    //    private String truth;
    private final double ms2Tolerance;
    private final InferPTM inferPTM;
    //    private final int minPepLen;
//    private final int maxPepLen;
    public MainSearch(int scanNum, BuildIndex buildIndex, MassTool massTool, double ms2Tolerance, double ms1Tolerance, double minPtmMass, double maxPtmMass
            , JMzReader spectraParser, ReentrantLock lock, String scanName, int precursorCharge, double precursorMass
            , SpecProcessor specProcessor) {

        this.buildIndex = buildIndex;
        this.massTool = massTool;
        this.ms1Tolerance = ms1Tolerance;
        this.minPtmMass = minPtmMass;
        this.maxPtmMass = maxPtmMass;
//        this.localMaxMs2Charge = localMaxMs2Charge;
        this.spectraParser = spectraParser;
        this.lock = lock;
        this.scanName = scanName;
        this.precursorCharge = precursorCharge;
        this.precursorMass = precursorMass;
        this.specProcessor = specProcessor;
        this.scanNum = scanNum;
        this.ms2Tolerance = ms2Tolerance;
        this.inferPTM = buildIndex.getInferPTM();
    }

    @Override
    public Entry call() throws Exception {
        Map<Double, Double> rawPLMap;
        try {
            lock.lock();
            rawPLMap = spectraParser.getSpectrumById(scanName.split("\\.")[1]).getPeakList();
//            spectraParser.
        } finally {
            lock.unlock();
        }

        double ms1TolAbs = Double.parseDouble(InferPTM.df3.format(precursorMass*ms1Tolerance/1000000));

        TreeMap<Double, Double> plMap = specProcessor.preSpectrumTopNStyleWithChargeLimit(rawPLMap, precursorMass, precursorCharge, DatasetReader.topN);

        if (plMap.isEmpty()) return null;
        // Coding
        InferSegment inferSegment = buildIndex.getInferSegment();
        TreeMap<Double, Double> finalPlMap = inferSegment.addVirtualPeaks(precursorMass, plMap);
        SparseVector expProcessedPL;
        if (PIPI.useXcorr) {
            expProcessedPL = specProcessor.prepareXCorr(finalPlMap, false);
        } else {
            expProcessedPL = specProcessor.digitizePL(finalPlMap);
        }
//        List<ExpTag> allLongTagList = inferSegment.getLongTag(finalPlMap, precursorMass - massTool.H2O + MassTool.PROTON, scanNum, minTagLenToExtract,maxTagLenToExtract);
        List<ExpTag> allLongTagList = inferSegment.getLongTag(finalPlMap,  minTagLenToExtract);

//        if (lszDebugScanNum.contains(scanNum)) {
//            int a = 1;
//            System.out.println("tagNum, "+ allLongTagList.size() +"," +finalPlMap.size());
//        }

        List<ExpTag> cleanedAllLongTagList = inferSegment.cleanAbundantTagsPrefix(allLongTagList, minTagLenToExtract);
        allLongTagList = cleanedAllLongTagList;
        if (allLongTagList.isEmpty())  return null;

        double totalMass = precursorMass + 2 * MassTool.PROTON;


        Map<String, PeptideInfo> peptideInfoMap = new HashMap<>(50000);

        int minTagLen = 3;

        GRBEnv env = new GRBEnv(true);
        if (java.lang.management.ManagementFactory.getRuntimeMXBean().getInputArguments().toString().indexOf("jdwp") < 0){
            env.set(GRB.IntParam.OutputFlag,0);
            env.set(GRB.IntParam.LogToConsole, 0);
        }
        env.start();


        int tagNumsToTest =  Math.min(20, allLongTagList.size());
        List<ExpTag> tagsToTest = new ArrayList<>();
        TreeSet<Peptide> resPeptTreeSet = new TreeSet<>(Collections.reverseOrder());
        ExpTag tagInfo;
        for (int i = 0; i < tagNumsToTest; i++) {
            tagInfo = allLongTagList.get(i);
            if (tagInfo.isNorC == N_TAG) { //n tag //todo can I relax this. can I check all N and C tags as non NC (both direction)
                tagsToTest.add(tagInfo);
            } else if (tagInfo.isNorC == C_TAG) { // c tag
                ExpTag revTagInfo = tagInfo.revTag(totalMass);
                tagsToTest.add(revTagInfo);
            } else { // non-nc tag
                tagsToTest.add(tagInfo);
                ExpTag revTagInfo = tagInfo.revTag(totalMass);
                tagsToTest.add(revTagInfo);
            }
        }
        List<OccGroup> occGroupList = new LinkedList<>();
        searchWithMultiTags(scanNum , tagsToTest, occGroupList, minTagLen, totalMass);
        occGroupList.sort(Comparator.comparing(o->o.totalScore, Comparator.reverseOrder()));

//        mergeTags(occGroupList);

        int count = 0;
//        double topScore = 0;
        if (occGroupList.isEmpty()) {
//            System.out.println(scanNum + ", occGroupList.isEmpty()," + occGroupList.size());
            return null;
        }
        double topOccScore = occGroupList.get(0).totalScore;
//        outLoop:
        for (OccGroup occG : occGroupList) {
            if (occG.totalScore < 0.7*topOccScore || count > 20) break;
            occG.tagPosList.sort(Comparator.comparing(o->o.getSecond()));
            for (Pair<ExpTag, Integer> tagPosPair : occG.tagPosList) {
                ExpTag tmpTag = tagPosPair.getFirst();
                if (Math.abs(tmpTag.getHeadLocation()-MassTool.PROTON) < ms2Tolerance) {
                    tmpTag.isNorC = N_TAG;
                } else if (Math.abs(tmpTag.getTailLocation()-precursorMass+ massTool.H2O-MassTool.PROTON) < ms2Tolerance) {
                    tmpTag.isNorC = C_TAG;
                }else {
                    tmpTag.isNorC = NON_NC_TAG;
                }
            }
//            occG.tagPosList.remove(2);// debug
//            occG.tagPosList.remove(1);// debug
            String protId = occG.protId;
            count++;
            addCandisWithMultiMaxPeakMILP(scanNum, protId, occG.tagPosList, ms1TolAbs, resPeptTreeSet, peptideInfoMap, expProcessedPL, plMap, env);
        }
        if (lszDebugScanNum.contains(scanNum)) {
            int a = 1;
        }
        if (resPeptTreeSet.isEmpty()) {
            env.dispose();
            return null;
        }
        List<Peptide> pepList = new ArrayList<>(resPeptTreeSet);
        Set<Integer> pepIdsToRemove = new HashSet<>();
        for (int j = pepList.size()-1; j >= 1; j--) {
            if (pepIdsToRemove.contains(j)) continue;
            for (int i = 0; i < j; i++) {
                if (pepIdsToRemove.contains(i)) continue;
                try {
                    if (isHomo(pepList.get(i), pepList.get(j), peptideInfoMap) || pepList.get(j).getScore() > 0.6*pepList.get(i).getScore()) {
                        int iPriority = pepList.get(i).getPriority();
                        int jPriority = pepList.get(j).getPriority();
                        if (iPriority < jPriority) {
                            pepIdsToRemove.add(i);
                        } else if (iPriority > jPriority) {
                            pepIdsToRemove.add(j);
                        } else {// iPriority == jPriority
                            if (onlyDifferUnsettledPtm(pepList.get(i), pepList.get(j))) {
                                pepIdsToRemove.add(j);
                            }
                        }
                    }
                } catch (Exception e) {
                    System.out.println(scanNum + ", isHomo");
                }
            }
        }
        List<Peptide> newPepList = new ArrayList<>();
        for (int id = 0; id < pepList.size(); id++){
            if (pepIdsToRemove.contains(id)) continue;
            newPepList.add(pepList.get(id));
        }
        if (newPepList.isEmpty()) {
            env.dispose();
            return null;
        }

        Peptide[] peptideArray = newPepList.toArray(new Peptide[0]);
        Peptide topPep = peptideArray[0];

        Entry entry;
        String pepSetString = "";
        for (Peptide peptide : peptideArray){
            PeptideInfo peptideInfo = peptideInfoMap.get(peptide.getFreeSeq());
            pepSetString += peptide.getVarPtmContainingSeqNow() + "," + peptide.getScore() + "," + String.join("_", peptideInfo.protIdSet) +",";
        }

        double deltaLCn = 1; // L means the last?
        if (peptideArray.length > candisNum - 1) {
            deltaLCn = (peptideArray[0].getScore() - peptideArray[candisNum - 1].getScore()) / peptideArray[0].getScore();
        }
        double deltaCn = 1;
        if (peptideArray.length > 1) {
            for(int i = 0; i < peptideArray.length; i++) {
                if (peptideArray[i].getScore() != peptideArray[0].getScore()){
                    deltaCn = (peptideArray[0].getScore() - peptideArray[i].getScore()) / peptideArray[0].getScore();
                    break;
                }
            }
        }
        String otherPtmPatterns = "-";
        entry = new MainSearch.Entry(
                scanNum, scanName, precursorCharge, precursorMass, buildIndex.getLabelling(), topPep.getPtmContainingSeq(buildIndex.returnFixModMap())
                , topPep.getTheoMass(), topPep.isDecoy() ? 1 : 0, topPep.getScore(), deltaLCn, deltaCn
                , topPep.getMatchedPeakNum(), topPep.getIonFrac(), topPep.getMatchedHighestIntensityFrac()
                , topPep.getExplainedAaFrac(), otherPtmPatterns, topPep.getaScore()
                , pepSetString.substring(0, pepSetString.length()-1)
        );
        entry.topPeptide = topPep;
        entry.topPeptide.precursorMass = precursorMass;
        entry.varPtmList.addAll(topPep.posVarPtmResMap.values());
        for (Peptide peptide : peptideArray) {
            entry.peptideInfoMapForRef.put(peptide.getFreeSeq(), peptideInfoMap.get(peptide.getFreeSeq()));
        }
//        System.out.println(scanNum + ","+entry.peptideInfoMapForRef.size() + "," + peptideArray.length);l
//        if (lszDebugScanNum.contains(scanNum)){
//            int a = 1;
//            System.out.println(scanNum + "," + pepSetString);
//        }
        int c = 1;
//        env.release();l
        env.dispose();
        return entry;
    }

    private boolean mergeTagPosList(List<Pair<ExpTag, Integer>> tagRelPosList, String protSeq) {
//        if (tagRelPosList.size() == 1) return;
//        Pair<ExpTag, Integer> newTagPosPair = tagRelPosList.get(tagRelPosList.size()-1);
        boolean successAdded = true;
        tagRelPosList.sort(Comparator.comparing(o->o.getSecond()));
        int lastTagId = 0;
        int contTagId = 1;
        List<Integer> tagIdAppended = new ArrayList<>();
//        List<Integer> badIds = new ArrayList<>();
        while (contTagId < tagRelPosList.size()) {
            int lastTagPos = tagRelPosList.get(lastTagId).getSecond();
            ExpTag lastTag = tagRelPosList.get(lastTagId).getFirst();
            int contTagPos = tagRelPosList.get(contTagId).getSecond();
            ExpTag contTag = tagRelPosList.get(contTagId).getFirst();

            if (contTagPos <= lastTagPos+lastTag.size()) {
                double lastFlagMass;
                double contFlagMass = contTag.getHeadLocation();
                if (lastTagPos == contTagPos) {
                    lastFlagMass = lastTag.getHeadLocation();
                } else {
                    lastFlagMass = lastTag.expAaList.get(contTagPos-lastTagPos-1).getTailLocation();
                }
                if (Math.abs(contFlagMass-lastFlagMass) < ms2Tolerance) {
                    lastTag.appendTag(contTagPos-lastTagPos, contTag);
                    tagIdAppended.add(contTagId);
                    contTagId++;
                } else {
                    int badTadId = -1;
                    if (lastTag.isNorC == N_TAG) {
                        badTadId = contTagId;
                    } else if (contTag.isNorC == C_TAG) {
                        badTadId = lastTagId;
                    } else {
                        badTadId = lastTag.getTotalIntensity() < contTag.getTotalIntensity() ? lastTagId : contTagId;
                    }
                    tagIdAppended.add(badTadId);
//                    badIds.add(badTadId);
                    if (badTadId == lastTagId) {
                        lastTagId = contTagId;
                        contTagId = lastTagId+1;
                    } else {
                        contTagId++;
                    }
                    successAdded = false;
                }
            } else {
                int badTadId = -1;
                if (lastTag.isNorC == N_TAG) {
                    badTadId = contTagId;
                } else if (contTag.isNorC == C_TAG) {
                    badTadId = lastTagId;
                } else {
                    badTadId = lastTag.getTotalIntensity() < contTag.getTotalIntensity() ? lastTagId : contTagId;
                }

                if ( isValidGap(lastTag, lastTagPos, contTag, contTagPos, protSeq))
                {
                    lastTagId = contTagId;
                    contTagId = lastTagId+1;
                } else  { // if they dont overlap but the mass diff is wrong, remove one
                    tagIdAppended.add(badTadId);
                    if (badTadId == lastTagId) {
                        lastTagId = contTagId;
                        contTagId = lastTagId+1;
                    } else {
                        contTagId++;
                    }
                    successAdded = false;
                }
            }
        }
        tagIdAppended.sort(Comparator.reverseOrder());
        for (int i : tagIdAppended) {
            tagRelPosList.remove(i);
        }
        Collections.sort(tagRelPosList, Comparator.comparing(o -> o.getSecond()));
        return successAdded;
    }

    private boolean isValidGap(ExpTag lastTag, int lastTagPos, ExpTag contTag, int contTagPos, String protSeq){
        double deltaMass = contTag.getHeadLocation() - lastTag.getTailLocation();
        double sumAaMass = massTool.calResidueMass(protSeq.substring(lastTagPos+lastTag.size(), contTagPos));
        if (deltaMass < sumAaMass + 250 && deltaMass > Math.max( 57*(contTagPos-lastTagPos-lastTag.size()), sumAaMass - 250)) {
            return true;
        }
        return false;
    }
    final class OccGroup implements Cloneable {
        public List<Pair<ExpTag, Integer>> tagPosList = new ArrayList<>();
        public int lPos;
        public int rPos;
        public int spanLen;
        public String protId;
        public double totalScore;
        private int hashcode;
        public OccGroup(ExpTag initialTag, String protId, int tagPosInProt) {
            this.tagPosList.add(new Pair<>(initialTag, tagPosInProt));
            this.spanLen = initialTag.size();
            this.lPos = tagPosInProt;
            this.rPos = tagPosInProt + this.spanLen;
            this.protId = protId;
            this.totalScore = initialTag.getTotalIntensity();

            StringBuilder sb = new StringBuilder();
            sb.append(protId);
            sb.append(lPos);
            sb.append(rPos);
            sb.append(totalScore);
            for (int i = 0; i < this.tagPosList.size(); ++i) {
                sb.append(tagPosList.get(i).getFirst().getFreeAaString());
                sb.append(tagPosList.get(i).getFirst().getHeadLocation());
            }
            hashcode = sb.toString().hashCode();
        }

        public OccGroup() {} // for clone


        public boolean addTag(ExpTag newTag, int newTagPos) {
            boolean shouldAdd = true;
            for (Pair<ExpTag, Integer> tagPosPair : this.tagPosList) {
                if (tagPosPair.getFirst().getFreeAaString().contentEquals(newTag.getFreeAaString())
                        && tagPosPair.getSecond() == newTagPos
                        && tagPosPair.getFirst().getHeadLocation() == newTag.getHeadLocation()) { //dont add exactly the same tag
                    shouldAdd = false;
                    break;
                }
            }
            if (! shouldAdd) return true;
            this.tagPosList.add(new Pair<>(newTag, newTagPos));
            boolean successAdded = mergeTagPosList(this.tagPosList, buildIndex.protSeqMap.get(protId));

            //update variables
            this.totalScore = 0;
            for (Pair<ExpTag, Integer> tagPosPair : this.tagPosList) this.totalScore += tagPosPair.getFirst().getTotalIntensity();
            this.lPos = this.tagPosList.get(0).getSecond();
            this.rPos = this.tagPosList.get(this.tagPosList.size()-1).getSecond() + this.tagPosList.get(this.tagPosList.size()-1).getFirst().size();
            this.spanLen = this.rPos - this.lPos;

            StringBuilder sb = new StringBuilder();
            sb.append(protId);
            sb.append(lPos);
            sb.append(rPos);
            sb.append(totalScore);
            for (int i = 0; i < this.tagPosList.size(); ++i) {
                sb.append(tagPosList.get(i).getFirst().getFreeAaString());
                sb.append(tagPosList.get(i).getFirst().getHeadLocation());
            }
            hashcode = sb.toString().hashCode(); //update hashcode

            return successAdded;
        }

        public OccGroup clone() {
            OccGroup newOccG = new OccGroup();
            for (Pair<ExpTag, Integer> tagPosPair : this.tagPosList) {
                newOccG.tagPosList.add(new Pair<>(tagPosPair.getFirst().clone(), tagPosPair.getSecond()));
            }
//            newOccG.tagPosList.addAll(this.tagPosList);
            newOccG.spanLen = this.spanLen;
            newOccG.lPos = this.lPos;
            newOccG.rPos = this.rPos;
            newOccG.protId = this.protId;
            newOccG.totalScore = this.totalScore;
            newOccG.hashcode = this.hashcode;
            return newOccG;
        }

        @Override
        public int hashCode() {
            return hashcode;
        }
        @Override
        public boolean equals(Object other) {
            if (this == other)
                return true;
            if (other == null || getClass() != other.getClass())
                return false;
            OccGroup tmp = (OccGroup) other;
            return hashCode() == tmp.hashCode();
        }

    }
    private void searchWithMultiTags(int scanNum, List<ExpTag> tagsToTest, List<OccGroup> occGroupList, int minTagLen, double totalMass) {

        char[] tagChar;
        FMRes fmRes;
        int matchedPos, tagLen;
        String protId;
        ExpTag tagInfo;
        int tagsNum = tagsToTest.size();
//        Set<Pair<String, Integer>> matchedProtIdRelPosSet = new HashSet<>();
        int occId = -1;
        for (int i = 0; i < tagsNum; i++) {
            boolean shouldTryForwardSearch = false;
            tagInfo = tagsToTest.get(i);
            tagLen = tagInfo.size();
            tagChar = tagInfo.getFreeAaString().toCharArray();
            fmRes = buildIndex.fmIndexNormal.fmSearchFuzzy(tagChar);
            if (fmRes != null) {
                matchedPos = fmRes.matchedPos;
//                int oriIsNorC = tagInfo.isNorC;
                if (matchedPos != 0) {
                    shouldTryForwardSearch = true;
//                    if (oriIsNorC == N_TAG) tagInfo.isNorC = NON_NC_TAG;
                }
                tagLen = tagInfo.size();
                if (matchedPos == 0 || tagLen-matchedPos >= minTagLen) {
                    ExpTag resTag = matchedPos == 0 ? tagInfo : tagInfo.subTag(matchedPos, tagLen);
                    int absTagPos, dotIndex, lPos, rPos;

                    Set<OccGroup> oldOccGToDel = new HashSet<>();
                    Set<OccGroup> newOccGToAdd = new HashSet<>();
                    for (int ii = fmRes.sp; ii <= fmRes.ep; ii++) {
                        absTagPos = buildIndex.fmIndexNormal.SA[ii];
                        dotIndex = -2-Arrays.binarySearch(buildIndex.dotPosArrNormal, absTagPos);
                        protId = buildIndex.posProtMapNormal.get(dotIndex);
                        lPos = absTagPos - buildIndex.dotPosArrNormal[dotIndex] - 1;
                        rPos = lPos + tagLen - matchedPos;
                        boolean belongToExistingOcc = false;

                        for (OccGroup occG : occGroupList) {
                            if ( occG.protId.contentEquals(protId) && (lPos <= occG.rPos + occG.spanLen) && (rPos >= occG.lPos - occG.spanLen) ) {
                                OccGroup newOccG = occG.clone();
                                boolean successAdded = newOccG.addTag(resTag.clone(), lPos);
                                if ( successAdded ) { //successfully added
                                    oldOccGToDel.add(occG);
                                    newOccGToAdd.add(newOccG);
                                    belongToExistingOcc = true;
                                }
                            }
                        }
                        if ( ! belongToExistingOcc) {  //new occurrence that does not belong to any pos range
                            OccGroup occG = new OccGroup(resTag.clone(), protId, lPos);
                            newOccGToAdd.add(occG);
                        }
                    }
                    for (OccGroup oldOccG : oldOccGToDel)  occGroupList.remove(oldOccG);
                    for (OccGroup newOccG : newOccGToAdd) {
                        occGroupList.add(newOccG);
                    }
                }
            }

            if (shouldTryForwardSearch) {
                int absTagPos;
                int dotIndex;
                int lPos;
                int rPos;

                ExpTag revTagInfo = tagInfo.revTag(totalMass);
                tagChar = revTagInfo.getFreeAaString().toCharArray();
                fmRes = buildIndex.fmIndexReverse.fmSearchFuzzy(tagChar);
                if (fmRes != null) {
                    matchedPos = fmRes.matchedPos;
                    if (matchedPos > 0 && tagInfo.size() - matchedPos < minTagLen) continue;

                    ExpTag resRevTag = matchedPos == 0 ? tagInfo : tagInfo.subTag(0, tagLen - matchedPos);
//                    Map<String, Set<Pair<Integer, Integer>>> protId_posRangesToDel = new HashMap<>();
                    Set<OccGroup> oldOccGToDel = new HashSet<>();
                    Set<OccGroup> newOccGToAdd = new HashSet<>();
                    for (int ii = fmRes.sp; ii <= fmRes.ep; ii++) {
                        absTagPos = buildIndex.fmIndexReverse.SA[ii];
                        absTagPos = buildIndex.textNormalLen - absTagPos - tagLen + matchedPos; // should I minus one more?
                        dotIndex = -2 - Arrays.binarySearch(buildIndex.dotPosArrNormal, absTagPos);
                        protId = buildIndex.posProtMapNormal.get(dotIndex);
                        lPos = absTagPos - buildIndex.dotPosArrNormal[dotIndex] - 1;
                        rPos = lPos + tagLen - matchedPos;

                        boolean foundRange = false;
                        for (OccGroup occG : occGroupList) {
                            if ( occG.protId.contentEquals(protId) && (lPos <= occG.rPos + occG.spanLen) && (rPos >= occG.lPos - occG.spanLen) ) {
                                OccGroup newOccG = occG.clone();
                                boolean successAdded = newOccG.addTag(resRevTag.clone(), lPos); // maybe in this func  just use merge tag
                                if (successAdded) {
                                    oldOccGToDel.add(occG);
                                    newOccGToAdd.add(newOccG);
                                    foundRange = true;
                                }
                            }
                        }
                        if ( ! foundRange) {  //new occurrence that does not belong to any pos range
                            OccGroup occG = new OccGroup(resRevTag.clone(), protId, lPos);
                            newOccGToAdd.add(occG);
                        }
                    }
                    for (OccGroup oldOccG : oldOccGToDel) occGroupList.remove(oldOccG);
                    for (OccGroup newOccG : newOccGToAdd) {
                        occGroupList.add(newOccG);
                    }
                }
            }
        }

    }
    private boolean isHomo(Peptide p1, Peptide p2, Map<String, PeptideInfo> peptideInfoMap) {
        Set<String> temp1 = peptideInfoMap.get(p1.getFreeSeq()).protIdSet;
        Set<String> temp2 = peptideInfoMap.get(p2.getFreeSeq()).protIdSet;

        Set<String> set = new HashSet<>(temp1);
        set.retainAll(temp2);
        SparseBooleanVector sbv1 = buildIndex.inferSegment.generateSegmentBooleanVector(p1.getFreeSeq());
        SparseBooleanVector sbv2 = buildIndex.inferSegment.generateSegmentBooleanVector(p2.getFreeSeq());
        return sbv1.dot(sbv2) > 0.3*Math.min(p1.getFreeSeq().length(), p2.getFreeSeq().length());
    }

    private boolean onlyDifferUnsettledPtm(Peptide p1, Peptide p2) {
        SparseBooleanVector sbv1 = buildIndex.inferSegment.generateSegmentBooleanVector(p1.getFreeSeq());
        SparseBooleanVector sbv2 = buildIndex.inferSegment.generateSegmentBooleanVector(p2.getFreeSeq());
        if (sbv1.dot(sbv2) < 0.3*Math.min(p1.getFreeSeq().length(), p2.getFreeSeq().length())){
            return false;
        }

        if (p1.posVarPtmResMap.size() != p2.posVarPtmResMap.size()) {
            return false;
        }
        if (!p1.posVarPtmResMap.keySet().containsAll(p2.posVarPtmResMap.keySet())){
            return false;
        }
        byte n_SameMass = 0;
        for (int pos : p2.posVarPtmResMap.keySet()) {
            if (p1.posVarPtmResMap.get(pos).mass == p2.posVarPtmResMap.get(pos).mass) {
                n_SameMass++;
            } else if (Math.abs(p1.posVarPtmResMap.get(pos).mass - p2.posVarPtmResMap.get(pos).mass) < 0.02) {
                if (p1.posVarPtmResMap.get(pos).classification.contains("PIPI_unsettled") || p2.posVarPtmResMap.get(pos).classification.contains("PIPI_unsettled")) {
                    n_SameMass++;
                }
            }
//            if (Math.abs(p1.posVarPtmResMap.get(pos).mass - p2.posVarPtmResMap.get(pos).mass) < 0.02
//                    && (p1.posVarPtmResMap.get(pos).classification.contains("PIPI_unsettled") || p2.posVarPtmResMap.get(pos).classification.contains("PIPI_unsettled"))
//            ) {
//                n_SameMass++;
//            }
        }
        return n_SameMass == p1.posVarPtmResMap.size();
    }

    private double addCandisWithMultiMaxPeakMILP(int scanNum, String protId, List<Pair<ExpTag, Integer>> tagRelPosList, double ms1TolAbs, TreeSet<Peptide> resPepTreeSet
            , Map<String, PeptideInfo> peptideInfoMap, SparseVector expProcessedPL, TreeMap<Double, Double> plMap, GRBEnv env) {

        int tagNum;
        String protSeq = buildIndex.protSeqMap.get(protId);
        int protLen = protSeq.length();
        if (tagRelPosList.get(0).getSecond() < 0
            || tagRelPosList.get(0).getSecond() == 0 && tagRelPosList.get(0).getFirst().isNorC != N_TAG) {
            tagRelPosList.remove(0);
        }
        tagNum = tagRelPosList.size();
        if (tagNum == 0) return 0;
        if (tagRelPosList.get(tagNum-1).getSecond() >= protLen
            || tagRelPosList.get(tagNum-1).getSecond() + tagRelPosList.get(tagNum-1).getFirst().size() == protLen && tagRelPosList.get(tagNum-1).getFirst().isNorC != C_TAG) {
            tagRelPosList.remove(tagNum-1);
        }
        tagNum = tagRelPosList.size();
        if (tagNum == 0) return 0;

        double maxScore = 0;
//        TreeMap<Double, Double> unUsedPlMap = new TreeMap<>(plMap);
//        for (Pair<ExpTag, Integer> tagRelPos : tagRelPosList) {
//            for (ExpAa expAa : tagRelPos.getFirst().expAaList) {
//                if (plMap.containsKey(expAa.getHeadLocation())) {
//                    unUsedPlMap.remove(expAa.getHeadLocation());
//                }
//                if (plMap.containsKey(expAa.getTailLocation())) {
//                    unUsedPlMap.remove(expAa.getTailLocation());
//                }
//            }
//        }
//        SparseVector unUsedExpPeakVec;
//        unUsedExpPeakVec = specProcessor.digitizePL(unUsedPlMap);

        Pair<ExpTag, Integer> lastTagPosPair = tagRelPosList.get(tagNum-1);
        int remainMC = massTool.missedCleavage -
                getNumOfMCFromStr(protSeq.substring(tagRelPosList.get(0).getSecond(), lastTagPosPair.getSecond()+ lastTagPosPair.getFirst().size()));
        if (lastTagPosPair.getFirst().isNorC == C_TAG && isKR(lastTagPosPair.getFirst().getFreeAaString().charAt(lastTagPosPair.getFirst().size()-1))) {
            remainMC++; // the last KR at pepC does not count as MC
        }

        List<TreeSet<Peptide>> segResList = new ArrayList<>();
        for (int gId = 0; gId < tagNum+1; ++gId) {
            ExpTag lTag = null;
            int lTagPos = 0;
            if (gId != 0) {
                lTag = tagRelPosList.get(gId-1).getFirst();
                lTagPos = tagRelPosList.get(gId-1).getSecond();
            }
            ExpTag rTag = null;
            int rTagPos = 0;
            if (gId != tagNum) {
                rTag = tagRelPosList.get(gId).getFirst();
                rTagPos = tagRelPosList.get(gId).getSecond();
            }
            boolean solved;
            if (gId == 0) {
                if (Math.abs(rTag.getHeadLocation()-MassTool.PROTON) < ms2Tolerance) {//rTag.isNorC == N_TAG
                    solved = true;
                } else {
                    TreeSet<Peptide> nModPepsSet = new TreeSet<>(Comparator.reverseOrder());
                    solved = solveGapN(scanNum, rTag, rTagPos, protSeq, remainMC, ms1TolAbs, expProcessedPL, plMap, env, nModPepsSet);
                    if (solved) segResList.add(nModPepsSet);
                }
                if (solved) segResList.add(getPepFromTag(rTag));
            } else if (gId == tagNum) {// solve a C gap
                if (Math.abs(lTag.getTailLocation()-precursorMass+ massTool.H2O-MassTool.PROTON) < ms2Tolerance) {// lTag.isNorC == C_TAG
                    solved = true;
                } else {
                    TreeSet<Peptide> cModPepsSet = new TreeSet<>(Comparator.reverseOrder());
                    solved = solveGapC(scanNum, lTag, lTagPos, protId, protSeq, remainMC,
                            ms1TolAbs, expProcessedPL, plMap, env, cModPepsSet);
                    if (solved) segResList.add(cModPepsSet);
                }
            } else {// solve a middle gap
                TreeSet<Peptide> midModPepsSet = new TreeSet<>(Comparator.reverseOrder());
//                try {
//                    protSeq.substring(tagRelPosList.get(gId-1).getSecond()+tagRelPosList.get(gId-1).getFirst().size(), tagRelPosList.get(gId).getSecond());
//                } catch (Exception e ) {
//                    System.out.println(scanNum + ", string index error, " + protId + "," + tagPosId);
//                }
                String gapSeq = protSeq.substring(tagRelPosList.get(gId-1).getSecond()+tagRelPosList.get(gId-1).getFirst().size(), tagRelPosList.get(gId).getSecond());
                double deltaMass = tagRelPosList.get(gId).getFirst().getHeadLocation() - tagRelPosList.get(gId-1).getFirst().getTailLocation() - massTool.calResidueMass(gapSeq);
                solved = solveGapMid(scanNum, lTag, lTagPos, rTag, rTagPos, gapSeq, deltaMass, ms1TolAbs, expProcessedPL, plMap, env, midModPepsSet);
                if (solved) {
                    segResList.add(midModPepsSet);
                    segResList.add(getPepFromTag(rTag));
                }
            }
            if (!solved) { // todo makeup solution
//                System.out.println(scanNum + ", falied not solved");
                return 0;
            }
        }
        //merge results
//        try {
        if ( ! segResList.isEmpty()) {
            collectResult(segResList, resPepTreeSet, expProcessedPL, plMap,  peptideInfoMap, protId, protSeq, tagRelPosList.get(0).getSecond(), tagRelPosList.get(0).getFirst().isNorC==N_TAG, tagRelPosList.get(tagNum-1).getSecond()+tagRelPosList.get(tagNum-1).getFirst().size(), tagRelPosList.get(tagNum-1).getFirst().isNorC==C_TAG);
        }
//        } catch (Exception e ) {
//            System.out.println(scanNum +",collectResult");
//        }
        return maxScore;
    }

    private void collectResult(List<TreeSet<Peptide>> segResList, TreeSet<Peptide> resPepTreeSet, SparseVector expProcessedPL, TreeMap<Double,Double> plMap,
                               Map<String, PeptideInfo> peptideInfoMap, String protId, String protSeq, int firstTagStartPos, boolean nIsTag, int lastTagEndPos, boolean cIsTag) {
        List<Set<Integer>> idSetList = new ArrayList<>(segResList.size());
        for (int segNum = 0; segNum < segResList.size(); segNum++) {
            Set<Integer> idSet = new HashSet<>();
            for (int id = 0; id < segResList.get(segNum).size(); ++id) {
                idSet.add(id);
            }
            idSetList.add(idSet);
        }

        for (List<Integer> resIdList : Sets.cartesianProduct(idSetList)) {
            StringBuilder pepSeqSB = new StringBuilder();
            PosMassMap posMassMap = new PosMassMap();
            TreeMap<Integer, VarPtm> posVarPtmResMap = new TreeMap<>();

            int curRelPos = 0;
            for (int i = 0; i < resIdList.size(); ++i) {
                int id = resIdList.get(i);
                Peptide[] segArr = new Peptide[segResList.get(i).size()];
                segResList.get(i).toArray(segArr);
                Peptide curPepSeg = segArr[id];
                pepSeqSB.append(curPepSeg.getFreeSeq());

                if (! curPepSeg.posVarPtmResMap.isEmpty()) {
                    for (int j : curPepSeg.posVarPtmResMap.keySet()) {
                        posMassMap.put(j+curRelPos, curPepSeg.posVarPtmResMap.get(j).mass);
                        posVarPtmResMap.put(j+curRelPos, curPepSeg.posVarPtmResMap.get(j));
                    }
                }
                curRelPos += segArr[id].getFreeSeq().length();
            }
            Peptide resPep = new Peptide(pepSeqSB.toString(), false, massTool);
            if (! posMassMap.isEmpty()) {
                resPep.setVarPTM(posMassMap);
                resPep.posVarPtmResMap.putAll(posVarPtmResMap);
            }
            double calScore = massTool.buildVectorAndCalXCorr(resPep.getIonMatrixNow(), 1, expProcessedPL, resPep.matchedBions, resPep.matchedYions);//todo decide the penalty

            resPep.setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, resPep.getIonMatrixNow(), ms2Tolerance));
            resPep.setScore(calScore*(1 - resPep.posVarPtmResMap.size() * 0.02));

            updatePeptideTreeSet(resPep, resPepTreeSet, peptideInfoMap, protId, protSeq, protSeq.indexOf(pepSeqSB.toString()), protSeq.indexOf(pepSeqSB.toString())+pepSeqSB.length()-1);
        }
    }

    private TreeSet<Peptide> getPepFromTag(ExpTag tag) {
        String freeSeq = tag.getFreeAaString();
        String modSeq  = tag.getPtmAaString();
        Peptide peptide = new Peptide(freeSeq, false, massTool);
        TreeSet<Peptide> modPepPool = new TreeSet<>();
        modPepPool.add(peptide);
        if (! freeSeq.contentEquals(modSeq)) {
            PosMassMap posMassMap = new PosMassMap();
            int idOfAa = 0;
            for (char aaChar : tag.getPtmAaString().toCharArray()) {
                if (Character.isUpperCase(aaChar)) {
                    idOfAa += 1;
                } else {
                    posMassMap.put(idOfAa, massTool.labelVarPtmMap.get(aaChar).mass);
                    peptide.posVarPtmResMap.put(idOfAa, massTool.labelVarPtmMap.get(aaChar));
                }
            }
        }
        return modPepPool;
    }
    private boolean solveGapC(int scanNum, ExpTag tag, int tagPosInProt, String protId, String protSeq, int remainMC, double ms1TolAbs,
                              SparseVector unUsedExpProcessedPL, TreeMap<Double, Double> unUsedPlMap, GRBEnv env, TreeSet<Peptide> cModPepsSet) {

        double tagCMass = tag.getTailLocation() + massTool.H2O - MassTool.PROTON; // +massTool.H2O-MassTool.PROTON  is to mimic the mass of the real neutral precursor mass
        if (tag.isNorC == C_TAG || Math.abs(tagCMass - precursorMass) < ms1TolAbs) return true;// C tag, tagCMass must be exactly settle. Otherwise no valid cPos will be, then no valid n pos.
        // only when oriTag is not C oriTag can it extend to c
        double cCutMass = tag.getTailLocation() - MassTool.PROTON;
        int tagLen = tag.size();
        int protLen = protSeq.length();
        List<Integer> cPoses = new ArrayList<>();
        List<Integer> krPoses = new ArrayList<>();
        int mcInC = 0;
        boolean startRecord = false;
        char aaChar;
        boolean isKR;

        for (int cPos = tagPosInProt+tagLen; cPos < protLen; cPos++) {  //
            aaChar = protSeq.charAt(cPos);
            if (isX(aaChar))  break;
            isKR = isKR(aaChar);
            tagCMass += massTool.getMassTable().get(aaChar); //even when tag is fuzzy, tagcmass wont be disturbed
            if (mcInC > remainMC || tagCMass > precursorMass+maxPtmMass) break;
            if (tagCMass >= precursorMass-maxPtmMass) {
                if (startRecord) {
                    cPoses.add(cPos);
                } else {
                    startRecord = true;
                }
                if (isKR)  krPoses.add(cPos);

                if (Math.abs(tagCMass - precursorMass) <= ms1TolAbs) {
                    String cPartSeq = protSeq.substring(tagPosInProt+tagLen, cPos+1);
                    storeCleanPartPeptides(cPartSeq, cCutMass, C_PART, unUsedExpProcessedPL, unUsedPlMap, cModPepsSet);
                    return true;
                }
            }
            if (isKR) mcInC++;
        }// finish collect all possible cposes and krPoses
        if ( !startRecord) return false;// cPoses is empty

        Map<Integer, Integer> posYIdMap = new HashMap<>();
        Map<Integer, Integer> yIdMaxAbsPosMap = new HashMap<>();
        int optStartPos = 1;
        int optEndPosP1 = 0;
        if (buildIndex.cTermSpecific) {
            if ( ! krPoses.isEmpty()) {
                optEndPosP1 = krPoses.get(krPoses.size()-1) + 1;  //max kr
                optStartPos = krPoses.get(0) + 1;  // good trick//min kr
                if ( ! cPoses.isEmpty() && cPoses.get(cPoses.size()-1) == protLen - 1) {
                    optEndPosP1 = protLen;
                }
                int yId = 0;
                for (int pos = optStartPos; pos < optEndPosP1; pos++) {
                    posYIdMap.put(pos, yId); // using rel pos in part seq i.e. start from 0 not the prot pos
                    int oldMaxPos = yIdMaxAbsPosMap.getOrDefault(yId, -99);
                    yIdMaxAbsPosMap.put(yId, pos > oldMaxPos ? pos : oldMaxPos);
                    if (krPoses.contains(pos)) {
                        yId++;
                    }
                }
            } else {
                if (cPoses.isEmpty()) {
                    optEndPosP1 = protLen;
                    optStartPos = protLen;  // good trick
                } else if (cPoses.get(cPoses.size()-1) == protLen-1) {
                    optEndPosP1 = protLen;
                    optStartPos = protLen;  // good trick
                }
            }
        } else {
            optEndPosP1 = cPoses.get(cPoses.size()-1) + 1; //max
            optStartPos = cPoses.get(0);  // good trick  //min
            int yId = 0;
            for (int pos = optStartPos; pos < optEndPosP1; pos++) {
                posYIdMap.put(pos, yId); // using rel pos in part seq i.e. start from 0 not the prot pos
                int oldMaxPos = yIdMaxAbsPosMap.getOrDefault(yId, -99);
                yIdMaxAbsPosMap.put(yId, pos > oldMaxPos ? pos : oldMaxPos);
                yId++;
            }
        }

//        byte isProtNorC_Term = NON_TERM_PROT;
//        if (optStartPos == 1 && protSeq.charAt(0) == 'M') {
//            isProtNorC_Term = N_TERM_PROT;
//        }
        if (optStartPos <= optEndPosP1) {
            int cPartStartPos = tagPosInProt + tagLen;
            String cPartSeq = protSeq.substring(cPartStartPos, optEndPosP1);
            int cPartSeqLen = cPartSeq.length();
            double flexiableMass = precursorMass - (tag.getTailLocation() + massTool.H2O-MassTool.PROTON) - massTool.calResidueMass(protSeq.substring(cPartStartPos, optStartPos));
            Map<Integer, TreeMap<Double, VarPtm>> absPos_MassVarPtm_Map = new HashMap<>(cPartSeqLen, 1);
            Map<Integer, List<Byte>> absPos_ptmPositions_Map = new HashMap<>(cPartSeqLen, 1);
            inferPTM.prepareInfoCTerm(scanNum, cPartSeq, absPos_MassVarPtm_Map,  yIdMaxAbsPosMap, optEndPosP1 == protLen, optStartPos, cPartStartPos, absPos_ptmPositions_Map);

            if (cPartSeqLen == 1) {
                findPtmOnOneAa(cModPepsSet, flexiableMass, absPos_MassVarPtm_Map, cPartSeq, cPartStartPos);
            } else {
                inferPTM.findBestPtmMIPExtCPWL(scanNum, env, flexiableMass, cPartStartPos,
                        cPartSeq, ms1TolAbs, posYIdMap, absPos_MassVarPtm_Map, optStartPos, optEndPosP1 == protLen
                        , yIdMaxAbsPosMap, unUsedPlMap, cModPepsSet, absPos_ptmPositions_Map);
            }
            if (yIdMaxAbsPosMap.isEmpty()) {
                inferPTM.findPosssible1Ptm(scanNum, cPartSeq, absPos_MassVarPtm_Map, cPartStartPos, cModPepsSet, flexiableMass);
            }
        }

        if (cModPepsSet.isEmpty()) return false;

        return true;
    }

    private boolean solveGapN(int scanNum, ExpTag tag, int tagPosInProt, String protSeq, int remainMC, double ms1TolAbs,
                              SparseVector unUsedExpProcessedPL, TreeMap<Double, Double> unUsedPlMap, GRBEnv env, TreeSet<Peptide> nModPepsSet) {
        boolean shouldSolveN = true;
        double nDeltaMass = tag.getHeadLocation() - MassTool.PROTON; //n term delta mass should be independent of cDeltaMass
        double nCutMass = precursorMass + MassTool.PROTON - tag.getHeadLocation() - massTool.H2O;
        if (tag.isNorC == N_TAG || Math.abs(tag.getHeadLocation() - MassTool.PROTON) < ms1TolAbs)  return true;

        List<Integer> nPoses = new ArrayList<>();
        List<Integer> krPoses = new ArrayList<>();
        boolean startRecord = false;
        int mcInN = 0;
        char aaChar;
        boolean isKR;
        for (int nPos = tagPosInProt-1; nPos >= 0; nPos--) {
            aaChar = protSeq.charAt(nPos);
            if (isX(aaChar)) break;
            isKR = isKR(aaChar);
            nDeltaMass -= massTool.getMassTable().get(aaChar);
            if (isKR) mcInN++;
            if (mcInN > remainMC || nDeltaMass < minPtmMass) {
                if (isKR || nPos == 0) {
                    krPoses.add(nPos);
                }
                break;
            }
            if (nDeltaMass < maxPtmMass) {
                if (startRecord) {
                    nPoses.add(nPos);
                    if (isKR)  krPoses.add(nPos); // only count the kr poses in the flexible zone
                } else {
                    startRecord = true;
                }
                if (Math.abs(nDeltaMass) <= ms1TolAbs) {
                    shouldSolveN = false;
                    String nPartSeq = protSeq.substring(nPos, tagPosInProt);
                    storeCleanPartPeptides(nPartSeq, nCutMass, N_PART, unUsedExpProcessedPL, unUsedPlMap, nModPepsSet);
                    break;
                }
            }
        }
        if ( ! startRecord) return false;

        if (shouldSolveN) {
            int optEndPosP1 = -10;
            int optStartPos = -2;
            Map<Integer, Integer> posYIdMap = new HashMap<>();
            Map<Integer, Integer> yIdMinAbsPosMap = new HashMap<>();
            if (buildIndex.nTermSpecific) {
                if (!krPoses.isEmpty()) {
                    optEndPosP1 = Collections.max(krPoses) + 1;
                    optStartPos = Collections.min(krPoses) + 1;

                    if (!nPoses.isEmpty() && nPoses.get(0) == 0) { // if so, no need to care isKR(-1) because there is no aa on -1
                        optStartPos = 0;
                    }
                    int yId = -1;
                    for (int pos = optEndPosP1 - 1; pos > optStartPos - 1; pos--) {
                        if (krPoses.contains(pos)) {
                            yId++;
                        }
                        int tmpMinPos = yIdMinAbsPosMap.getOrDefault(yId, 999999);
                        yIdMinAbsPosMap.put(yId, pos < tmpMinPos ? pos : tmpMinPos);
                        posYIdMap.put(pos, yId);
                    }
                } else {
                    if (nPoses.isEmpty()) {
                        optEndPosP1 = 0;
                        optStartPos = 0;
                    } else if (Collections.min(nPoses) == 0) { //if krPoses is empty, the only feasible situation is min(nPoses) == 0 where there is no digestion at N term
                        optEndPosP1 = 0;
                        optStartPos = 0;
                    }
                }
            } else {
                optEndPosP1 = Collections.max(nPoses) + 1;
                optStartPos = Collections.min(nPoses);
                int yId = -1;
                for (int pos = optEndPosP1 - 1; pos > optStartPos - 1; pos--) {
                    yId++;
                    posYIdMap.put(pos, yId);
                    int tmpMinPos = yIdMinAbsPosMap.getOrDefault(yId, 999999);
                    yIdMinAbsPosMap.put(yId, pos < tmpMinPos ? pos : tmpMinPos);
                }
            }

            if (optStartPos <= optEndPosP1 && shouldSolveN) {
                String nPartSeq = protSeq.substring(optStartPos, tagPosInProt);
                int nPartSeqLen = tagPosInProt - optStartPos;
                double flexiableMass = tag.getHeadLocation() - MassTool.PROTON - massTool.calResidueMass(protSeq.substring(optEndPosP1, tagPosInProt));
                Map<Integer, TreeMap<Double, VarPtm>> absPos_MassVarPtm_Map = new HashMap<>(nPartSeqLen, 1);
                Map<Integer, List<Byte>> absPos_ptmPositions_Map = new HashMap<>(nPartSeqLen, 1);

                inferPTM.prepareInfoNTerm(scanNum, nPartSeq,
                        absPos_MassVarPtm_Map, yIdMinAbsPosMap, optStartPos == 0, optEndPosP1, tagPosInProt, protSeq, absPos_ptmPositions_Map);

                if (nPartSeqLen == 1) {
                    findPtmOnOneAa(nModPepsSet, flexiableMass, absPos_MassVarPtm_Map, nPartSeq, tagPosInProt-1);
                } else {
//                    inferPTM.findBestPtmMIPExtN(scanNum, env, flexiableMass, tagPosInProt, nPartSeq, ms1TolAbs, posYIdMap, absPos_MassVarPtm_Map
//                            , optEndPosP1, optStartPos == 0, yIdMinAbsPosMap, protSeq, unUsedPlMap, nModPepsSet, absPos_ptmPositions_Map);
                    inferPTM.findBestPtmMIPExtNPWL(scanNum, env, flexiableMass, tagPosInProt, nPartSeq, ms1TolAbs, posYIdMap, absPos_MassVarPtm_Map
                            , optEndPosP1, optStartPos == 0, yIdMinAbsPosMap, protSeq, unUsedPlMap, nModPepsSet, absPos_ptmPositions_Map);
                }
                if (yIdMinAbsPosMap.isEmpty()) {
                    inferPTM.findPosssible1Ptm(scanNum, nPartSeq, absPos_MassVarPtm_Map, optStartPos, nModPepsSet, flexiableMass);
                }
            }
        }
        if (nModPepsSet.isEmpty()) return false;

        return true;
    }

    private boolean solveGapMid(int scanNum, ExpTag lTag, int lTagPos, ExpTag rTag, int rTagPos, String gapSeq, double deltaMass, double ms1TolAbs,
                              SparseVector unUsedExpProcessedPL, TreeMap<Double, Double> unUsedPlMap, GRBEnv env, TreeSet<Peptide> midModPepsSet) {

        if (Math.abs(deltaMass) < 0.02) {
            Peptide tmpPeptide = new Peptide(gapSeq, false, massTool);
            tmpPeptide.setScore(0);
            midModPepsSet.add(tmpPeptide);
            return true;
        }
        int gapLen = gapSeq.length();
        Map<Integer, TreeMap<Double, VarPtm>> absPos_MassVarPtm_Map = new HashMap<>(gapLen, 1);
        inferPTM.prepareInfoMid(scanNum, gapSeq, absPos_MassVarPtm_Map, lTagPos+lTag.size());
        if (gapLen == 1) {
            findPtmOnOneAa(midModPepsSet, deltaMass, absPos_MassVarPtm_Map, gapSeq, lTagPos+lTag.size());
        } else {
            // note: these two refMass is similar to cCutMass nCutMass, but different. CutMass are pure sum of aa res masses, the Proton and H2O mass will be added in subsequent Xcorr step
            double bIonRefMass = lTag.getTailLocation();
            double yIonRefMass = precursorMass + 2*MassTool.PROTON - rTag.getHeadLocation();
//            inferPTM.findBestPtmInGap(scanNum, env, deltaMass, lTagPos+lTag.size(), rTagPos,
//                    gapSeq, ms1TolAbs, absPos_MassVarPtm_Map, unUsedPlMap, midModPepsSet, bIonRefMass, yIonRefMass);
            inferPTM.findBestPtmInGapPWL(scanNum, env, deltaMass, lTagPos+lTag.size(), rTagPos,
                    gapSeq, ms1TolAbs, absPos_MassVarPtm_Map, unUsedPlMap, midModPepsSet, bIonRefMass, yIonRefMass, precursorMass);
        }

        // a manual way to find possibly missing 1-ptm solutions.
        inferPTM.findPosssible1Ptm(scanNum, gapSeq, absPos_MassVarPtm_Map, lTagPos+lTag.size(), midModPepsSet, deltaMass);

        if (midModPepsSet.isEmpty()) return false;
        return true;
    }

    private void findPtmOnOneAa(TreeSet<Peptide> midModPepsSet, double deltaMass, Map<Integer, TreeMap<Double, VarPtm>> absPos_MassVarPtm_Map, String gapSeq, int absPos) {
        TreeMap<Double, VarPtm> massPtmMap = absPos_MassVarPtm_Map.get(absPos);
        if (massPtmMap == null) return;
        for (double ptmMass : massPtmMap.keySet()) {
            if (Math.abs(ptmMass - deltaMass) < 2*ms2Tolerance) {
                Peptide tmpPeptide = new Peptide(gapSeq, false, massTool);
                PosMassMap posMassMap = new PosMassMap();
                posMassMap.put(0, ptmMass);
                tmpPeptide.posVarPtmResMap.put(0, massPtmMap.get(ptmMass));
                tmpPeptide.setVarPTM(posMassMap);
                tmpPeptide.setScore(0);
                midModPepsSet.add(tmpPeptide);
            }
        }
    }

    private void updatePeptideTreeSet(Peptide newPeptide, TreeSet<Peptide> peptideTreeSet, Map<String, PeptideInfo> peptideInfoMap,
                                      String protId, String protSeq, int pepStartPos, int pepEndPos) {

        String shortProtId = protId.split(" ")[0];
        boolean added = false;
        if (peptideTreeSet.size() < candisNum) {
            peptideTreeSet.add(newPeptide);
            added = true;
        } else if (newPeptide.compareTo( peptideTreeSet.last() ) == 1) { //use the override '>' to compare because I am considering priority
            added = true;
            peptideTreeSet.pollLast();
            peptideTreeSet.add(newPeptide);
        }
        if (added) {
            char leftFlank  = pepStartPos==0 ? '-' : protSeq.charAt(pepStartPos-1);
            char rightFlank = pepEndPos==protSeq.length()-1 ? '-' : protSeq.charAt(pepEndPos+1);
            String freePepSeq = newPeptide.getFreeSeq();
            PeptideInfo peptideInfo = peptideInfoMap.get(freePepSeq);
            if (peptideInfo != null) {
                if ( ! peptideInfo.protIdSet.contains(shortProtId)) { //if this pep with prot is not recorded, add this prot
                    peptideInfo.leftFlank = leftFlank;
                    peptideInfo.rightFlank = rightFlank;
                    peptideInfo.protIdSet.add(shortProtId);
                    if (!shortProtId.startsWith("DECOY_")) {
                        peptideInfo.isDecoy = false;
                    }
                }
            } else {
                peptideInfo = new PeptideInfo(freePepSeq, shortProtId.startsWith("DECOY_"), leftFlank, rightFlank);
                peptideInfo.protIdSet.add(shortProtId);
                peptideInfoMap.put(freePepSeq, peptideInfo);
            }
        }
    }

    private void storeCleanPartPeptides(String partSeq, double cutMass, byte isNCPart, SparseVector expProcessedPL, TreeMap<Double,Double> plMap, TreeSet<Peptide> modPepsSet) {

        Peptide partPeptide = new Peptide( partSeq, false, massTool);
        double[][] ionMatrix = partPeptide.getIonMatrixNow();
        inferPTM.updateIonMatrix(ionMatrix, cutMass, isNCPart);
        Set<Integer> jRange = IntStream.rangeClosed(0, partSeq.length()-1).boxed().collect(Collectors.toSet());
        double score = massTool.buildVectorAndCalXCorr(ionMatrix, 1, expProcessedPL, partPeptide.matchedBions, partPeptide.matchedYions, jRange) ;
        partPeptide.setScore(score*(1-partPeptide.posVarPtmResMap.size()*0.01));
        partPeptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, ionMatrix, ms2Tolerance));
        if (modPepsSet.size() < 2) { //max restore 2 patterns for one peptide  //todo make this avaliable to differentiate priority
            modPepsSet.add(partPeptide);
        } else if (modPepsSet.last().compareTo(partPeptide) < 0) {
            modPepsSet.pollLast();
            modPepsSet.add(partPeptide);
        }
    }

    private int getNumOfMCFromStr(String tagStr) {
        String str2 = tagStr.replaceAll("[KR]","");
        int numMC = tagStr.length()-str2.length();
        return numMC;
    }
    private boolean isKR(char aa){
        return aa == 'K' || aa == 'R';
    }

    private boolean isX(char aa){
        return aa == 'X';
    }
    public class Entry { //copied from PtmSearch
        public Map<String, PeptideInfo> peptideInfoMapForRef = new HashMap<>();
        public Peptide topPeptide = null;
        final int scanNum;
        final String scanName;
        final int precursorCharge;
        final double precursorMass;
        final String labelling;
        public final String peptide;
        final double theoMass;
        final int isDecoy;
        public final double score;
        final double deltaLCn;
        final double deltaCn;
        final int matchedPeakNum;
        final double ionFrac;
        final double matchedHighestIntensityFrac;
        final double explainedAaFrac;
        final String otherPtmPatterns; // It has 4 decimal because it is write the the result file for checking. It is not used in scoring or other purpose.
        final String aScore;
        final String peptideSet;
        List<VarPtm> varPtmList = new ArrayList<>();
        Entry(int scanNum, String scanName, int precursorCharge, double precursorMass
                ,String labelling, String peptide, double theoMass, int isDecoy
                , double score, double deltaLCn, double deltaCn, int matchedPeakNum, double ionFrac, double matchedHighestIntensityFrac
                , double explainedAaFrac, String otherPtmPatterns, String aScore, String peptideSet) {
            this.scanNum = scanNum;
            this.scanName = scanName;
            this.precursorCharge = precursorCharge;
            this.precursorMass = precursorMass;
            this.labelling = labelling;
            this.peptide = peptide;
            this.theoMass = theoMass;
            this.isDecoy = isDecoy;
            this.score = score;
            this.deltaLCn = deltaLCn;
            this.deltaCn = deltaCn;
            this.matchedPeakNum = matchedPeakNum;
            this.ionFrac = ionFrac;
            this.matchedHighestIntensityFrac = matchedHighestIntensityFrac;
            this.explainedAaFrac = explainedAaFrac;
            this.otherPtmPatterns = otherPtmPatterns;
            this.aScore = aScore;
            this.peptideSet = peptideSet;
        }
    }
}
