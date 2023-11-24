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
import gurobi.GRB;
import gurobi.GRBEnv;
import gurobi.GRBException;
import org.apache.commons.math3.util.Pair;
import proteomics.FM.FMIndex;
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

public final class PreSearch implements Callable<PreSearch.Entry> {
//    private static final Logger logger = LoggerFactory.getLogger(PreSearch.class);
    private static final int candisNum = 20;
    private final BuildIndex buildIndex;
    private final MassTool massTool;
    private final double ms1Tolerance;
//    private final double leftInverseMs1Tolerance;
//    private final double rightInverseMs1Tolerance;
//    private final int ms1ToleranceUnit;
    private final double minPtmMass;
    private final double maxPtmMass;
//    private final int localMaxMs2Charge;
    private final JMzReader spectraParser;
    private final double minClear;
    private final double maxClear;
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
    public PreSearch(int scanNum, BuildIndex buildIndex, MassTool massTool, double ms2Tolerance, double ms1Tolerance, double minPtmMass, double maxPtmMass
            , JMzReader spectraParser, double minClear, double maxClear, ReentrantLock lock, String scanName, int precursorCharge, double precursorMass
            , SpecProcessor specProcessor) {

        this.buildIndex = buildIndex;
        this.massTool = massTool;
        this.ms1Tolerance = ms1Tolerance;
        this.minPtmMass = minPtmMass;
        this.maxPtmMass = maxPtmMass;
//        this.localMaxMs2Charge = localMaxMs2Charge;
        this.spectraParser = spectraParser;
        this.minClear = minClear;
        this.maxClear = maxClear;
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
        } finally {
            lock.unlock();
        }

        double ms1TolAbs = Double.parseDouble(InferPTM.df3.format(precursorMass*ms1Tolerance/1000000));

        TreeMap<Double, Double> plMap = specProcessor.preSpectrumTopNStyleWithChargeLimit(rawPLMap, precursorMass, precursorCharge, minClear, maxClear, DatasetReader.topN, ms2Tolerance);

        if (plMap.isEmpty()) return null;
        // Coding
        SparseVector expProcessedPL;
        if (PIPI.useXcorr) {
            expProcessedPL = specProcessor.prepareXCorr(plMap, false);
        } else {
            expProcessedPL = specProcessor.digitizePL(plMap);
        }
        InferSegment inferSegment = buildIndex.getInferSegment();
        TreeMap<Double, Double> finalPlMap = inferSegment.addVirtualPeaks(precursorMass, plMap);

//        List<ExpTag> allLongTagList = inferSegment.getLongTag(finalPlMap, precursorMass - massTool.H2O + MassTool.PROTON, scanNum, minTagLenToExtract,maxTagLenToExtract);
        List<ExpTag> allLongTagList = inferSegment.getLongTag(finalPlMap, precursorMass - massTool.H2O + MassTool.PROTON, scanNum, minTagLenToExtract,maxTagLenToExtract);


        List<ExpTag> cleanedAllLongTagList = inferSegment.cleanAbundantTagsPrefix(allLongTagList, minTagLenToExtract);

        allLongTagList = cleanedAllLongTagList;
        double s = massTool.getScoreWithPtmSeq("SEDTISK(14.016)M(15.995)NDFM(15.995)R", expProcessedPL);
        if (lszDebugScanNum.contains(this.scanNum)) {
            int a = 1;
        }
        if (allLongTagList.isEmpty())  return null;

        double totalMass = precursorMass + 2 * MassTool.PROTON;


        Map<String, PeptideInfo> peptideInfoMap = new HashMap<>(50000);

        int minTagLen = 3;

        GRBEnv env = new GRBEnv(true);
//        env.set(GRB.IntParam.OutputFlag,0);
//        env.set(GRB.IntParam.LogToConsole, 0);
        env.start();


        int tagNumsToTest =  Math.min(20, allLongTagList.size());
        List<ExpTag> tagsToTest = new ArrayList<>();
        Set<Integer> distinctTagIds = new HashSet<>();
        TreeSet<Peptide> peptideTreeSet = new TreeSet<>(Collections.reverseOrder());
        ExpTag tagInfo;
        int t = 0;
        for (int i = 0; i < tagNumsToTest; i++) {
            tagInfo = allLongTagList.get(i);
            if (tagInfo.isNorC == N_TAG) { //n tag //todo can I relax this. can I check all N and C tags as non NC (both direction)
                tagsToTest.add(tagInfo);
                distinctTagIds.add(t);
                t++;
            } else if (tagInfo.isNorC == C_TAG) { // c tag
                ExpTag revTagInfo = tagInfo.revTag(totalMass);
                tagsToTest.add(revTagInfo);
                distinctTagIds.add(t);
                t++;
            } else { // non-nc tag
                tagsToTest.add(tagInfo);
                distinctTagIds.add(t);
                t++;
                ExpTag revTagInfo = tagInfo.revTag(totalMass);
                tagsToTest.add(revTagInfo);
                distinctTagIds.add(t);
                t++;
            }
        }
        List<OccGroup> occGroupList = new ArrayList<>(100);
//        Map<Pair<Integer, Integer>, Integer> posRangeOccIdMap = new HashMap<>();
        int numResSub = searchAndSaveFuzzy1(scanNum , tagsToTest, distinctTagIds, occGroupList, minTagLen, totalMass);
        occGroupList.sort(Comparator.comparing(o->o.totalScore, Comparator.reverseOrder()));

        int count = 0;
        double topScore = 0;
        outLoop:
        for (OccGroup occG : occGroupList) {
            List<Pair<ExpTag, Integer>> tagRelPosList = occG.tagRelPosList;
            tagRelPosList.sort(Comparator.comparing(o->o.getFirst().getTotalIntensity(),Comparator.reverseOrder()));

            String protId = occG.protId;

            for (Pair<ExpTag, Integer> tagRelPos : tagRelPosList) {
                double thisScore = tagRelPos.getFirst().getTotalIntensity();
                if (count == 0) topScore = thisScore;
                count++;
                if (count > 10 && thisScore < topScore) break outLoop;
                updateCandiListMaxP(scanNum, protId, tagRelPos.getSecond(), tagRelPos.getFirst(), ms1TolAbs, peptideTreeSet, peptideInfoMap, expProcessedPL, plMap, env);
            }
        }
        if (lszDebugScanNum.contains(scanNum)){
            int a = 1;
        }
//        StringBuilder a = new StringBuilder("sadsa");
//        a.in
        if (peptideTreeSet.isEmpty()) {
            env.dispose();
            return null;
        }
        List<Peptide> pepList = new ArrayList<>(peptideTreeSet);
        Set<Integer> pepIdsToRemove = new HashSet<>();
        for (int j = pepList.size()-1; j >= 1; j--) {
            if (pepIdsToRemove.contains(j)) continue;
            for (int i = 0; i < j; i++) {
                if (pepIdsToRemove.contains(i)) continue;
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

        boolean shouldPtm = Math.abs(precursorMass-massTool.calResidueMass(topPep.getFreeSeq()) - massTool.H2O) > 0.01;
        boolean hasPTM = topPep.hasVarPTM();
        int ptmNum = 0;
        boolean isSettled = true;
        double totalPtmMass = 0;
        if (hasPTM) {
            ptmNum = topPep.getVarPTMs().size();
            for (double mass : topPep.getVarPTMs().values()){
                totalPtmMass += mass;
            }
        }
        isSettled = Math.abs(totalPtmMass-(precursorMass-massTool.calResidueMass(topPep.getFreeSeq()) - massTool.H2O)) <= 0.01;

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
        entry = new PreSearch.Entry(
                scanNum, scanName, shouldPtm ? 1 : 0, hasPTM ? 1 : 0, ptmNum, isSettled ? 1 : 0
                , precursorCharge, precursorMass, buildIndex.getLabelling(), topPep.getPtmContainingSeq(buildIndex.returnFixModMap())
                , topPep.getTheoMass(), topPep.isDecoy() ? 1 : 0, topPep.getScore(), deltaLCn, deltaCn
                , topPep.getMatchedPeakNum(), topPep.getIonFrac(), topPep.getMatchedHighestIntensityFrac()
                , topPep.getExplainedAaFrac(), otherPtmPatterns, topPep.getaScore(), ""
                , pepSetString.substring(0, pepSetString.length()-1)
        );
        entry.topPeptide = topPep;
        entry.topPeptide.precursorMass = precursorMass;
        entry.varPtmList.addAll(topPep.posVarPtmResMap.values());
        for (Peptide peptide : peptideArray) {
            entry.peptideInfoMapForRef.put(peptide.getFreeSeq(), peptideInfoMap.get(peptide.getFreeSeq()));
        }
//        System.out.println(scanNum + ","+entry.peptideInfoMapForRef.size() + "," + peptideArray.length);l
        if (lszDebugScanNum.contains(scanNum)){
            int a = 1;
            System.out.println(scanNum + "," + pepSetString);
        }
        int c = 1;
//        env.release();l
        env.dispose();
        return entry;
    }

    final class OccGroup {
        public List<Pair<ExpTag, Integer>> tagRelPosList = new ArrayList<>();
        public int lPos;
        public int rPos;
        public int spanLen;
        final public String protId;
        public double totalScore;
        public OccGroup(ExpTag initialTag, String protId, int tagPosInProt) {
            this.tagRelPosList.add(new Pair<>(initialTag, tagPosInProt));
            this.spanLen = initialTag.size();
            this.lPos = tagPosInProt;
            this.rPos = tagPosInProt + this.spanLen;
            this.protId = protId;
            this.totalScore = initialTag.getTotalIntensity();
        }
        public void addTag(ExpTag newTag, int tagPosInProt) {
            boolean shouldAdd = true;
            for (Pair<ExpTag, Integer> tagRelPosPair : this.tagRelPosList) {
                if (tagRelPosPair.getFirst().getFreeAaString().contentEquals(newTag.getFreeAaString())
                        && tagRelPosPair.getSecond() == tagPosInProt
                        && tagRelPosPair.getFirst().getHeadLocation() == newTag.getHeadLocation()) {
                    shouldAdd = false;
                    break;
                }
            }
            if (! shouldAdd) return;
            int tagLen = newTag.size();
            for (int i = 0; i < tagLen; i++) {   // for exptag out of the current range, add there intens into totalScore
                if ( i+tagPosInProt < this.lPos) {
                    this.totalScore += newTag.getExpAAs()[i].getHeadIntensity();
                    if (i == tagLen-1) {
                        this.totalScore += newTag.getExpAAs()[i].getTailIntensity();
                    }
                }
                if ( i+tagPosInProt > this.rPos) {
                    this.totalScore += newTag.getExpAAs()[i].getTailIntensity();
                    if (i == 0) {
                        this.totalScore += newTag.getExpAAs()[i].getHeadIntensity();
                    }
                }
            }
            this.tagRelPosList.add(new Pair<>(newTag, tagPosInProt));
            this.lPos = Math.min(tagPosInProt, this.lPos);
            this.rPos = Math.max(tagPosInProt + tagLen, this.rPos);
            this.spanLen = Math.min(this.rPos - this.lPos, this.spanLen+ tagLen);
        }
    }
    private int searchAndSaveFuzzy1(int scanNum, List<ExpTag> tagsToTest, Set<Integer> distinctTagIds
            , List<OccGroup> occGroupList, int minTagLen, double totalMass) {

        char[] tagChar;
        int numRes;
        int solCount = 0;
        FMRes fmRes;
        int matchedPos;
        String protId;
        int tagLen;
        boolean isSubLast;
        ExpTag tagInfo;
        Map<Pair<Integer, Integer>, Integer> posRangeOccIdMap = new HashMap<>();
        int tagsNum = tagsToTest.size();
        Set<Pair<String, Integer>> matchedProtIdRelPosSet = new HashSet<>();
        int occId = -1;
        for (int i = 0; i < tagsNum; i++) {
            boolean shouldTryForwardSearch = false;
            tagInfo = tagsToTest.get(i);
            if (distinctTagIds.contains(i)) matchedProtIdRelPosSet.clear();
            tagChar = tagInfo.getFreeAaString().toCharArray();

            fmRes = buildIndex.fmIndexNormal.fmSearchFuzzy(tagChar);
            solCount += fmRes.ep-fmRes.sp+1;
            matchedPos = fmRes.matchedPos;
            int oriIsNorC = tagInfo.isNorC;
            if (matchedPos != 0) {
                shouldTryForwardSearch = true;
                if (oriIsNorC == N_TAG) {
                    tagInfo.isNorC = NON_NC_TAG;
                }
            }
            tagLen = tagInfo.size();
            if (matchedPos > 0 && tagInfo.size()-matchedPos < minTagLen) continue;

            ExpTag resTagInfo = matchedPos == 0 ? tagInfo : tagInfo.subTag(matchedPos, tagLen);
            int absTagPos;
            int dotIndex;
            int lPos;
            int rPos;
            for (int ii = fmRes.sp; ii <= fmRes.ep; ii++) {
                absTagPos = buildIndex.fmIndexNormal.SA[ii];
                dotIndex = -2-Arrays.binarySearch(buildIndex.dotPosArrNormal, absTagPos);
                protId = buildIndex.posProtMapNormal.get(dotIndex);
                lPos = absTagPos - buildIndex.dotPosArrNormal[dotIndex] - 1;
                rPos = lPos + tagLen - matchedPos;
                Pair<String, Integer> protIdRelPos = new Pair<>(protId, lPos);
                if (matchedProtIdRelPosSet.contains(protIdRelPos)) {
                    continue;
                } else {
                    matchedProtIdRelPosSet.add(protIdRelPos);
                }

                boolean foundRange = false;
                Iterator<Map.Entry<Pair<Integer, Integer>,Integer>> iterator = posRangeOccIdMap.entrySet().iterator();
                while (iterator.hasNext()) {
                    Map.Entry<Pair<Integer, Integer>,Integer> entry = iterator.next();
                    int lSpan = entry.getKey().getFirst();
                    int rSpan = entry.getKey().getSecond();
                    if (lPos >= lSpan && rPos <= rSpan) {
                        if (occGroupList.get(entry.getValue()).protId.contentEquals(protId)) {
                            occId = entry.getValue();
                            OccGroup occG = occGroupList.get(occId);
                            occG.addTag(resTagInfo, lPos);
                            iterator.remove(); // remove it for later update because it will be modified
                            posRangeOccIdMap.put(new Pair<>(Math.min(occG.lPos-occG.spanLen, lSpan), Math.max(occG.rPos+occG.spanLen, rSpan)), occId);
                            foundRange = true;
                            break;
                        }
                    }
                }

                if ( ! foundRange) {  //new occurrence that does not belong to any pos range
                    occId++;
                    OccGroup occG = new OccGroup(resTagInfo, protId, lPos);
                    occGroupList.add(occG);
                    posRangeOccIdMap.put(new Pair<>(occG.lPos-occG.spanLen, occG.rPos+occG.spanLen), occId);
                }
            }
            if (shouldTryForwardSearch) {
                ExpTag revTagInfo = tagInfo.revTag(totalMass);
                if (distinctTagIds.contains(i)) matchedProtIdRelPosSet.clear();
                tagChar = revTagInfo.getFreeAaString().toCharArray();
                fmRes = buildIndex.fmIndexReverse.fmSearchFuzzy(tagChar);
                solCount += fmRes.ep-fmRes.sp+1;
                matchedPos = fmRes.matchedPos;

                if (matchedPos != 0) {
                    if (oriIsNorC == C_TAG) {
                        tagInfo.isNorC = NON_NC_TAG;
                    }
                    if (oriIsNorC == N_TAG) {
                        tagInfo.isNorC = N_TAG;
                    }
                }
                if (matchedPos > 0 && tagInfo.size()-matchedPos < minTagLen) continue;

//                ExpTag resRevTagInfo = matchedPos == 0 ? revTagInfo : revTagInfo.subTag(matchedPos, tagLen);
                ExpTag resRevTagInfo = matchedPos == 0 ? tagInfo : tagInfo.subTag(0, tagLen-matchedPos);

                for (int ii = fmRes.sp; ii <= fmRes.ep; ii++) {
                    absTagPos = buildIndex.fmIndexReverse.SA[ii];
                    absTagPos = buildIndex.textNormalLen - absTagPos - tagLen + matchedPos; // should I minus one more?
                    dotIndex = -2-Arrays.binarySearch(buildIndex.dotPosArrNormal, absTagPos);
                    protId = buildIndex.posProtMapNormal.get(dotIndex);
                    lPos = absTagPos - buildIndex.dotPosArrNormal[dotIndex] - 1;
                    rPos = lPos + tagLen - matchedPos;
                    Pair<String, Integer> protIdRelPos = new Pair<>(protId, lPos);
                    if (matchedProtIdRelPosSet.contains(protIdRelPos)) {
                        continue;
                    } else {
                        matchedProtIdRelPosSet.add(protIdRelPos);
                    }

                    boolean foundRange = false;
                    Iterator<Map.Entry<Pair<Integer, Integer>,Integer>> iterator = posRangeOccIdMap.entrySet().iterator();
                    while (iterator.hasNext()) {
                        Map.Entry<Pair<Integer, Integer>,Integer> entry = iterator.next();
                        int lSpan = entry.getKey().getFirst();
                        int rSpan = entry.getKey().getSecond();
                        if (lPos >= lSpan && rPos <= rSpan) {
                            if (occGroupList.get(entry.getValue()).protId.contentEquals(protId)) {
                                occId = entry.getValue();
                                OccGroup occG = occGroupList.get(occId);
                                occG.addTag(resRevTagInfo, lPos);
                                iterator.remove(); // remove it for later update because it will be modified
                                posRangeOccIdMap.put(new Pair<>(Math.min(occG.lPos-occG.spanLen, lSpan), Math.max(occG.rPos+occG.spanLen, rSpan)), occId);
                                foundRange = true;
                                break;
                            }
                        }
                    }

                    if ( ! foundRange) {  //new occurrence that does not belong to any pos range
                        occId++;
                        OccGroup occG = new OccGroup(resRevTagInfo, protId, lPos);
                        occGroupList.add(occG);
                        posRangeOccIdMap.put(new Pair<>(occG.lPos-occG.spanLen, occG.rPos+occG.spanLen), occId);
                    }
                }
            }
        }
        return solCount;
    }
//    private int searchAndSaveFuzzy1(int scanNum, List<ExpTag> tagsToTest, Set<Integer> distinctTagIds, double ms1TolAbs, TreeSet<Peptide> peptideTreeSet, Map<String, PeptideInfo> peptideInfoMap
//            , List<OccGroup> occGroupList, FMIndex fmIndex, int minTagLen, SparseVector expProcessedPL,TreeMap<Double, Double> plMap) throws CloneNotSupportedException {
//
//        char[] tagChar;
//        int numRes;
//        int solCount = 0;
//        FMRes fmRes;
//        int matchedPos;
//        String protId;
//        int tagLen;
//        boolean isSubLast;
//        ExpTag tagInfo;
//        Map<Pair<Integer, Integer>, Integer> posRangeOccIdMap = new HashMap<>();
//        int tagsNum = tagsToTest.size();
//        Set<Pair<String, Integer>> matchedProtIdRelPosSet = new HashSet<>();
//        int occId = -1;
//        for (int i = 0; i < tagsNum; i++) {
//            tagInfo = tagsToTest.get(i);
//            if (distinctTagIds.contains(i)) matchedProtIdRelPosSet.clear();
//            forwardSearch()
//
//        }
//        return solCount;
//    }

//    private void forwardSearch(ExpTag tagInfo, FMIndex fmIndex, int minTagLen, Set<Pair<String, Integer>> matchedProtIdRelPosSet){
//        char[] tagChar = tagInfo.getFreeAaString().toCharArray();
//        FMRes fmRes = fmIndex.fmSearchFuzzy(tagChar);
//        if (fmRes == null) return;
//
//        int matchedPos = fmRes.matchedPos;
//        int tagLen = tagChar.length;
//        if (matchedPos > 0 && tagLen-matchedPos < minTagLen) return;
//
//        ExpTag resTagInfo = matchedPos == 0 ? tagInfo : tagInfo.subTag(matchedPos, tagLen);
//        int absTagPos;
//        int dotIndex;
//        int lPos;
//        int rPos;
//        String protId;
//        int occId = -1;
//        for (int ii = fmRes.sp; ii <= fmRes.ep; ii++) {
//            absTagPos = fmIndex.SA[ii];
//            dotIndex = -2-Arrays.binarySearch(buildIndex.dotPosArrNormal, absTagPos);
//            protId = buildIndex.posProtMapNormal.get(dotIndex);
//            lPos = absTagPos - buildIndex.dotPosArrNormal[dotIndex] - 1;
//            rPos = lPos + tagLen;
//            Pair<String, Integer> protIdRelPos = new Pair<>(protId, lPos);
//            if (matchedProtIdRelPosSet.contains(protIdRelPos)) {
//                continue;
//            } else {
//                matchedProtIdRelPosSet.add(protIdRelPos);
//            }
//
//            boolean foundRange = false;
//            Iterator<Map.Entry<Pair<Integer, Integer>,Integer>> iterator = posRangeOccIdMap.entrySet().iterator();
//            while (iterator.hasNext()) {
//                Map.Entry<Pair<Integer, Integer>,Integer> entry = iterator.next();
//                int lSpan = entry.getKey().getFirst();
//                int rSpan = entry.getKey().getSecond();
//                if (lPos >= lSpan && rPos <= rSpan) {
//                    if (occGroupList.get(entry.getValue()).protId.contentEquals(protId)) {
//                        occId = entry.getValue();
//                        OccGroup occG = occGroupList.get(occId);
//                        occG.addTag(resTagInfo, lPos);
//                        iterator.remove(); // remove it for later update because it will be modified
//                        posRangeOccIdMap.put(new Pair<>(Math.min(occG.lPos-occG.spanLen, lSpan), Math.max(occG.rPos+occG.spanLen, rSpan)), occId);
//                        foundRange = true;
//                        break;
//                    }
//                }
//            }
//
//            if ( ! foundRange) {  //new occurrence that does not belong to any pos range
//                occId++;
//                OccGroup occG = new OccGroup(resTagInfo, protId, lPos);
//                occGroupList.add(occG);
//                posRangeOccIdMap.put(new Pair<>(occG.lPos-occG.spanLen, occG.rPos+occG.spanLen), occId);
//            }
//        }
//    }
    private boolean isHomo(Peptide p1, Peptide p2, Map<String, PeptideInfo> peptideInfoMap) {
//        int a = 1;
        Set<String> temp1 = peptideInfoMap.get(p1.getFreeSeq()).protIdSet;
        Set<String> temp2 = peptideInfoMap.get(p2.getFreeSeq()).protIdSet;

        Set<String> set = new HashSet<>(temp1);
        set.retainAll(temp2);
//        if (set.isEmpty()) return false;
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
    private boolean alreadySearched(ExpTag tagInfo, Set<ExpTag> usedTagInfoSet) {
        boolean res = false;
        for (ExpTag tmpTag : usedTagInfoSet) {
            if (tmpTag.size() < tagInfo.size()) continue;

            int alignedPos = tmpTag.getFreeAaString().indexOf(tagInfo.getFreeAaString());
            int endPos = alignedPos + tagInfo.size()-1;

            if (alignedPos >= 0
                    && Math.abs(tmpTag.expAaList.get(alignedPos).getHeadLocation() - tagInfo.getHeadLocation()) < 0.02){
                res = true;
                break;
            } else {
                int n_sharedAa = LCS(tagInfo.getFreeAaString(), tmpTag.getFreeAaString()).length();
                if (n_sharedAa > tagInfo.size()*0.7){
                    res = true;
                    break;
                }
            }
        }
        return res;
    }

//    private boolean sameTagStrUsed(String tagStr, Set<String> usedTagStrSet) {
//        return res;
//    }

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

    private boolean isCMassValid(double mass) {
        double deltaMass = mass - precursorMass;
        if (deltaMass > minPtmMass - 530 && mass < maxPtmMass) {
            return true;
        }
        return false;
    }

    private boolean isNMassValid(double mass) {
        double deltaMass = mass - precursorMass;
        if (deltaMass > minPtmMass - 530 && mass < maxPtmMass) {
            return true;
        }
        return false;
    }

    private void updateCandiListMaxP(int scanNum, String protId, int tagPosInProt, ExpTag finderTag, double ms1TolAbs, TreeSet<Peptide> peptideTreeSet
            , Map<String, PeptideInfo> peptideInfoMap, SparseVector expProcessedPL, TreeMap<Double, Double> plMap, GRBEnv env) {

        if (lszDebugScanNum.contains(scanNum) ) {
            int a = 1;
        }
        TreeMap<Double, Double> unUsedPlMap = new TreeMap<>(plMap);
        for (ExpAa expAa : finderTag.expAaList) {
            if (plMap.containsKey(expAa.getHeadLocation())) {
                unUsedPlMap.remove(expAa.getHeadLocation());
            }
            if (plMap.containsKey(expAa.getTailLocation())) {
                unUsedPlMap.remove(expAa.getTailLocation());
            }
        }
        SparseVector unUsedExpProcessedPL;
        unUsedExpProcessedPL = specProcessor.digitizePL(unUsedPlMap);
        double tagCMass = finderTag.getTailLocation() + massTool.H2O-MassTool.PROTON; // +massTool.H2O-MassTool.PROTON  is to mimic the mass of the real neutral precursor mass
        double cCutMass = finderTag.getTailLocation() - MassTool.PROTON;
        String protSeq = buildIndex.protSeqMap.get(protId);
        int protLen = protSeq.length();
        int finderTagSize = finderTag.size();
        if ( (tagPosInProt+finderTag.size() == protLen && finderTag.isNorC != C_TAG)  // not C tag but found at prot Cterm, impossible
                || (tagPosInProt == 0 && finderTag.isNorC != N_TAG) // not N tag but found at prot Nterm, impossible
                || tagPosInProt < 0
                || tagPosInProt+finderTag.size() > protLen
        ) {
            return;
        }
        int remainMC = massTool.missedCleavage - getNumOfMissCleavSite(finderTag);
        if (lszDebugScanNum.contains(scanNum) && finderTag.getFreeAaString().contentEquals("TLD")){ // SLEEEGAA
            int a = 1;
        }
        boolean shouldSolveC = true;
        TreeSet<Peptide> cModPepsSet = new TreeSet<>(Comparator.reverseOrder());
        if (finderTag.isNorC == C_TAG || Math.abs(tagCMass - precursorMass) < ms1TolAbs){ // C tag, tagCMass must be exactly settle. Otherwise no valid cPos will be, then no valid n pos.
            shouldSolveC = false;
        } else { // only when oriTag is not C oriTag can it extend to c
            List<Integer> cPoses = new ArrayList<>();
            List<Integer> krPoses = new ArrayList<>();
            int mcInC = 0;
            boolean startRecord = false;
            char aaChar;
            boolean isKR;
            for (int cPos = tagPosInProt+finderTagSize; cPos < protLen; cPos++) {  //
                aaChar = protSeq.charAt(cPos);
                if (isX(aaChar))  break;
                isKR = isKR(aaChar);
                tagCMass += massTool.getMassTable().get(aaChar); //even when tag is fuzzy, tagcmass wont be disturbed
                if (mcInC > remainMC || tagCMass > precursorMass+maxPtmMass) break; // the total max minus ptm is -250. although the max minus ptm for single ptm is -156
                if (tagCMass >= precursorMass-maxPtmMass) { // the total max plus ptm is 600. though the max plus ptm for single ptm is 527
                    if (startRecord) {
                        cPoses.add(cPos);
                    } else {
                        startRecord = true;
                    }
                    if (isKR)  krPoses.add(cPos);

                    if (Math.abs(tagCMass - precursorMass) <= ms1TolAbs) {
                        shouldSolveC = false;
                        String cPartSeq = protSeq.substring(tagPosInProt+finderTagSize, cPos+1);
                        storeCleanPartPeptides(cPartSeq, cCutMass, C_PART, unUsedExpProcessedPL, unUsedPlMap, cModPepsSet);
                        break;
                    }
                }
                if (isKR) mcInC++;
            }// finish collect all possible cposes
            if ( shouldSolveC && !startRecord) {
                return;   // shouldSolveC but cPoses is empty, just return
            }
            if (lszDebugScanNum.contains(scanNum) && finderTag.getFreeAaString().contentEquals("VAV")){ // SLEEEGAA
                int a = 1;
            }
            if (shouldSolveC) {
                Map<Integer, Integer> posYIdMap = new HashMap<>();
                Map<Integer, Integer> yIdMaxAbsPosMap = new HashMap<>();
                int optStartPos = 1;
                int optEndPosP1 = 0;
                if (cTermSpecific) {
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

                byte isProtNorC_Term = NON_TERM_PROT;
                if (optStartPos == 1 && protSeq.charAt(0) == 'M') {
                    isProtNorC_Term = N_TERM_PROT;
                }
                if (optStartPos <= optEndPosP1) {
                    int cPartStartPos = tagPosInProt + finderTagSize;
                    String cPartSeq = protSeq.substring(cPartStartPos, optEndPosP1);
                    int cPartSeqLen = cPartSeq.length();
                    if (lszDebugScanNum.contains(scanNum)) {
                        if (cPartSeq.contentEquals("SSAQQVAVSR")) { // SLEEEGAA
                            int a = 1;
                        }
                    }
                    double flexiableMass = precursorMass - (finderTag.getTailLocation() + massTool.H2O-MassTool.PROTON) - massTool.calResidueMass(protSeq.substring(cPartStartPos, optStartPos));
//                    Map<Integer, Set<Double>> oneTimeMassGroups = new HashMap<>(cPartSeqLen);
//                    Map<Set<Integer>, Set<Double>> posComb_multiMassSet_Map = new HashMap<>(200);
                    Map<Integer, TreeMap<Double, VarPtm>> absPos_MassVarPtm_Map = new HashMap<>(cPartSeqLen, 1);
//                    Map<Double, Set<Integer>> allMassAllPosesMap = new HashMap<>(cPartSeqLen * 30, 1); //Set<pos>

                    if (cPartSeq.contentEquals("KLYTGAGYDEVK")) {
                        int a =1;
                    }
                    inferPTM.prepareInfoCTerm(scanNum, cPartSeq, absPos_MassVarPtm_Map,  yIdMaxAbsPosMap, optEndPosP1 == protLen, optStartPos, cPartStartPos);

                    Set<Pair<Integer, Map<Double, Integer>>> resList = new HashSet<>(100);
                    Set<Map<Integer, Double>> posPtmIdResList = new HashSet<>();

                    if (!yIdMaxAbsPosMap.isEmpty()) {
                        int a= 1;
                    }
                    inferPTM.findBestPtmMIPExtC(scanNum, env, flexiableMass, cPartStartPos,
                            cPartSeq, ms1TolAbs, posYIdMap, absPos_MassVarPtm_Map, resList, isProtNorC_Term, optStartPos, optEndPosP1 == protLen
                            , yIdMaxAbsPosMap, posPtmIdResList, unUsedPlMap, cModPepsSet);

                    if (lszDebugScanNum.contains(scanNum) && cPartSeq.contentEquals("SSAQQVAVSR")){ // SLEEEGAA
                        int a = 1;
                    }
//                    inferPTM.getFeasibleMassPosMapC(scanNum, resList, unUsedPlMap, cPartSeq, cCutMass, C_PART,
//                            expProcessedPL, false, allMassAllPosesMap, absPos_MassVarPtm_Map, cModPepsSet, cPartStartPos, yIdMaxAbsPosMap, optStartPos);
                }
            }
        } // end settle C section
        if (shouldSolveC && cModPepsSet.isEmpty()) {
            return;
        }

        if (lszDebugScanNum.contains(scanNum) && finderTag.getFreeAaString().contentEquals("VDLDAPDVSL")){ //0,VDLDAPDVSL, SLEEEGAA
            int a = 1;
        }
//        protSeq = "PEPTLDEKLVFSGNLFQHQEDSKKLQDELQNMKEEMARLVSGKDYNVTANSKLGSPSGKTGTCRMSFVKGWGAEYRLLEQKTQESQKGFCFLTYTDEEPVKKTQAYQDQKPGTSGLRLLNEPTAAALAYGLDKKLLDLKTGTVKAANMLQQSGSKNTGAKAGPNASLLSLKSDKGADFLVTEVENGGSLGSKKDASTLQSQKAEGTGDAKGTVPDDAVEALADSLGKKSKHLQEQLNELKLLEEEDSKLKLVEKYGYTHLSAGELLRVNVEAPDVNLEGLGGKLKSEASSEFAKKDTPTSAGPNSFNKGKGTLSAPGKVVTAAAQAKAVVKTDLQLALPSGCYGRLSALQGKLSKSLFFGGKGAPLQRATVASSTQKFQDLGVKASNAAVKLAESKDHLSDKESKAAQASDLEKLHLDEKLGNDFHTNKRDQEALMKSVKVNLGQGSHPQKVKSQNTDMVQKSVSKVGLDTPDLDLHGPEGKLKLLTWDVKDTLLRVDLDVPDVNLEGPDAKLKFLQEFYQDDELGKKDALHFYNKSLAEHRVYLASSSGSTALKKGWLKSNVSDAVAQSTRLFQSDTNAMLGKKDLGKFQVATDALKSGNLTEDDKHNNAKLSSNFSSLLAEKLRAMEAASSLKSQAATKKYEELDNAPEERTLQKQSVVYGGKGETASKLQSELSRTKMSAAQALKSPSDSGYSYETLGKTTKLPQEKCLLQTDVKLTQDLMAKVRLEEVTGKLQVARLQTKAGEATVKLDLPAENSNSETLLLTGKRVNALKNLQVKHELQANCYEEVKDRLTASLCDLKSRSSFAYKDQNENRFSMPGFKGEGPDVDVSLPKLMDQNLKCLSAAEEKAGQKLLDVNHYAKLDHLAEKFRALASLDSKLNQAKVEDCKMESTETEERSSGPGGQNVNKVNSKDSLNFLDHLNGKKYLELYVQKVNPSRTGGTQTDLFTCGKCKSLKLNTAEDAKSPDEAYALAKKSDSGKPYYYNSQTKVDVEVPDVSLEGPEGKLKVFLGNLNTLVVKKKNLSQELNMEQKLTSLNVKYNNDKLLLPGELAKHAVSEGTKSGFGKFERLDEKENLSAKVQGGALEDSQLVAGVAFKKPGGEGPQCEKTTDVKANLEKAQAELVGTADEATRPANEKATDDYHYEKSGNFSAAMKDLSGKAALNQKLLETGERTEVLSPNSKVESKPSVSTKLLSKGSLGAQKLANTCFNELEKGTLTVSAQELKDNRVFSLLSSEKELKFSMPGFKGEGPDVDVNLPKDDSFLGKLGGTLARVGDDLAKATGDWKDNSTMGYMAAKKDLSTVEALQNLKNLKFDDTNPEKEEAKDAEAAEATAEGALKAEKYEELVKEVSTYLKMESEEGKEARALLQAKASGKASDTSSETVFGKRMAQELKYGDLVERDFDTALKHYDKLDEMPEAAVKSTANKVGLQVVAVKAPGFGDNRDKALPLEAEVAPVKASEALLKQLKLQLELDQKKDYSAPVNFLSAGLKKFGTLNLVHPKLGSYTKLQAAYAGDKADDLQKTPGQNAQKWLPARLFVGGLKEDTEEHHLRLQEKVESAQSEQKLEKTLDDLEEKDSSPTTNSKLSALQRSLGEPAAGEGSAKRLLSSLEQKEENKTTEETSKPKLHNAENLQPGEQKYEYKYEQEHAALQDKLFQVAKSEDTLSKMNDFMRLTESVAETAQTLKKVSPASSVDSNLPSSQGYKKTGASWTDNLMAQKCSKDESFLGKLGGTLARDYSSGFGGKYGVQADRTFGKAAGPSLSHTSGGTQSKVETTVTSLKTKTTGFYSGFSEVAEKRLTLSKHQNVQLPRVGTASVLQPVKKVLTANSNPSSPSAAKRANEKTESSSAQQVAVSRLSSAMSAAKALCDHVRLQDVSGQLSSSKKLTGKNQVTATKTQDQLSNLKYHEEFEKDVNAALAALKTKSPPLPEHQKLPCNSAEPKLAEEENKAKSSVACKWNLAEAQQKPSVGSQSNQAGQGKRTEDEVLTSKGDAWAKAQQHWGSGVGVKKKFGVLSDNFKVDLDAPDVSLEGPDAKLKLLFSNTAAQKLRVLDSGAPLKLPVGPETLGRSNEEGSEEKGPEVRSGDWKGYTGKLTMLNTVSKLRFSAYLKNSNPALNDNLEKDEPSVAAMVYPFTGDHKQKMGSKSPGNTSQPPAFFSKVTLEVGKVLQQGRVTAAAMAGNKSTPRDTMGLADKTENTLERGAEKTEVAEPRSMSTEGLMKFVDSKVGEFSGANKEKASALSLKLGSSKGSLEKGSPEDKGEVGAGAGPGAQAGPSAKRHAEMVHTGLKLERSLTHSTVGSKGEKGGNLGAKWTLDLKLGNSYFKEEKSDDESSQKDLKDLVAALQHNYKMSAFKLNNFSADLKDSKLEQVEKELLRVPTTEKPTVTVNFRVGVEVPDVNLEGPEGKLKLLLVSTTPYSEKDTKFELGKLMELHGEGSSSGKSTGLELETPSLVPVKKLQQELQTAKKPSEGPLQSVQVFGRSLEEEGAAHSGKRSKLLFSNTAAQKTPNLYLYSKKDTGEKLTVAENEAETKLLEEELQAPTSSKRLGVATVAKGVCGACKMMSKPQTSGAYVLNKGGPGSTLSFVGKRAANDAGYFNDEMAPLEVKTKHAVSEGTKAVTKTVSKVDDFLANEAKSGKFTQQDLDEAKGFALVGVGSEASSKKMQQLKEQYRGEVQVSDKERTLDELLQEKKLTEKELAEAASKLHVQSSSDSSDEPAEKRSQPATKTEDYGMGPGRLPLANTEKYMADKPASGAGAGAGAGKRKPEPTLDE";
        TreeSet<Peptide> nModPepsSet = new TreeSet<>(Comparator.reverseOrder());
        boolean shouldSolveN = true;
        double nDeltaMass = finderTag.getHeadLocation() - MassTool.PROTON; //n term delta mass should be independent of cDeltaMass
        double nCutMass = precursorMass + MassTool.PROTON - finderTag.getHeadLocation() - massTool.H2O;
        if (finderTag.isNorC == N_TAG || Math.abs(finderTag.getHeadLocation() - MassTool.PROTON) < ms1TolAbs) {
            shouldSolveN = false;
        } else {
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
            if (lszDebugScanNum.contains(scanNum) ){ //0,VDLDAPDVSL, SLEEEGAA
                int a = 1;
            }
            if ( shouldSolveN && ! startRecord) {
                return;   // shouldSolveN but nPoses is empty, just return
            }
            if (shouldSolveN) {
                int optEndPosP1 = -10;
                int optStartPos = -2;
                Map<Integer, Integer> posYIdMap = new HashMap<>();
                Map<Integer, Integer> yIdMinAbsPosMap = new HashMap<>();
                if (nTermSpecific) {
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
                    double flexiableMass = finderTag.getHeadLocation() - MassTool.PROTON - massTool.calResidueMass(protSeq.substring(optEndPosP1, tagPosInProt));

                    Map<Integer, TreeMap<Double, VarPtm>> pos_MassVarPtm_Map = new HashMap<>(nPartSeqLen, 1);

                    inferPTM.prepareInfoNTerm(scanNum, nPartSeq,
                            pos_MassVarPtm_Map, yIdMinAbsPosMap, optStartPos == 0, optEndPosP1, tagPosInProt, protSeq);

                    Set<Pair<Integer, Map<Double, Integer>>> resList = new HashSet<>(100);
                    Set<Map<Integer, Double>> posPtmIdResSet = new HashSet<>();
//                    if (nPartSeq.contentEquals("LTE")) {
//                        int a = 1;
//                    }
                    if (!yIdMinAbsPosMap.isEmpty()) {
                        int a= 1;
                    }
                    inferPTM.findBestPtmMIPExtN(scanNum, env, flexiableMass, tagPosInProt,
                            nPartSeq, ms1TolAbs , posYIdMap, pos_MassVarPtm_Map, optEndPosP1,optStartPos == 0, yIdMinAbsPosMap, posPtmIdResSet,protSeq, unUsedPlMap, nModPepsSet, unUsedExpProcessedPL, nCutMass);
//                    inferPTM.getFeasibleMassPosMapN(scanNum, resList, plMap, nPartSeq, nCutMass, N_PART,
//                            expProcessedPL, false, allMassAllPosesMap, pos_MassVarPtm_Map, nModPepsSet, tagPosInProt, yIdMinAbsPosMap, optEndPosP1);
                }
            }
        }// end settle N sectionl
        if (shouldSolveN && nModPepsSet.isEmpty()) {
            return;
        }
        if (lszDebugScanNum.contains(scanNum) && finderTag.getFreeAaString().contentEquals("SLEEEGAA")){ //0,VDLDAPDVSL, SLEEEGAA
            int a = 1;
        }

        if ( cModPepsSet.isEmpty() ) {
            //only n
            double lastScore = -1;
            for (Peptide nPartpeptide : nModPepsSet) {
                if (nPartpeptide.getScore() < lastScore) {
                    continue;   // only takes top1 or binglie top1
                }
                String freePepSeq = nPartpeptide.getFreeSeq()+finderTag.getFreeAaString();
                Peptide fullPeptide = new Peptide(freePepSeq, nPartpeptide.isDecoy, massTool);
                fullPeptide.finderTag = finderTag;

                PosMassMap fullPosMassMap;
                if (nPartpeptide.hasVarPTM()){
                    fullPosMassMap = nPartpeptide.getVarPTMs().clone(); // maybe replace PosMass at all
                    fullPeptide.posVarPtmResMap.putAll(nPartpeptide.posVarPtmResMap);
                } else {
                    fullPosMassMap = new PosMassMap();
                }
                int idOfAa = nPartpeptide.length()-1;
                for (char aaChar : finderTag.getPtmAaString().toCharArray()) {
                    if (Character.isUpperCase(aaChar)) {
                        idOfAa += 1;
                    } else {
                        fullPosMassMap.put(idOfAa, massTool.labelVarPtmMap.get(aaChar).mass);
                        fullPeptide.posVarPtmResMap.put(idOfAa, massTool.labelVarPtmMap.get(aaChar));
                    }
                }
                if (! fullPosMassMap.isEmpty()) fullPeptide.setVarPTM(fullPosMassMap);
                double calScore = massTool.buildVectorAndCalXCorr(fullPeptide.getIonMatrixNow(), 1, expProcessedPL, fullPeptide.matchedBions, fullPeptide.matchedYions);//todo decide the penalty
                fullPeptide.setScore(calScore*(1 - fullPeptide.posVarPtmResMap.size() * 0.05));
                fullPeptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, fullPeptide.getIonMatrixNow(), ms2Tolerance));
                updatePeptideTreeSet(fullPeptide, peptideTreeSet, peptideInfoMap, protId, protSeq, tagPosInProt-nPartpeptide.length(), tagPosInProt+finderTag.size()-1);
            }
        } else if ( nModPepsSet.isEmpty() ) {
            //only c
            double lastScore = -1;
            for (Peptide cPartpeptide : cModPepsSet) {
                if (cPartpeptide.getScore() < lastScore) {
                    continue;   // only takes top1 or binglie top1
                }
                String freePepSeq = finderTag.getFreeAaString() + cPartpeptide.getFreeSeq();
                Peptide fullPeptide = new Peptide(freePepSeq, cPartpeptide.isDecoy, massTool);
                fullPeptide.finderTag = finderTag;

                PosMassMap fullPosMassMap = new PosMassMap(); // maybe replace PosMass at all
                int idOfAa = cPartpeptide.length()-1;
                for (char aaChar : finderTag.getPtmAaString().toCharArray()) {
                    if (Character.isUpperCase(aaChar)) {
                        idOfAa += 1;
                    } else {
                        fullPosMassMap.put(idOfAa, massTool.labelVarPtmMap.get(aaChar).mass);
                        fullPeptide.posVarPtmResMap.put(idOfAa, massTool.labelVarPtmMap.get(aaChar));
                    }
                }
                for (int pos : cPartpeptide.posVarPtmResMap.keySet()) { //copy the top 1 ptm pattern in n part // whhat if also choose the largeset priority
                    fullPosMassMap.put(pos + finderTag.size(), cPartpeptide.posVarPtmResMap.get(pos).mass); // copy the ptms from partModPepsUnsettled
                    fullPeptide.posVarPtmResMap.put(pos + finderTag.size(), cPartpeptide.posVarPtmResMap.get(pos));
                }

                if (! fullPosMassMap.isEmpty()) fullPeptide.setVarPTM(fullPosMassMap);
                double calScore = massTool.buildVectorAndCalXCorr(fullPeptide.getIonMatrixNow(), 1, expProcessedPL, fullPeptide.matchedBions, fullPeptide.matchedYions);//todo decide the penalty
                fullPeptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, fullPeptide.getIonMatrixNow(), ms2Tolerance));
                fullPeptide.setScore(calScore*(1 - fullPeptide.posVarPtmResMap.size() * 0.05));

                updatePeptideTreeSet(fullPeptide, peptideTreeSet, peptideInfoMap, protId, protSeq, tagPosInProt, tagPosInProt+finderTag.size()+cPartpeptide.length()-1);
            }
        } else if ( ! cModPepsSet.isEmpty() && ! nModPepsSet.isEmpty() ){
            List<NCid> ncIdList = new ArrayList<>(4);
            int n = 0;
            int c;
            for (Peptide nPeptide : nModPepsSet) {
                c = 0;
                for (Peptide cPeptide : cModPepsSet) {
                    ncIdList.add(new NCid(n,c,nPeptide.getScore()+cPeptide.getScore()));
                    c++;
                }
                n++;
            }
            ncIdList.sort(Comparator.comparing(o->o.score, Comparator.reverseOrder()));

            Peptide[] nPeptideArray = new Peptide[nModPepsSet.size()];
            nModPepsSet.toArray(nPeptideArray);
            Peptide[] cPeptideArray = new Peptide[cModPepsSet.size()];
            cModPepsSet.toArray(cPeptideArray);
//            double lastScore = -1;
            for (NCid ncId : ncIdList.subList(0,Math.min(2, ncIdList.size()))) {
                Peptide nPartpeptide = nPeptideArray[ncId.nId];
                Peptide cPartpeptide = cPeptideArray[ncId.cId];
//                if (nPartpeptide.getScore()+cPartpeptide.getScore() < lastScore) {
//                    continue;
//                }
                String freePepSeq = nPartpeptide.getFreeSeq()+finderTag.getFreeAaString()+cPartpeptide.getFreeSeq();
                Peptide fullPeptide = new Peptide(freePepSeq, nPartpeptide.isDecoy, massTool);// todo
                fullPeptide.finderTag = finderTag;

                PosMassMap fullPosMassMap;
                if (nPartpeptide.hasVarPTM()){
                    fullPosMassMap = nPartpeptide.getVarPTMs().clone(); // maybe replace PosMass at all
                    fullPeptide.posVarPtmResMap.putAll(nPartpeptide.posVarPtmResMap);
                } else {
                    fullPosMassMap = new PosMassMap();
                }
                int idOfAa = nPartpeptide.length()-1;
                for (char aaChar : finderTag.getPtmAaString().toCharArray()) {
                    if (Character.isUpperCase(aaChar)) {
                        idOfAa += 1;
                    } else {
                        fullPosMassMap.put(idOfAa, massTool.labelVarPtmMap.get(aaChar).mass);
                        fullPeptide.posVarPtmResMap.put(idOfAa, massTool.labelVarPtmMap.get(aaChar));
                    }
                }
                for (int pos : cPartpeptide.posVarPtmResMap.keySet()) { //copy the top 1 ptm pattern in n part // whhat if also choose the largeset priority
                    fullPosMassMap.put(pos + finderTag.size() + nPartpeptide.length(), cPartpeptide.posVarPtmResMap.get(pos).mass); // copy the ptms from partModPepsUnsettled
                    fullPeptide.posVarPtmResMap.put(pos + finderTag.size() + nPartpeptide.length(), cPartpeptide.posVarPtmResMap.get(pos));
                }

                if (! fullPosMassMap.isEmpty()) fullPeptide.setVarPTM(fullPosMassMap);
                double calScore = massTool.buildVectorAndCalXCorr(fullPeptide.getIonMatrixNow(), 1, expProcessedPL, fullPeptide.matchedBions, fullPeptide.matchedYions);//todo decide the penalty
                fullPeptide.setScore(calScore*(1 - fullPeptide.posVarPtmResMap.size() * 0.05));
                fullPeptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, fullPeptide.getIonMatrixNow(), ms2Tolerance));
                updatePeptideTreeSet(fullPeptide, peptideTreeSet, peptideInfoMap, protId, protSeq, tagPosInProt-nPartpeptide.length(), tagPosInProt+finderTag.size()+cPartpeptide.length()-1);
            }
        }
    }

//    private void updateCandiList(int scanNum, String protId, int tagPosInProt, ExpTag finderTag, double ms1TolAbs, TreeSet<Peptide> peptideTreeSet
//            , Map<String, PeptideInfo> peptideInfoMap, SparseVector expProcessedPL, TreeMap<Double, Double> plMap, GRBEnv env) {
//
//        if (lszDebugScanNum.contains(scanNum) && protId.contentEquals("Packet: PTM_Trainingskit_Kmod")) {
//            int a = 1;
//        }
//        TreeMap<Double, Double> unUsedPlMap = new TreeMap<>(plMap);
//        for (ExpAa expAa : finderTag.expAaList) {
//            if (plMap.containsKey(expAa.getHeadLocation())) {
//                unUsedPlMap.remove(expAa.getHeadLocation());
//            }
//            if (plMap.containsKey(expAa.getTailLocation())) {
//                unUsedPlMap.remove(expAa.getTailLocation());
//            }
//        }
//        SparseVector unUsedExpProcessedPL;
//        unUsedExpProcessedPL = specProcessor.digitizePL(unUsedPlMap);
//        double tagCMass = finderTag.getTailLocation() + massTool.H2O-MassTool.PROTON; // +massTool.H2O-MassTool.PROTON  is to mimic the mass of the real neutral precursor mass
//        double cCutMass = finderTag.getTailLocation() - MassTool.PROTON;
//        String protSeq = buildIndex.protSeqMap.get(protId);
//        int protLen = protSeq.length();
//        int finderTagSize = finderTag.size();
//        if ( (tagPosInProt+finderTag.size() == protLen && finderTag.isNorC != C_TAG)  // not C tag but found at prot Cterm, impossible
//                || (tagPosInProt == 0 && finderTag.isNorC != N_TAG) // not N tag but found at prot Nterm, impossible
//                || tagPosInProt < 0
//                || tagPosInProt+finderTag.size() > protLen
//        ) {
//            return;
//        }
//        int remainMC = massTool.missedCleavage - getNumOfMissCleavSite(finderTag);
//        if (lszDebugScanNum.contains(scanNum) && finderTag.getFreeAaString().contentEquals("TLD")){ // SLEEEGAA
//            int a = 1;
//        }
//        boolean shouldSolveC = true;
//        TreeSet<Peptide> cModPepsSet = new TreeSet<>(Comparator.reverseOrder());
//        if (finderTag.isNorC == C_TAG || Math.abs(tagCMass - precursorMass) < ms1TolAbs){ // C tag, tagCMass must be exactly settle. Otherwise no valid cPos will be, then no valid n pos.
//            shouldSolveC = false;
//        } else { // only when oriTag is not C oriTag can it extend to c
//            List<Integer> cPoses = new ArrayList<>();
//            List<Integer> krPoses = new ArrayList<>();
//            int mcInC = 0;
//            boolean startRecord = false;
//            char aaChar;
//            boolean isKR;
//            for (int cPos = tagPosInProt+finderTagSize; cPos < protLen; cPos++) {  //
//                aaChar = protSeq.charAt(cPos);
//                if (isX(aaChar))  break;
//                isKR = isKR(aaChar);
//                tagCMass += massTool.getMassTable().get(aaChar); //even when tag is fuzzy, tagcmass wont be disturbed
//                if (mcInC > remainMC || tagCMass > precursorMass+maxPtmMass) break; // the total max minus ptm is -250. although the max minus ptm for single ptm is -156
//                if (tagCMass >= precursorMass-maxPtmMass) { // the total max plus ptm is 600. though the max plus ptm for single ptm is 527
//                    if (startRecord) {
//                        cPoses.add(cPos);
//                    } else {
//                        startRecord = true;
//                    }
//                    if (isKR)  krPoses.add(cPos);
//
//                    if (Math.abs(tagCMass - precursorMass) <= ms1TolAbs) {
//                        shouldSolveC = false;
//                        String cPartSeq = protSeq.substring(tagPosInProt+finderTagSize, cPos+1);
//                        storeCleanPartPeptides(cPartSeq, cCutMass, C_PART, unUsedExpProcessedPL, unUsedPlMap, cModPepsSet);
//                        break;
//                    }
//                }
//                if (isKR) mcInC++;
//            }// finish collect all possible cposes
//            if ( shouldSolveC && !startRecord) {
//                return;   // shouldSolveC but cPoses is empty, just return
//            }
//            if (lszDebugScanNum.contains(scanNum) && finderTag.getFreeAaString().contentEquals("TQAYQ")){ // SLEEEGAA
//                int a = 1;
//            }
//            if (shouldSolveC) {
//                Map<Integer, Integer> posYIdMap = new HashMap<>();
//                Map<Integer, Integer> yIdMaxAbsPosMap = new HashMap<>();
//                int optStartPos = 1;
//                int optEndPosP1 = 0;
//                if (cTermSpecific) {
//                    if ( ! krPoses.isEmpty()) {
//                        optEndPosP1 = krPoses.get(krPoses.size()-1) + 1;  //max kr
//                        optStartPos = krPoses.get(0) + 1;  // good trick//min kr
//                        if ( ! cPoses.isEmpty() && cPoses.get(cPoses.size()-1) == protLen - 1) {
//                            optEndPosP1 = protLen;
//                        }
//                        int yId = 0;
//                        for (int pos = optStartPos; pos < optEndPosP1; pos++) {
//                            posYIdMap.put(pos, yId); // using rel pos in part seq i.e. start from 0 not the prot pos
//                            int oldMaxPos = yIdMaxAbsPosMap.getOrDefault(yId, -99);
//                            yIdMaxAbsPosMap.put(yId, pos > oldMaxPos ? pos : oldMaxPos);
//                            if (krPoses.contains(pos)) {
//                                yId++;
//                            }
//                        }
//                    } else {
//                        if (cPoses.isEmpty()) {
//                            optEndPosP1 = protLen;
//                            optStartPos = protLen;  // good trick
//                        } else if (cPoses.get(cPoses.size()-1) == protLen-1) {
//                            optEndPosP1 = protLen;
//                            optStartPos = protLen;  // good trick
//                        }
//                    }
//                } else {
//                    optEndPosP1 = cPoses.get(cPoses.size()-1) + 1; //max
//                    optStartPos = cPoses.get(0);  // good trick  //min
//                    int yId = 0;
//                    for (int pos = optStartPos; pos < optEndPosP1; pos++) {
//                        posYIdMap.put(pos, yId); // using rel pos in part seq i.e. start from 0 not the prot pos
//                        int oldMaxPos = yIdMaxAbsPosMap.getOrDefault(yId, -99);
//                        yIdMaxAbsPosMap.put(yId, pos > oldMaxPos ? pos : oldMaxPos);
//                        yId++;
//                    }
//                }
//                if (optStartPos <= optEndPosP1) {
//                    int cPartStartPos = tagPosInProt+finderTagSize;
//                    String cPartSeq = protSeq.substring(cPartStartPos, optEndPosP1);
//                    int cPartSeqLen = cPartSeq.length();
//                    if (lszDebugScanNum.contains(scanNum) ){ // SLEEEGAA
//                        int a = 1;
//                    }
//                    double flexiableMass = precursorMass - (finderTag.getTailLocation() + massTool.H2O-MassTool.PROTON) - massTool.calResidueMass(protSeq.substring(cPartStartPos, optStartPos));
//                    Map<Integer, Set<Double>> oneTimeMassGroups = new HashMap<>(cPartSeqLen);
//                    Map<Set<Integer>, Set<Double>> posComb_multiMassSet_Map = new HashMap<>(200);
//                    Map<Integer, Map<Double, VarPtm>> pos_MassVarPtm_Map = new HashMap<>(cPartSeqLen, 1);
//                    Map<Double, Set<Integer>> allMassAllPosesMap = new HashMap<>(cPartSeqLen * 30, 1); //Set<pos>
//
//
//                    inferPTM.prepareInfoCTerm(scanNum, cPartSeq, oneTimeMassGroups, posComb_multiMassSet_Map,
//                            pos_MassVarPtm_Map, allMassAllPosesMap, yIdMaxAbsPosMap, optEndPosP1 == protLen, optStartPos, cPartStartPos);
//
//                    Set<Pair<Integer, Map<Double, Integer>>> resList = new HashSet<>(100);
//                    inferPTM.findBestPtmMIPExtC(scanNum, env, allMassAllPosesMap, flexiableMass, (cPartStartPos),
//                            cPartSeq, ms1TolAbs, oneTimeMassGroups, posComb_multiMassSet_Map, posYIdMap, pos_MassVarPtm_Map, resList);
//
//                    inferPTM.getFeasibleMassPosMapC(scanNum, resList, unUsedPlMap, cPartSeq, cCutMass, C_PART,
//                            expProcessedPL, false, allMassAllPosesMap, pos_MassVarPtm_Map, cModPepsSet, cPartStartPos, yIdMaxAbsPosMap, optStartPos);
//
//
//
//                }
//            }
//        } // end settle C section
//        if (shouldSolveC && cModPepsSet.isEmpty()) {
//            return;
//        }
//
//        if (lszDebugScanNum.contains(scanNum) && finderTag.getFreeAaString().contentEquals("VDLDAPDVSL")){ //0,VDLDAPDVSL, SLEEEGAA
//            int a = 1;
//        }
////        protSeq = "PEPTLDEKLVFSGNLFQHQEDSKKLQDELQNMKEEMARLVSGKDYNVTANSKLGSPSGKTGTCRMSFVKGWGAEYRLLEQKTQESQKGFCFLTYTDEEPVKKTQAYQDQKPGTSGLRLLNEPTAAALAYGLDKKLLDLKTGTVKAANMLQQSGSKNTGAKAGPNASLLSLKSDKGADFLVTEVENGGSLGSKKDASTLQSQKAEGTGDAKGTVPDDAVEALADSLGKKSKHLQEQLNELKLLEEEDSKLKLVEKYGYTHLSAGELLRVNVEAPDVNLEGLGGKLKSEASSEFAKKDTPTSAGPNSFNKGKGTLSAPGKVVTAAAQAKAVVKTDLQLALPSGCYGRLSALQGKLSKSLFFGGKGAPLQRATVASSTQKFQDLGVKASNAAVKLAESKDHLSDKESKAAQASDLEKLHLDEKLGNDFHTNKRDQEALMKSVKVNLGQGSHPQKVKSQNTDMVQKSVSKVGLDTPDLDLHGPEGKLKLLTWDVKDTLLRVDLDVPDVNLEGPDAKLKFLQEFYQDDELGKKDALHFYNKSLAEHRVYLASSSGSTALKKGWLKSNVSDAVAQSTRLFQSDTNAMLGKKDLGKFQVATDALKSGNLTEDDKHNNAKLSSNFSSLLAEKLRAMEAASSLKSQAATKKYEELDNAPEERTLQKQSVVYGGKGETASKLQSELSRTKMSAAQALKSPSDSGYSYETLGKTTKLPQEKCLLQTDVKLTQDLMAKVRLEEVTGKLQVARLQTKAGEATVKLDLPAENSNSETLLLTGKRVNALKNLQVKHELQANCYEEVKDRLTASLCDLKSRSSFAYKDQNENRFSMPGFKGEGPDVDVSLPKLMDQNLKCLSAAEEKAGQKLLDVNHYAKLDHLAEKFRALASLDSKLNQAKVEDCKMESTETEERSSGPGGQNVNKVNSKDSLNFLDHLNGKKYLELYVQKVNPSRTGGTQTDLFTCGKCKSLKLNTAEDAKSPDEAYALAKKSDSGKPYYYNSQTKVDVEVPDVSLEGPEGKLKVFLGNLNTLVVKKKNLSQELNMEQKLTSLNVKYNNDKLLLPGELAKHAVSEGTKSGFGKFERLDEKENLSAKVQGGALEDSQLVAGVAFKKPGGEGPQCEKTTDVKANLEKAQAELVGTADEATRPANEKATDDYHYEKSGNFSAAMKDLSGKAALNQKLLETGERTEVLSPNSKVESKPSVSTKLLSKGSLGAQKLANTCFNELEKGTLTVSAQELKDNRVFSLLSSEKELKFSMPGFKGEGPDVDVNLPKDDSFLGKLGGTLARVGDDLAKATGDWKDNSTMGYMAAKKDLSTVEALQNLKNLKFDDTNPEKEEAKDAEAAEATAEGALKAEKYEELVKEVSTYLKMESEEGKEARALLQAKASGKASDTSSETVFGKRMAQELKYGDLVERDFDTALKHYDKLDEMPEAAVKSTANKVGLQVVAVKAPGFGDNRDKALPLEAEVAPVKASEALLKQLKLQLELDQKKDYSAPVNFLSAGLKKFGTLNLVHPKLGSYTKLQAAYAGDKADDLQKTPGQNAQKWLPARLFVGGLKEDTEEHHLRLQEKVESAQSEQKLEKTLDDLEEKDSSPTTNSKLSALQRSLGEPAAGEGSAKRLLSSLEQKEENKTTEETSKPKLHNAENLQPGEQKYEYKYEQEHAALQDKLFQVAKSEDTLSKMNDFMRLTESVAETAQTLKKVSPASSVDSNLPSSQGYKKTGASWTDNLMAQKCSKDESFLGKLGGTLARDYSSGFGGKYGVQADRTFGKAAGPSLSHTSGGTQSKVETTVTSLKTKTTGFYSGFSEVAEKRLTLSKHQNVQLPRVGTASVLQPVKKVLTANSNPSSPSAAKRANEKTESSSAQQVAVSRLSSAMSAAKALCDHVRLQDVSGQLSSSKKLTGKNQVTATKTQDQLSNLKYHEEFEKDVNAALAALKTKSPPLPEHQKLPCNSAEPKLAEEENKAKSSVACKWNLAEAQQKPSVGSQSNQAGQGKRTEDEVLTSKGDAWAKAQQHWGSGVGVKKKFGVLSDNFKVDLDAPDVSLEGPDAKLKLLFSNTAAQKLRVLDSGAPLKLPVGPETLGRSNEEGSEEKGPEVRSGDWKGYTGKLTMLNTVSKLRFSAYLKNSNPALNDNLEKDEPSVAAMVYPFTGDHKQKMGSKSPGNTSQPPAFFSKVTLEVGKVLQQGRVTAAAMAGNKSTPRDTMGLADKTENTLERGAEKTEVAEPRSMSTEGLMKFVDSKVGEFSGANKEKASALSLKLGSSKGSLEKGSPEDKGEVGAGAGPGAQAGPSAKRHAEMVHTGLKLERSLTHSTVGSKGEKGGNLGAKWTLDLKLGNSYFKEEKSDDESSQKDLKDLVAALQHNYKMSAFKLNNFSADLKDSKLEQVEKELLRVPTTEKPTVTVNFRVGVEVPDVNLEGPEGKLKLLLVSTTPYSEKDTKFELGKLMELHGEGSSSGKSTGLELETPSLVPVKKLQQELQTAKKPSEGPLQSVQVFGRSLEEEGAAHSGKRSKLLFSNTAAQKTPNLYLYSKKDTGEKLTVAENEAETKLLEEELQAPTSSKRLGVATVAKGVCGACKMMSKPQTSGAYVLNKGGPGSTLSFVGKRAANDAGYFNDEMAPLEVKTKHAVSEGTKAVTKTVSKVDDFLANEAKSGKFTQQDLDEAKGFALVGVGSEASSKKMQQLKEQYRGEVQVSDKERTLDELLQEKKLTEKELAEAASKLHVQSSSDSSDEPAEKRSQPATKTEDYGMGPGRLPLANTEKYMADKPASGAGAGAGAGKRKPEPTLDE";
//        TreeSet<Peptide> nModPepsSet = new TreeSet<>(Comparator.reverseOrder());
//        boolean shouldSolveN = true;
//        double nDeltaMass = finderTag.getHeadLocation() - MassTool.PROTON; //n term delta mass should be independent of cDeltaMass
//        double nCutMass = precursorMass + MassTool.PROTON - finderTag.getHeadLocation() - massTool.H2O;
//        if (finderTag.isNorC == N_TAG || Math.abs(finderTag.getHeadLocation() - MassTool.PROTON) < ms1TolAbs) {
//            shouldSolveN = false;
//        } else {
//            List<Integer> nPoses = new ArrayList<>();
//            List<Integer> krPoses = new ArrayList<>();
//            boolean startRecord = false;
//            int mcInN = 0;
//            char aaChar;
//            boolean isKR;
//            for (int nPos = tagPosInProt-1; nPos >= 0; nPos--) {
//                aaChar = protSeq.charAt(nPos);
//                if (isX(aaChar)) break;
//                isKR = isKR(aaChar);
//
//                nDeltaMass -= massTool.getMassTable().get(aaChar);
//                if (isKR) mcInN++;
//                if (mcInN > remainMC || nDeltaMass < minPtmMass) {
//                    if (isKR || nPos == 0) {
//                        krPoses.add(nPos);
//                    }
//                    break;
//                }
//
//                if (nDeltaMass < maxPtmMass) {
//                    if (startRecord) {
//                        nPoses.add(nPos);
//                        if (isKR)  krPoses.add(nPos); // only count the kr poses in the flexible zone
//                    } else {
//                        startRecord = true;
//                    }
//
//                    if (Math.abs(nDeltaMass) <= ms1TolAbs) {
//                        shouldSolveN = false;
//                        String nPartSeq = protSeq.substring(nPos, tagPosInProt);
//                        storeCleanPartPeptides(nPartSeq, nCutMass, N_PART, unUsedExpProcessedPL, unUsedPlMap, nModPepsSet);
//                        break;
//                    }
//                }
//            }
//            if (lszDebugScanNum.contains(scanNum) ){ //0,VDLDAPDVSL, SLEEEGAA
//                int a = 1;
//            }
//            if ( shouldSolveN && ! startRecord) {
//                return;   // shouldSolveN but nPoses is empty, just return
//            }
//            if (shouldSolveN) {
//                int optEndPosP1 = -10;
//                int optStartPos = -2;
//                Map<Integer, Integer> posYIdMap = new HashMap<>();
//                Map<Integer, Integer> yIdMinAbsPosMap = new HashMap<>();
//                if (nTermSpecific) {
//                    if (!krPoses.isEmpty()) {
//                        optEndPosP1 = Collections.max(krPoses) + 1;
//                        optStartPos = Collections.min(krPoses) + 1;
//
//                        if (!nPoses.isEmpty() && nPoses.get(0) == 0) { // if so, no need to care isKR(-1) because there is no aa on -1
//                            optStartPos = 0;
//                        }
//                        int yId = -1;
//                        for (int pos = optEndPosP1 - 1; pos > optStartPos - 1; pos--) {
//                            if (krPoses.contains(pos)) {
//                                yId++;
//                            }
//                            int tmpMinPos = yIdMinAbsPosMap.getOrDefault(yId, 999999);
//                            yIdMinAbsPosMap.put(yId, pos < tmpMinPos ? pos : tmpMinPos);
//                            posYIdMap.put(pos, yId);
//                        }
//                    } else {
//                        if (nPoses.isEmpty()) {
//                            optEndPosP1 = 0;
//                            optStartPos = 0;
//                        } else if (Collections.min(nPoses) == 0) { //if krPoses is empty, the only feasible situation is min(nPoses) == 0 where there is no digestion at N term
//                            optEndPosP1 = 0;
//                            optStartPos = 0;
//                        }
//                    }
//                } else {
//                    optEndPosP1 = Collections.max(nPoses) + 1;
//                    optStartPos = Collections.min(nPoses);
//                    int yId = -1;
//                    for (int pos = optEndPosP1 - 1; pos > optStartPos - 1; pos--) {
//                        yId++;
//                        posYIdMap.put(pos, yId);
//                        int tmpMinPos = yIdMinAbsPosMap.getOrDefault(yId, 999999);
//                        yIdMinAbsPosMap.put(yId, pos < tmpMinPos ? pos : tmpMinPos);
//                    }
//                }
//
//                if (optStartPos <= optEndPosP1 && shouldSolveN) {
//                    String nPartSeq = protSeq.substring(optStartPos, tagPosInProt);
//                    int nPartSeqLen = tagPosInProt - optStartPos;
//                    double flexiableMass = finderTag.getHeadLocation() - MassTool.PROTON - massTool.calResidueMass(protSeq.substring(optEndPosP1, tagPosInProt));
//
//                    Map<Integer, Set<Double>> oneTimeMassGroups = new HashMap<>(nPartSeqLen);
//                    Map<Set<Integer>, Set<Double>> posComb_multiMassSet_Map = new HashMap<>(200);
//                    Map<Integer, Map<Double, VarPtm>> pos_MassVarPtm_Map = new HashMap<>(nPartSeqLen, 1);
//                    Map<Double, Set<Integer>> allMassAllPosesMap = new HashMap<>(nPartSeqLen * 20, 1); //Set<pos>
//
//                    inferPTM.prepareInfoNTerm(scanNum, nPartSeq, oneTimeMassGroups, posComb_multiMassSet_Map,
//                            pos_MassVarPtm_Map, allMassAllPosesMap, yIdMinAbsPosMap, optStartPos == 0, optEndPosP1, tagPosInProt, protSeq);
//
//                    Set<Pair<Integer, Map<Double, Integer>>> resList = new HashSet<>(100);
//                    inferPTM.findBestPtmMIPExtN(scanNum, env, allMassAllPosesMap, flexiableMass, tagPosInProt,
//                            nPartSeq, ms1TolAbs, oneTimeMassGroups, posComb_multiMassSet_Map, posYIdMap, pos_MassVarPtm_Map, resList);
//                    try {
//                        inferPTM.getFeasibleMassPosMapN(scanNum, resList, plMap, nPartSeq, nCutMass, N_PART,
//                                expProcessedPL, false, allMassAllPosesMap, pos_MassVarPtm_Map, nModPepsSet, tagPosInProt, yIdMinAbsPosMap, optEndPosP1);
//                    } catch (Exception e) {
//                        System.out.println(scanNum + " ," + finderTag.getFreeAaString() + "," + protId);
//                    }
//                }
//            }
//        }// end settle N sectionl
//        if (shouldSolveN && nModPepsSet.isEmpty()) {
//            return;
//        }
//        if (lszDebugScanNum.contains(scanNum) && finderTag.getFreeAaString().contentEquals("SLEEEGAA")){ //0,VDLDAPDVSL, SLEEEGAA
//            int a = 1;
//        }
//
//        if ( cModPepsSet.isEmpty() ) {
//            //only n
//            double lastScore = -1;
//            for (Peptide nPartpeptide : nModPepsSet) {
//                if (nPartpeptide.getScore() < lastScore) {
//                    continue;   // only takes top1 or binglie top1
//                }
//                String freePepSeq = nPartpeptide.getFreeSeq()+finderTag.getFreeAaString();
//                Peptide fullPeptide = new Peptide(freePepSeq, nPartpeptide.isDecoy, massTool);
//                fullPeptide.finderTag = finderTag;
//
//                PosMassMap fullPosMassMap;
//                if (nPartpeptide.hasVarPTM()){
//                    fullPosMassMap = nPartpeptide.getVarPTMs().clone(); // maybe replace PosMass at all
//                    fullPeptide.posVarPtmResMap.putAll(nPartpeptide.posVarPtmResMap);
//                } else {
//                    fullPosMassMap = new PosMassMap();
//                }
//                int idOfAa = nPartpeptide.length()-1;
//                for (char aaChar : finderTag.getPtmAaString().toCharArray()) {
//                    if (Character.isUpperCase(aaChar)) {
//                        idOfAa += 1;
//                    } else {
//                        fullPosMassMap.put(idOfAa, massTool.labelVarPtmMap.get(aaChar).mass);
//                        fullPeptide.posVarPtmResMap.put(idOfAa, massTool.labelVarPtmMap.get(aaChar));
//                    }
//                }
//                if (! fullPosMassMap.isEmpty()) fullPeptide.setVarPTM(fullPosMassMap);
//                double calScore = massTool.buildVectorAndCalXCorr(fullPeptide.getIonMatrixNow(), 1, expProcessedPL, fullPeptide.matchedBions, fullPeptide.matchedYions);//todo decide the penalty
//                fullPeptide.setScore(calScore*(1 - fullPeptide.posVarPtmResMap.size() * 0.05));
//                fullPeptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, fullPeptide.getIonMatrixNow(), ms2Tolerance));
//                updatePeptideTreeSet(fullPeptide, peptideTreeSet, peptideInfoMap, protId, protSeq, tagPosInProt-nPartpeptide.length(), tagPosInProt+finderTag.size()-1);
//            }
//        } else if ( nModPepsSet.isEmpty() ) {
//            //only c
//            double lastScore = -1;
//            for (Peptide cPartpeptide : cModPepsSet) {
//                if (cPartpeptide.getScore() < lastScore) {
//                    continue;   // only takes top1 or binglie top1
//                }
//                String freePepSeq = finderTag.getFreeAaString() + cPartpeptide.getFreeSeq();
//                Peptide fullPeptide = new Peptide(freePepSeq, cPartpeptide.isDecoy, massTool);
//                fullPeptide.finderTag = finderTag;
//
//                PosMassMap fullPosMassMap = new PosMassMap(); // maybe replace PosMass at all
//                int idOfAa = cPartpeptide.length()-1;
//                for (char aaChar : finderTag.getPtmAaString().toCharArray()) {
//                    if (Character.isUpperCase(aaChar)) {
//                        idOfAa += 1;
//                    } else {
//                        fullPosMassMap.put(idOfAa, massTool.labelVarPtmMap.get(aaChar).mass);
//                        fullPeptide.posVarPtmResMap.put(idOfAa, massTool.labelVarPtmMap.get(aaChar));
//                    }
//                }
//                for (int pos : cPartpeptide.posVarPtmResMap.keySet()) { //copy the top 1 ptm pattern in n part // whhat if also choose the largeset priority
//                    fullPosMassMap.put(pos + finderTag.size(), cPartpeptide.posVarPtmResMap.get(pos).mass); // copy the ptms from partModPepsUnsettled
//                    fullPeptide.posVarPtmResMap.put(pos + finderTag.size(), cPartpeptide.posVarPtmResMap.get(pos));
//                }
//
//                if (! fullPosMassMap.isEmpty()) fullPeptide.setVarPTM(fullPosMassMap);
//                double calScore = massTool.buildVectorAndCalXCorr(fullPeptide.getIonMatrixNow(), 1, expProcessedPL, fullPeptide.matchedBions, fullPeptide.matchedYions);//todo decide the penalty
//                fullPeptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, fullPeptide.getIonMatrixNow(), ms2Tolerance));
//                fullPeptide.setScore(calScore*(1 - fullPeptide.posVarPtmResMap.size() * 0.05));
//
//                updatePeptideTreeSet(fullPeptide, peptideTreeSet, peptideInfoMap, protId, protSeq, tagPosInProt, tagPosInProt+finderTag.size()+cPartpeptide.length()-1);
//            }
//        } else if ( ! cModPepsSet.isEmpty() && ! nModPepsSet.isEmpty() ){
//            List<NCid> ncIdList = new ArrayList<>(4);
//            int n = 0;
//            int c;
//            for (Peptide nPeptide : nModPepsSet) {
//                c = 0;
//                for (Peptide cPeptide : cModPepsSet) {
//                    ncIdList.add(new NCid(n,c,nPeptide.getScore()+cPeptide.getScore()));
//                    c++;
//                }
//                n++;
//            }
//            ncIdList.sort(Comparator.comparing(o->o.score, Comparator.reverseOrder()));
//
//            Peptide[] nPeptideArray = new Peptide[2];
//            nModPepsSet.toArray(nPeptideArray);
//            Peptide[] cPeptideArray = new Peptide[2];
//            cModPepsSet.toArray(cPeptideArray);
//            double lastScore = -1;
//            for (NCid ncId : ncIdList.subList(0,Math.min(2, ncIdList.size()))) {
//                Peptide nPartpeptide = nPeptideArray[ncId.nId];
//                Peptide cPartpeptide = cPeptideArray[ncId.cId];
//                if (nPartpeptide.getScore()+cPartpeptide.getScore() < lastScore) {
//                    continue;
//                }
//                String freePepSeq = nPartpeptide.getFreeSeq()+finderTag.getFreeAaString()+cPartpeptide.getFreeSeq();
//                Peptide fullPeptide = new Peptide(freePepSeq, nPartpeptide.isDecoy, massTool);// todo
//                fullPeptide.finderTag = finderTag;
//
//                PosMassMap fullPosMassMap;
//                if (nPartpeptide.hasVarPTM()){
//                    fullPosMassMap = nPartpeptide.getVarPTMs().clone(); // maybe replace PosMass at all
//                    fullPeptide.posVarPtmResMap.putAll(nPartpeptide.posVarPtmResMap);
//                } else {
//                    fullPosMassMap = new PosMassMap();
//                }
//                int idOfAa = nPartpeptide.length()-1;
//                for (char aaChar : finderTag.getPtmAaString().toCharArray()) {
//                    if (Character.isUpperCase(aaChar)) {
//                        idOfAa += 1;
//                    } else {
//                        fullPosMassMap.put(idOfAa, massTool.labelVarPtmMap.get(aaChar).mass);
//                        fullPeptide.posVarPtmResMap.put(idOfAa, massTool.labelVarPtmMap.get(aaChar));
//                    }
//                }
//                for (int pos : cPartpeptide.posVarPtmResMap.keySet()) { //copy the top 1 ptm pattern in n part // whhat if also choose the largeset priority
//                    fullPosMassMap.put(pos + finderTag.size() + nPartpeptide.length(), cPartpeptide.posVarPtmResMap.get(pos).mass); // copy the ptms from partModPepsUnsettled
//                    fullPeptide.posVarPtmResMap.put(pos + finderTag.size() + nPartpeptide.length(), cPartpeptide.posVarPtmResMap.get(pos));
//                }
//
//                if (! fullPosMassMap.isEmpty()) fullPeptide.setVarPTM(fullPosMassMap);
//                double calScore = massTool.buildVectorAndCalXCorr(fullPeptide.getIonMatrixNow(), 1, expProcessedPL, fullPeptide.matchedBions, fullPeptide.matchedYions);//todo decide the penalty
//                fullPeptide.setScore(calScore*(1 - fullPeptide.posVarPtmResMap.size() * 0.05));
//                fullPeptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, fullPeptide.getIonMatrixNow(), ms2Tolerance));
//                updatePeptideTreeSet(fullPeptide, peptideTreeSet, peptideInfoMap, protId, protSeq, tagPosInProt-nPartpeptide.length(), tagPosInProt+finderTag.size()+cPartpeptide.length()-1);
//            }
//        }
//    }


    private void updatePeptideTreeSet(Peptide newPeptide, TreeSet<Peptide> peptideTreeSet, Map<String, PeptideInfo> peptideInfoMap, String protId, String protSeq, int pepStartPos, int pepEndPos) {

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
                if ( ! peptideInfo.protIdSet.contains(protId)) { //if this pep with prot is not recorded, add this prot
                    peptideInfo.leftFlank = leftFlank;
                    peptideInfo.rightFlank = rightFlank;
                    peptideInfo.protIdSet.add(protId);
                    if (!protId.startsWith("DECOY_")) {
                        peptideInfo.isDecoy = false;
                    }
                }
            } else {
                peptideInfo = new PeptideInfo(freePepSeq, protId.startsWith("DECOY_"), leftFlank, rightFlank);
                peptideInfo.protIdSet.add(protId);
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
//         {
            partPeptide.setScore(score*(1-partPeptide.posVarPtmResMap.size()*0.05));
            partPeptide.setMatchedPeakNum(Score.getMatchedIonNum(plMap, 1, ionMatrix, ms2Tolerance));
            if (modPepsSet.size() < 2) { //max restore 2 patterns for one peptide  //todo make this avaliable to differentiate priority
                modPepsSet.add(partPeptide);
            } else if (modPepsSet.last().compareTo(partPeptide) < 0) {
                modPepsSet.pollLast();
                modPepsSet.add(partPeptide);
            }
//        }
    }
    public class NCid {
        public NCid(int nId, int cId, double score) {
            this.nId = nId;
            this.cId = cId;
            this.score = score;
        }
        public int nId;
        public int cId;
        public double score;
//        private double getScore(){
//            return score;
//        }
    }
//    private void updateCandiList(int scanNum, String protId, int tagPosInProt, ExpTag finderTag, double ms1TolAbs, Map<String, List<Peptide>> resPeptideListMap
//            , Map<String, PeptideInfo> peptideInfoMap, SparseVector expProcessedPL, TreeMap<Double, Double> plMap, GRBEnv env) throws CloneNotSupportedException, GRBException {
//
//        double tagCMass = finderTag.getTailLocation() + massTool.H2O-MassTool.PROTON; // +massTool.H2O-MassTool.PROTON  is to mimic the mass of the real neutral precursor mass
//        //tagCMass is correct even the tag is fuzzy. because FM starts searching from C
//        String protSeq = buildIndex.protSeqMap.get(protId);
//        Map<Integer, Double> cPoscMassMap = new HashMap<>();
//        if ( (tagPosInProt+finderTag.size() == protSeq.length() && finderTag.isNorC != C_TAG)  // not C tag but found at prot Cterm, impossible
//                || (tagPosInProt == 0 && finderTag.isNorC != N_TAG) // not N tag but found at prot Nterm, impossible
//                || tagPosInProt < 0
//                || tagPosInProt+finderTag.size() > protSeq.length()
//        ) {
//            return;
//        }
//
//        int remainMC = massTool.missedCleavage - getNumOfMissCleavSite(finderTag.getFreeAaString());
//
//        List<Map<Integer, Integer>> cResPosPtmIdList = new ArrayList<>();
//        if (finderTag.isNorC == C_TAG || Math.abs(tagCMass - precursorMass) < ms1TolAbs){ // C tag, tagCMass must be exactly settle. Otherwise no valid cPos will be, then no valid n pos.
//            cPoscMassMap.put(tagPosInProt+finderTag.size()-1, tagCMass); // amino acid at cPos is also counted
//            System.out.println("lsz No C section");
//        } else { // only when oriTag is not C oriTag can it extend to c
//            List<Integer> cPoses = new ArrayList<>();
//            List<Integer> krPoses = new ArrayList<>();
//            int mcInC = 0;
//            for (int i = tagPosInProt+finderTag.size(); i < protSeq.length(); i++) {  //
//                //the if-block is in order, dont change the order.
//                if (isX(protSeq.charAt(i))) {
//                    break;
//                }
//                tagCMass += massTool.getMassTable().get(protSeq.charAt(i)); //even when tag is fuzzy, tagcmass wont be disturbed
//                if (tagCMass > precursorMass+maxPtmMass) break; // the total max minus ptm is -250. although the max minus ptm for single ptm is -156
//
//                if (isKR(protSeq.charAt(i))) {
//                    mcInC++; //missed cleavage
//                }
//                if (tagCMass >= precursorMass-maxPtmMass) { // the total max plus ptm is 600. though the max plus ptm for single ptm is 527
//                    cPoses.add(i);
//                    cPoscMassMap.put(i, tagCMass);
//                    if (isKR(protSeq.charAt(i))) {
//                        krPoses.add(i);
//                    }
//                }
//                if (mcInC > remainMC) { // if current num of KR is max, dont need to extend to c because it is impossible to take one more KR
//                    cPoses.add(i);
//                    cPoscMassMap.put(i, tagCMass);
//                    break;         // stop extend to C term
//                }
//            }// finish collect all possible cpos
//            if (lszDebugScanNum.contains(scanNum) && finderTag.getFreeAaString().contentEquals("VDLDAPDVSL")){
//                int a = 1;
//            }
//            int maxCPos = Collections.max(cPoses);
//            int minCPos = Collections.min(cPoses);
//            byte isProtNorC_Term = NON_TERM_PROT;
//            if (maxCPos == protSeq.length()-1) {
//                isProtNorC_Term = C_TERM_PROT;
//            }
//            Map<Integer, Integer> posYIdMap = new HashMap<>();
//            Map<Integer, Integer> yIdMaxPosMap = new HashMap<>();
//            if (cTermSpecific) {
////                if ( ! krPoses.isEmpty()) { // if no kr poses but cTermSpecific, all poses will be trimed , no need to continue
//                Set<Integer> posIdToTrim = new HashSet<>();
//                for (int k = cPoses.size() - 1; k >= 0; k--) {
//                    if (!isKR(protSeq.charAt(cPoses.get(k)))) {
//                        posIdToTrim.add(cPoses.get(k));
//                    } else {
//                        break; //break at the first time we see a KR
//                    }
//                }
//                cPoses.removeAll(posIdToTrim); // to ensure c specific
//                if ( ! cPoses.isEmpty() && ! krPoses.isEmpty()) {
//                    System.out.println("lsz cPoses empty");
//                    maxCPos = Collections.max(cPoses);
//                    minCPos = Collections.min(cPoses);
//                    if (krPoses.isEmpty()) {
//                        System.out.println("lsz krPoses empty");
//                    }
//                    minCPos = krPoses.get(0) + 1;
//                    krPoses.remove(0);
//                    int yId = 0;
//                    for (int pos = minCPos; pos < maxCPos + 1; pos++) {
////                        int relPos = pos-(tagPosInProt+finderTag.size());
//                        posYIdMap.put(pos, yId); // using rel pos in part seq i.e. start from 0 not the prot pos
//                        int tmpMaxPos = yIdMaxPosMap.getOrDefault(yId, -99);
//                        yIdMaxPosMap.put(yId, pos > tmpMaxPos ? pos : tmpMaxPos);
//                        if (krPoses.contains(pos)) {
//                            yId++;
//                        }
//                    }
//                }
//            } else {
//                int yId = 0;
//                for (int pos = minCPos; pos < maxCPos+1; pos++) {
////                    int relPos = pos-(tagPosInProt+finderTag.size());
//                    posYIdMap.put(pos, yId); // using rel pos in part seq i.e. start from 0 not the prot pos
//                    int tmpMaxPos = yIdMaxPosMap.getOrDefault(yId, -1);
//                    yIdMaxPosMap.put(yId, pos > tmpMaxPos ? pos : tmpMaxPos);
//                    yId++;
//                }
//            }
//            String cPartSeq = protSeq.substring(tagPosInProt+finderTag.size(), maxCPos+1);
////                double cCutMass = finderTag.getTailLocation() - MassTool.PROTON;
//            double flexiableMass = precursorMass - (finderTag.getTailLocation() + massTool.H2O-MassTool.PROTON) - massTool.calResidueMass(protSeq.substring(tagPosInProt+finderTag.size(), minCPos));
//            cResPosPtmIdList = inferPTM.findBestPtmInC_Spec(scanNum, env, cPartSeq, (tagPosInProt+finderTag.size()), ms1TolAbs, flexiableMass, posYIdMap, yIdMaxPosMap, isProtNorC_Term);
//        } // end settle C section
//
//
////        protSeq = "PEPTLDEKLVFSGNLFQHQEDSKKLQDELQNMKEEMARLVSGKDYNVTANSKLGSPSGKTGTCRMSFVKGWGAEYRLLEQKTQESQKGFCFLTYTDEEPVKKTQAYQDQKPGTSGLRLLNEPTAAALAYGLDKKLLDLKTGTVKAANMLQQSGSKNTGAKAGPNASLLSLKSDKGADFLVTEVENGGSLGSKKDASTLQSQKAEGTGDAKGTVPDDAVEALADSLGKKSKHLQEQLNELKLLEEEDSKLKLVEKYGYTHLSAGELLRVNVEAPDVNLEGLGGKLKSEASSEFAKKDTPTSAGPNSFNKGKGTLSAPGKVVTAAAQAKAVVKTDLQLALPSGCYGRLSALQGKLSKSLFFGGKGAPLQRATVASSTQKFQDLGVKASNAAVKLAESKDHLSDKESKAAQASDLEKLHLDEKLGNDFHTNKRDQEALMKSVKVNLGQGSHPQKVKSQNTDMVQKSVSKVGLDTPDLDLHGPEGKLKLLTWDVKDTLLRVDLDVPDVNLEGPDAKLKFLQEFYQDDELGKKDALHFYNKSLAEHRVYLASSSGSTALKKGWLKSNVSDAVAQSTRLFQSDTNAMLGKKDLGKFQVATDALKSGNLTEDDKHNNAKLSSNFSSLLAEKLRAMEAASSLKSQAATKKYEELDNAPEERTLQKQSVVYGGKGETASKLQSELSRTKMSAAQALKSPSDSGYSYETLGKTTKLPQEKCLLQTDVKLTQDLMAKVRLEEVTGKLQVARLQTKAGEATVKLDLPAENSNSETLLLTGKRVNALKNLQVKHELQANCYEEVKDRLTASLCDLKSRSSFAYKDQNENRFSMPGFKGEGPDVDVSLPKLMDQNLKCLSAAEEKAGQKLLDVNHYAKLDHLAEKFRALASLDSKLNQAKVEDCKMESTETEERSSGPGGQNVNKVNSKDSLNFLDHLNGKKYLELYVQKVNPSRTGGTQTDLFTCGKCKSLKLNTAEDAKSPDEAYALAKKSDSGKPYYYNSQTKVDVEVPDVSLEGPEGKLKVFLGNLNTLVVKKKNLSQELNMEQKLTSLNVKYNNDKLLLPGELAKHAVSEGTKSGFGKFERLDEKENLSAKVQGGALEDSQLVAGVAFKKPGGEGPQCEKTTDVKANLEKAQAELVGTADEATRPANEKATDDYHYEKSGNFSAAMKDLSGKAALNQKLLETGERTEVLSPNSKVESKPSVSTKLLSKGSLGAQKLANTCFNELEKGTLTVSAQELKDNRVFSLLSSEKELKFSMPGFKGEGPDVDVNLPKDDSFLGKLGGTLARVGDDLAKATGDWKDNSTMGYMAAKKDLSTVEALQNLKNLKFDDTNPEKEEAKDAEAAEATAEGALKAEKYEELVKEVSTYLKMESEEGKEARALLQAKASGKASDTSSETVFGKRMAQELKYGDLVERDFDTALKHYDKLDEMPEAAVKSTANKVGLQVVAVKAPGFGDNRDKALPLEAEVAPVKASEALLKQLKLQLELDQKKDYSAPVNFLSAGLKKFGTLNLVHPKLGSYTKLQAAYAGDKADDLQKTPGQNAQKWLPARLFVGGLKEDTEEHHLRLQEKVESAQSEQKLEKTLDDLEEKDSSPTTNSKLSALQRSLGEPAAGEGSAKRLLSSLEQKEENKTTEETSKPKLHNAENLQPGEQKYEYKYEQEHAALQDKLFQVAKSEDTLSKMNDFMRLTESVAETAQTLKKVSPASSVDSNLPSSQGYKKTGASWTDNLMAQKCSKDESFLGKLGGTLARDYSSGFGGKYGVQADRTFGKAAGPSLSHTSGGTQSKVETTVTSLKTKTTGFYSGFSEVAEKRLTLSKHQNVQLPRVGTASVLQPVKKVLTANSNPSSPSAAKRANEKTESSSAQQVAVSRLSSAMSAAKALCDHVRLQDVSGQLSSSKKLTGKNQVTATKTQDQLSNLKYHEEFEKDVNAALAALKTKSPPLPEHQKLPCNSAEPKLAEEENKAKSSVACKWNLAEAQQKPSVGSQSNQAGQGKRTEDEVLTSKGDAWAKAQQHWGSGVGVKKKFGVLSDNFKVDLDAPDVSLEGPDAKLKLLFSNTAAQKLRVLDSGAPLKLPVGPETLGRSNEEGSEEKGPEVRSGDWKGYTGKLTMLNTVSKLRFSAYLKNSNPALNDNLEKDEPSVAAMVYPFTGDHKQKMGSKSPGNTSQPPAFFSKVTLEVGKVLQQGRVTAAAMAGNKSTPRDTMGLADKTENTLERGAEKTEVAEPRSMSTEGLMKFVDSKVGEFSGANKEKASALSLKLGSSKGSLEKGSPEDKGEVGAGAGPGAQAGPSAKRHAEMVHTGLKLERSLTHSTVGSKGEKGGNLGAKWTLDLKLGNSYFKEEKSDDESSQKDLKDLVAALQHNYKMSAFKLNNFSADLKDSKLEQVEKELLRVPTTEKPTVTVNFRVGVEVPDVNLEGPEGKLKLLLVSTTPYSEKDTKFELGKLMELHGEGSSSGKSTGLELETPSLVPVKKLQQELQTAKKPSEGPLQSVQVFGRSLEEEGAAHSGKRSKLLFSNTAAQKTPNLYLYSKKDTGEKLTVAENEAETKLLEEELQAPTSSKRLGVATVAKGVCGACKMMSKPQTSGAYVLNKGGPGSTLSFVGKRAANDAGYFNDEMAPLEVKTKHAVSEGTKAVTKTVSKVDDFLANEAKSGKFTQQDLDEAKGFALVGVGSEASSKKMQQLKEQYRGEVQVSDKERTLDELLQEKKLTEKELAEAASKLHVQSSSDSSDEPAEKRSQPATKTEDYGMGPGRLPLANTEKYMADKPASGAGAGAGAGKRKPEPTLDE";
//        List<Map<Integer, Integer>> nResPosPtmIdList = new ArrayList<>();
//        double nDeltaMass = finderTag.getHeadLocation() - MassTool.PROTON; //n term delta mass should be independent of cDeltaMass
//        if (finderTag.isNorC == N_TAG || Math.abs(finderTag.getHeadLocation() - MassTool.PROTON) < ms1TolAbs) { //fixme, should use MS1 tolerance not MS2 tolerance
//            System.out.println("lsz No N section"); //todo
//        } else {
//            List<Integer> nPoses = new ArrayList<>();
//            List<Integer> krPosesInFlex = new ArrayList<>();
//            int mcInN = 0;
//            for (int nPos = tagPosInProt-1; nPos >= 0; nPos--) {
//                if (isX(protSeq.charAt(nPos))) break;
//                if (nPos < tagPosInProt) {
//                    nDeltaMass -= massTool.getMassTable().get(protSeq.charAt(nPos));
//                }
//                if (nDeltaMass < minPtmMass) break;
//                if (nDeltaMass < maxPtmMass) {
//                    nPoses.add(nPos);
//                    if (isKR(protSeq.charAt(nPos))) {
//                        krPosesInFlex.add(nPos); // only count the kr poses in the flexible zone
//                    }
//                }
//                if (isKR(protSeq.charAt(nPos))) {
//                    mcInN++; //missed cleavage
//                }
//                if (mcInN > remainMC) { // if current num of KR is max, dont need to extend to c because it is impossible to take one more KR
//                    break;         // diff NC
//                }
//                double nCutMass = precursorMass + MassTool.PROTON - finderTag.getHeadLocation() - massTool.H2O; //                double nCutMass = finderTag.getTailLocation() - MassTool.PROTON;
//            }
//            if (lszDebugScanNum.contains(scanNum) && finderTag.getFreeAaString().contentEquals("VDLDAPDVSL")){
//                int a = 1;
//            }
//            int maxNPos = Collections.max(nPoses);
//            int minNPos = Collections.min(nPoses);
//            Map<Integer, Integer> posYIdMap = new HashMap<>();
//            byte isProtNorC_Term = NON_TERM_PROT;
//            if (minNPos == 1 && protSeq.charAt(0) == 'M') {
//                isProtNorC_Term = N_TERM_PROT;
//            }
//            Map<Integer, Integer> yIdMinPosMap = new HashMap<>();
//            if (nTermSpecific) {
//                Set<Integer> posIdToTrim = new HashSet<>();
//                if (nPoses.get(nPoses.size()-1) > 1 && ! isKR(protSeq.charAt(nPoses.get(nPoses.size()-1)-1))){  //if min pos is at prot head, dont trim; else if min pos-1 is KR dont trim
//                    for (int k = nPoses.size()-1; k >= 0; k--){
//                        posIdToTrim.add(nPoses.get(k));
//                        if (isKR(protSeq.charAt(nPoses.get(k)))) {
//                            break; //break at the first time we see a KR
//                        }
//                    }
//                }
//                nPoses.removeAll(posIdToTrim); // to ensure c specific
//                if ( ! nPoses.isEmpty()) {
//                    System.out.println("lsz nPoses empty");
//                    maxNPos = Collections.max(nPoses);
//                    minNPos = Collections.min(nPoses);
//                    if (krPosesInFlex.isEmpty()) {
//                        System.out.println(scanNum+"lsz krPoses empty");
//                    }
////                    maxNPos = krPosesInFlex.get(0);
//                    maxNPos = krPosesInFlex.isEmpty() ? minNPos - 1 : krPosesInFlex.get(0);
//
//
////                    krPosesInFlex.remove(0);
//                    int yId = 0;
//                    for (int pos = maxNPos; pos > minNPos - 1; pos--) {
////                        int relPos = pos - minNPos;
//                        posYIdMap.put(pos, yId);
//                        if (krPosesInFlex.contains(pos)) {
//                            int tmpMinPos = yIdMinPosMap.getOrDefault(yId, 999999);
//                            yIdMinPosMap.put(yId, pos < tmpMinPos ? pos : tmpMinPos);
//                            yId++;
//                        }
//                    }
//                }
//            } else {
//                maxNPos = Collections.max(nPoses);
//                minNPos = Collections.min(nPoses);
//                int yId = 0;
//                for (int pos = maxNPos; pos > minNPos-1; pos--) {
////                    int relPos = pos - minNPos;
//                    posYIdMap.put(pos, yId);
//                    int tmpMinPos = yIdMinPosMap.getOrDefault(yId, 999999);
//                    yIdMinPosMap.put(yId, pos < tmpMinPos ? pos : tmpMinPos);
//                    yId++;
//                }
//            }
//            String nPartSeq = protSeq.substring(minNPos, tagPosInProt);
////                double nCutMass = precursorMass+MassTool.PROTON - finderTag.getHeadLocation() - massTool.H2O;
//            double flexiableMass = finderTag.getHeadLocation() - MassTool.PROTON - massTool.calResidueMass(protSeq.substring(maxNPos+1, tagPosInProt));
////            Set<Integer> fixModIdxes = inferPTM.getFixModIdxes(nPartSeq);
////            Map<Integer, VarPtm[]> posAllVarPtmMap = inferPTM.getIdxVarModMapNew(nPartSeq, fixModIdxes, N_PART, isProtNorC_Term);
//            nResPosPtmIdList = inferPTM.findBestPtmInN_Spec(scanNum, env, nPartSeq, minNPos, ms1TolAbs, flexiableMass, posYIdMap, yIdMinPosMap,isProtNorC_Term);
//        }// end settle N section
//        return;
//    }

//    private void updateCandiList(int scanNum, String protId, int tagPosInProt, ExpTag finderTag, double ms1TolAbs, Map<String, List<Peptide>> resPeptideListMap
//            , Map<String, PeptideInfo> peptideInfoMap, SparseVector expProcessedPL, TreeMap<Double, Double> plMap, GRBEnv env) throws CloneNotSupportedException, GRBException {
//true
//        double tagCMass = finderTag.getTailLocation() + massTool.H2O-MassTool.PROTON; // +massTool.H2O-MassTool.PROTON  is to mimic the mass of the real neutral precursor mass
//        //tagCMass is correct even the tag is fuzzy. because FM starts searching from C
//        String protSeq = buildIndex.protSeqMap.get(protId);
//        Map<Integer, Double> cPoscMassMap = new HashMap<>();
//        if ( (tagPosInProt+finderTag.size() == protSeq.length() && finderTag.isNorC != C_TAG)  // not C tag but found at prot Cterm, impossible
//                || (tagPosInProt == 0 && finderTag.isNorC != N_TAG) // not N tag but found at prot Cterm, impossible
//                || tagPosInProt < 0
//                || tagPosInProt+finderTag.size() > protSeq.length()
//        ) {
//            return;
//        }
//
//        int missCleav = getNumOfMissCleavSite(finderTag.getFreeAaString());
//        if (finderTag.isNorC == C_TAG || Math.abs(tagCMass - precursorMass) < ms1TolAbs){ // C tag, the end of tag must be KR and the tagCMass must be exactly settle. Otherwise no valid cPos will be, then no valid n pos.
//            if (isKR(protSeq.charAt(tagPosInProt+finderTag.size()-1)) && Math.abs(tagCMass - precursorMass) < ms1TolAbs) { //fixme, should use MS1 tolerance not MS2 tolerance
//                cPoscMassMap.put(tagPosInProt+finderTag.size()-1, tagCMass); // amino acid at cPos is also counted
//            }
//        } else { // only when oriTag is not C oriTag can it extend to c
//            for (int i = tagPosInProt+finderTag.size(); i < protSeq.length(); i++) {  //
//                //the if-block is in order, dont change the order.
//                if (isX(protSeq.charAt(i))) break;
//                tagCMass += massTool.getMassTable().get(protSeq.charAt(i)); //even when tag is fuzzy, tagcmass wont be disturbed
//                if (tagCMass > precursorMass+maxPtmMass) break; // the total max minus ptm is -250. although the max minus ptm for single ptm is -156
//
//                if (isKR(protSeq.charAt(i))) {
//                    missCleav++; //current num of missed cleavage
//                }
//                if (tagCMass >= precursorMass-maxPtmMass && isKR(protSeq.charAt(i))) { // the total max plus ptm is 600. though the max plus ptm for single ptm is 527
//                    cPoscMassMap.put(i, tagCMass);
//                }
//                if (missCleav > massTool.missedCleavage  && !isPtmSimuTest) { // if current num of KR is max, dont need to extend to c because it is impossible to take one more KR
//                    break;         // stop extend to C term
//                }
//            }
//        }
//
//        int n_extension = 0;
//        for (int cPos : cPoscMassMap.keySet()) {
//            //check and find valid solution for cDeltaMass
//            double cDeltaMass = precursorMass - cPoscMassMap.get(cPos);
//            List<TreeMap<Integer, VarPtm>> cPosVarPtmResMap_List = new LinkedList<>();
//            char rightFlank;
//            byte isProtNorC_Term = NON_TERM_PROT;
//            if (cPos == protSeq.length()-1) {
//                rightFlank = '-';
//                isProtNorC_Term = C_TERM_PROT;
//            } else {
//                rightFlank = protSeq.charAt(cPos+1);
//            }
//            String cPartSeq = protSeq.substring(tagPosInProt+finderTag.size(), cPos+1);
//            double cCutMass = finderTag.getTailLocation() - MassTool.PROTON;
//            boolean isCTermFree = false;
//            if (Math.abs(cDeltaMass) > 0.1) { // not ptm free in c part
//                Set<Integer> fixModIdxes = inferPTM.getFixModIdxes(cPartSeq);
//
//                Map<Integer, Set<Double>> oneTimeMassGroups = new HashMap<>(cPartSeq.length());
//                List<Pair<Integer, Set<Double>>> multiTimeMassGroups = new ArrayList<>(cPartSeq.length());
//                Map<Integer, Map<Double, VarPtm>> pos_MassVarPtm_Map = new HashMap<>(cPartSeq.length(), 1);
//                Map<Double, List<Integer>> allMassAllPosesMapC = inferPTM.getMassPosInfoMap(scanNum, cPartSeq, fixModIdxes, C_PART, isProtNorC_Term, oneTimeMassGroups, multiTimeMassGroups, pos_MassVarPtm_Map);
//                ModPepPool cPartPepsWithPtm = inferPTM.settlePtmOnSide(scanNum, expProcessedPL, plMap, precursorMass, cPartSeq, false,
//                        allMassAllPosesMapC, cCutMass, cDeltaMass, precursorCharge, C_PART,ms1TolAbs, oneTimeMassGroups, multiTimeMassGroups, pos_MassVarPtm_Map,env);
//
//                if (cPartPepsWithPtm.peptideTreeSet.isEmpty()) {
//                    continue;
//                }
//                double topScore = cPartPepsWithPtm.getTopPepPtn().getScore();
//                boolean shouldSkip = true;
//                for (Peptide peptide : cPartPepsWithPtm.peptideTreeSet){
//                    if (peptide.getScore() >= 0.5 *topScore){
////                        if (!peptide.posVarPtmResMap.isEmpty()){
//                        cPosVarPtmResMap_List.add(peptide.posVarPtmResMap);
//                        shouldSkip = false;
////                        }
//                    }
//                }
//                if (shouldSkip) {
//                    continue;  // c part should has PTM but unsettled.
//                }
//            } else {
//                isCTermFree = true;
//            }
//
//
//            //extend and settle N part
//            double nDeltaMass = finderTag.getHeadLocation() -MassTool.PROTON; //n term delta mass should be independent of cDeltaMass
//            missCleav = getNumOfMissCleavSite(protSeq.substring(tagPosInProt, cPos)) ; //dont contain the c-term K, it does not count as miss cleav
//
//            int min_nPos = (finderTag.isNorC == N_TAG && Math.abs(finderTag.getHeadLocation()-MassTool.PROTON) <= 0.02 ) ? tagPosInProt : 0;  //if it is from N tag and (when fuzzy) it is still N tag
//            int max_nPos = (finderTag.isNorC == N_TAG && Math.abs(finderTag.getHeadLocation()-MassTool.PROTON) <= 0.02 ) ? tagPosInProt : tagPosInProt-1;
//            for (int nPos = max_nPos; nPos >= min_nPos; nPos--) {
//                if (cPos+1-nPos < minPepLen){
//                    continue;
//                }
//                else if ( cPos+1-nPos > maxPepLen) {
//                    break;
//                }
//                if (isX(protSeq.charAt(nPos))) break;
//                if (nPos < tagPosInProt) {
//                    nDeltaMass -= massTool.getMassTable().get(protSeq.charAt(nPos));
//                }
//                if (nDeltaMass < minPtmMass) break;
//                if (nDeltaMass > maxPtmMass)
//                    continue;// the total max plus ptm is 600. though the max plus ptm for single ptm is 527
//
//                if (nTermSpecific) {
//                    if (nPos != 0 && !isKR(protSeq.charAt(nPos - 1))) {// n term must be specific
//                        continue;
//                    }
//                }
//
////                if (isKR(protSeq.charAt(nPos))) {  //this is not working because the continue above
////                    missCleav++; //current num of missed cleavage
////                }
//                missCleav = getNumOfMissCleavSite(protSeq.substring(nPos, cPos));
//                if (missCleav > massTool.missedCleavage && !isPtmSimuTest) {
//                    break;         // stop extend to n term
//                }
//                if (cPos + 1 - nPos < 6) { // min length of pep
//                    continue;
//                }
//
//                // check N part
//                String nPartSeq = protSeq.substring(nPos, tagPosInProt);// should not plus 1
//                char leftFlank;
//                isProtNorC_Term = NON_TERM_PROT;
//                if (nPos == 0 || (nPos == 1 && protSeq.charAt(0) == 'M')) {
//                    leftFlank = '-';
//                    isProtNorC_Term = N_TERM_PROT;
//                } else {
//                    leftFlank = protSeq.charAt(nPos - 1);
//                }
//                double nCutMass = precursorMass+MassTool.PROTON - finderTag.getHeadLocation() - massTool.H2O; //                double nCutMass = finderTag.getTailLocation() - MassTool.PROTON;
//
//                List<TreeMap<Integer, VarPtm>> nPosVarPtmResMap_List = new LinkedList<>();
//
//                boolean isNTermFree = false;
//                if (Math.abs(nDeltaMass) > 0.1) {
//                    Set<Integer> fixModIdxes = inferPTM.getFixModIdxes(nPartSeq);
//                    Map<Integer, Set<Double>> oneTimeMassGroups = new HashMap<>(nPartSeq.length());
//                    List<Pair<Integer, Set<Double>>> multiTimeMassGroups = new ArrayList<>(nPartSeq.length());
//                    Map<Integer, Map<Double, VarPtm>> pos_MassVarPtm_Map = new HashMap<>(nPartSeq.length(), 1);
//                    Map<Double, List<Integer>> allMassAllPosesMapN = inferPTM.getMassPosInfoMap(scanNum, nPartSeq, fixModIdxes, N_PART, isProtNorC_Term, oneTimeMassGroups, multiTimeMassGroups, pos_MassVarPtm_Map); //todo no need to generate var mod list for aa again and again, make it stored.
//                    ModPepPool nPartpepsWithPtm = inferPTM.settlePtmOnSide(scanNum, expProcessedPL, plMap, precursorMass, nPartSeq, false,
//                            allMassAllPosesMapN, nCutMass, nDeltaMass, precursorCharge, N_PART, ms1TolAbs, oneTimeMassGroups, multiTimeMassGroups, pos_MassVarPtm_Map, env);
//
//                    if (nPartpepsWithPtm.peptideTreeSet.isEmpty()) {
//                        continue; // how could it be return!!!
//                    }
//                    double topScore = nPartpepsWithPtm.getTopPepPtn().getScore();
//                    boolean shouldSkip = true;
//                    for (Peptide peptide : nPartpepsWithPtm.peptideTreeSet) {
//                        if (peptide.getScore() >= 0.5 * topScore) {
////                            if (!peptide.posVarPtmResMap.isEmpty()) {
//                            nPosVarPtmResMap_List.add(peptide.posVarPtmResMap);
//                            shouldSkip = false;
////                            }
//                        }
//                    }
//                    if (shouldSkip) {
//                        continue;  // n part should has PTM but unsettled.
//                    }
//                } else {
//                    isNTermFree = true;
//                }
//                n_extension++;
//                int tagPosInPep = tagPosInProt - nPos;
//                // store this peptide
//                // to combine the top 10 cPartSeq and top 10 nPartSeq.
//                String freePepSeq = protSeq.substring(nPos, cPos + 1);
//                StringBuilder ptmPepSeqSB = new StringBuilder(freePepSeq);
//                ptmPepSeqSB.replace(tagPosInPep, tagPosInPep + finderTag.size(), finderTag.getPtmAaString());
//                List<PosMassMap> posMassMap_List = new LinkedList<>();
//                List<TreeMap<Integer, VarPtm>> posVarPtmResMap_List = new LinkedList<>();
//                if (!cPosVarPtmResMap_List.isEmpty() && !nPosVarPtmResMap_List.isEmpty()) {
//                    for (TreeMap<Integer, VarPtm> cPosVarPtmResMap : cPosVarPtmResMap_List){
//                        for (TreeMap<Integer, VarPtm> nPosVarPtmResMap : nPosVarPtmResMap_List) {
//                            PosMassMap fullPosMassMap = new PosMassMap(freePepSeq.length());
//                            TreeMap<Integer, VarPtm> posVarPtmResMap = new TreeMap<>();
//                            for (int pos : nPosVarPtmResMap.keySet()) { //copy the top 1 ptm pattern in n part // whhat if also choose the largeset priority
//                                fullPosMassMap.put(pos, nPosVarPtmResMap.get(pos).mass); // copy the ptms from partModPepsUnsettled
//                                posVarPtmResMap.put(pos, nPosVarPtmResMap.get(pos));
//                            }
//                            for (int pos : cPosVarPtmResMap.keySet()) { //copy the top 1 ptm pattern in n part // whhat if also choose the largeset priority
//                                fullPosMassMap.put(pos + tagPosInPep + finderTag.size(), cPosVarPtmResMap.get(pos).mass); // copy the ptms from partModPepsUnsettled
//                                posVarPtmResMap.put(pos + tagPosInPep + finderTag.size(), cPosVarPtmResMap.get(pos));
//                            }
//                            int idOfAa = -1;
//                            for (char aaChar : finderTag.getPtmAaString().toCharArray()) {
//                                if (Character.isUpperCase(aaChar)) {
//                                    idOfAa += 1;
//                                } else {
//                                    fullPosMassMap.put(tagPosInPep + idOfAa, massTool.labelVarPtmMap.get(aaChar).mass);
//                                    posVarPtmResMap.put(tagPosInPep + idOfAa, massTool.labelVarPtmMap.get(aaChar));
//                                }
//                            }
//                            posVarPtmResMap_List.add(posVarPtmResMap);
//                            posMassMap_List.add(fullPosMassMap);
//                        }
//                    }
//                } else if (cPosVarPtmResMap_List.isEmpty() && !nPosVarPtmResMap_List.isEmpty()) {
//                    for (TreeMap<Integer, VarPtm> nPosVarPtmResMap : nPosVarPtmResMap_List) {
//                        PosMassMap fullPosMassMap = new PosMassMap(freePepSeq.length());
//                        TreeMap<Integer, VarPtm> posVarPtmResMap = new TreeMap<>();
//                        for (int pos : nPosVarPtmResMap.keySet()) { //copy the top 1 ptm pattern in n part // whhat if also choose the largeset priority
//                            fullPosMassMap.put(pos, nPosVarPtmResMap.get(pos).mass); // copy the ptms from partModPepsUnsettled
//                            posVarPtmResMap.put(pos, nPosVarPtmResMap.get(pos));
//                        }
//                        int idOfAa = -1;
//                        for (char aaChar : finderTag.getPtmAaString().toCharArray()) {
//                            if (Character.isUpperCase(aaChar)) {
//                                idOfAa += 1;
//                            } else {
//                                fullPosMassMap.put(tagPosInPep + idOfAa, massTool.labelVarPtmMap.get(aaChar).mass);
//                                posVarPtmResMap.put(tagPosInPep + idOfAa, massTool.labelVarPtmMap.get(aaChar));
//                            }
//                        }
//                        posVarPtmResMap_List.add(posVarPtmResMap);
//                        posMassMap_List.add(fullPosMassMap);
//                    }
//                }else if (!cPosVarPtmResMap_List.isEmpty() && nPosVarPtmResMap_List.isEmpty()){
//                    for (TreeMap<Integer, VarPtm> cPosVarPtmResMap : cPosVarPtmResMap_List){
//                        PosMassMap fullPosMassMap = new PosMassMap(freePepSeq.length());
//                        TreeMap<Integer, VarPtm> posVarPtmResMap = new TreeMap<>();
//                        for (int pos : cPosVarPtmResMap.keySet()) { //copy the top 1 ptm pattern in n part // whhat if also choose the largeset priority
//                            fullPosMassMap.put(pos + tagPosInPep + finderTag.size(), cPosVarPtmResMap.get(pos).mass); // copy the ptms from partModPepsUnsettled
//                            posVarPtmResMap.put(pos + tagPosInPep + finderTag.size(), cPosVarPtmResMap.get(pos));
//                        }
//                        int idOfAa = -1;
//                        for (char aaChar : finderTag.getPtmAaString().toCharArray()) {
//                            if (Character.isUpperCase(aaChar)) {
//                                idOfAa += 1;
//                            } else {
//                                fullPosMassMap.put(tagPosInPep + idOfAa, massTool.labelVarPtmMap.get(aaChar).mass);
//                                posVarPtmResMap.put(tagPosInPep + idOfAa, massTool.labelVarPtmMap.get(aaChar));
//                            }
//                        }
//                        posVarPtmResMap_List.add(posVarPtmResMap);
//                        posMassMap_List.add(fullPosMassMap);
//                    }
//                }
//
//                if (posVarPtmResMap_List.isEmpty()){
//                    Peptide peptide = new Peptide(freePepSeq, true, massTool);// these paras are dummy answer will be deleted
//                    peptide.tagPosInPep = tagPosInPep;
//                    peptide.ptmSeq = ptmPepSeqSB.toString();
//                    peptide.finderTag = finderTag;
//                    peptide.cDeltaMass = cDeltaMass;
//                    peptide.nDeltaMass = nDeltaMass;
//
//                    peptide.setScore(massTool.buildVectorAndCalXCorr(peptide.getIonMatrixNow(), 1, expProcessedPL, peptide.matchedBions, peptide.matchedYions));//todo decide the penalty
//
//                    List<Peptide> peptideList = resPeptideListMap.get(freePepSeq);
//                    if (peptideList == null) {
//                        peptideList = new ArrayList<>();
//                        peptideList.add(peptide);
//                        resPeptideListMap.put(freePepSeq, peptideList);
//                    } else {
//                        peptideList.add(peptide);
//                    }
//                }else {
//                    for (int i = 0; i < posVarPtmResMap_List.size(); i++) {
//                        Peptide peptide = new Peptide(freePepSeq, true, massTool);// these paras are dummy answer will be deleted
//                        peptide.tagPosInPep = tagPosInPep;
//                        peptide.ptmSeq = ptmPepSeqSB.toString();
//                        peptide.finderTag = finderTag;
//                        peptide.cDeltaMass = cDeltaMass;
//                        peptide.nDeltaMass = nDeltaMass;
//
//                        PosMassMap fullPosMassMap = posMassMap_List.get(i);
//                        TreeMap<Integer, VarPtm> posVarPtmResMap = posVarPtmResMap_List.get(i);
//                        if (!fullPosMassMap.isEmpty()) { // has ptm in n , c or on the tag, it is a ptm-containing peptide
//                            peptide.setVarPTM(fullPosMassMap);
//                            peptide.posVarPtmResMap = posVarPtmResMap;
//                        }
//                        peptide.setScore(massTool.buildVectorAndCalXCorr(peptide.getIonMatrixNow(), 1, expProcessedPL, peptide.matchedBions, peptide.matchedYions));//todo decide the penalty
//
//                        List<Peptide> peptideList = resPeptideListMap.get(freePepSeq);
//                        if (peptideList == null) {
//                            peptideList = new ArrayList<>();
//                            peptideList.add(peptide);
//                            resPeptideListMap.put(freePepSeq, peptideList);
//                        } else {
//                            peptideList.add(peptide);
//                        }
//                    }
//                }
//                //update peptideInfoMapForRef,  just the flanks and proteins, these are standard info should be independent of peptide candidates.
//                if (peptideInfoMap.containsKey(freePepSeq)) {
//                    PeptideInfo peptideInfo = peptideInfoMap.get(freePepSeq);
//                    if (!peptideInfo.protIdSet.contains(protId)) { //if this pep with prot is not recorded, add this prot
//                        if (peptideInfo.leftFlank != '-' && peptideInfo.rightFlank != '-') {
//                            if (rightFlank == '-' || leftFlank == '-') {
//                                peptideInfo.leftFlank = leftFlank;
//                                peptideInfo.rightFlank = rightFlank;
//                            }
//                        }
//                        peptideInfo.protIdSet.add(protId);
//                        if (!protId.startsWith("DECOY_")) {
//                            peptideInfo.isDecoy = false;
//                        }
//                    }
//                } else {
//                    PeptideInfo peptideInfo = new PeptideInfo(freePepSeq, protId.startsWith("DECOY_"), leftFlank, rightFlank);
//                    peptideInfo.protIdSet.add(protId);
//                    peptideInfoMap.put(freePepSeq, peptideInfo);
//                }
//
//            }
//        }
//        return;
//    }
    private int getNumOfMissCleavSite(ExpTag tag) {
        String str1 = tag.getFreeAaString();
        String str2 = str1.replaceAll("[KR]","");
        int numMC = str1.length()-str2.length();
        if (tag.isNorC == C_TAG && isKR(str1.charAt(str1.length()-1))) {
            numMC--; // the last KR at pepC does not count as MC
        }
        return numMC;
    }
    private boolean isKR(char aa){
        return aa == 'K' || aa == 'R';
    }

    private boolean isX(char aa){
        return aa == 'X';
    }
    public class PresearchEntry {

        public Map<String, PeptideInfo> peptideInfoMapForRef = new HashMap<>();
        public double precursorMass = PreSearch.this.precursorMass;
        public int precursorCharge = PreSearch.this.precursorCharge;
        public List<Peptide> ptmCandiList = new ArrayList<>();
        public List<Peptide> freeCandiList = new ArrayList<>();

        public String scanName;
        PresearchEntry() {
        }
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
        final String candidates;
        final String peptideSet;
        final int hasPTM;
        final int ptmNum;
        final int isSettled;
        final int shouldPtm;

        List<VarPtm> varPtmList = new ArrayList<>();
        Entry(int scanNum, String scanName, int shouldPtm, int hasPTM, int ptmNum, int isSetteld, int precursorCharge, double precursorMass
                ,String labelling, String peptide, double theoMass, int isDecoy
                , double score, double deltaLCn, double deltaCn, int matchedPeakNum, double ionFrac, double matchedHighestIntensityFrac
                , double explainedAaFrac, String otherPtmPatterns, String aScore, String candidates, String peptideSet) {
            this.scanNum = scanNum;
            this.scanName = scanName;
            this.shouldPtm = shouldPtm;
            this.hasPTM = hasPTM;
            this.ptmNum = ptmNum;
            this.isSettled = isSetteld;
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
            this.candidates = candidates;
            this.peptideSet = peptideSet;
        }
    }
}
