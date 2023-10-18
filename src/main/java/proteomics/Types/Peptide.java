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

package proteomics.Types;

import ProteomicsLibrary.Types.Coordinate;
import proteomics.Segment.InferSegment;
import ProteomicsLibrary.MassTool;
import java.util.*;

public class Peptide implements Comparable<Peptide>, Cloneable{
    public double nDeltaMass = -0.99;
    public double cDeltaMass = -0.99;
//    public boolean isTarget;
    public ExpTag finderTag = null;
    public int tagPosInPep = -1;

    public int scanNum = 0;
    private final String freeSeq;
    public String ptmSeq;
    public boolean isDecoy;
    private final String normalizedPeptideString;
    private final MassTool massTool;
    public double absDeltaMass = 0d;
    public Peptide bestPep = null;
    private int hashCode;

    // these fields need to be changed every time PTM changed.
    private PosMassMap varPtmMap = null;
    public TreeMap<Integer, VarPtm> posVarPtmResMap = new TreeMap<>();
    private double theoMass = -1;
    private double[][] ionMatrix = null;
    private String varPtmContainingSeq = null;

    private String ptmContainingSeq = null;

    // score part
    private double score = -1;
    private int matchedPeakNum = -1;
    private double ionFrac = -1;
    private double matchedHighestIntensityFrac = -1;
    private double explainedAaFrac = -1;
    private double qValue = -1;
    private String aScore = "-";
    public Map<Integer, Double> matchedBions = new HashMap<>();
    public Map<Integer, Double> matchedYions = new HashMap<>();
    public double precursorMass;
    public double calculatedMass;

    public Peptide(String freeSeq, boolean isDecoy, MassTool massTool) {
        this.freeSeq = freeSeq;
        this.ptmSeq = freeSeq;
        this.isDecoy = isDecoy;
        this.normalizedPeptideString = InferSegment.normalizeSequence(freeSeq);
        this.massTool = massTool;

        hashCode = freeSeq.hashCode();
    }

    public void setScanNum(int scanNum) {this.scanNum = scanNum;}


    public int getPriority(){
        int res = 0;
        for (VarPtm varPtm : posVarPtmResMap.values()){
            res += varPtm.priority;
        }
        return res;
    }
    public double[][] getIonMatrix() {
        if (ionMatrix == null) {
            varPtmContainingSeq = getVarPtmContainingSeq();
            ionMatrix = massTool.buildIonArray(varPtmContainingSeq);
            theoMass = massTool.calResidueMass(varPtmContainingSeq) + massTool.H2O;
        }
        return ionMatrix;
    }

    public double[][] getIonMatrixNow() {
        varPtmContainingSeq = getVarPtmContainingSeqNow();
        ionMatrix = massTool.buildIonArray(varPtmContainingSeq);
        theoMass = massTool.calResidueMass(varPtmContainingSeq) + massTool.H2O;
        return ionMatrix;
    }

    public String getVarPtmContainingSeqNow() {
        if (varPtmMap != null) {
            StringBuilder sb = new StringBuilder(freeSeq.length() * 5);
            int tempIdx = varPtmMap.firstKey()+1;
            if (tempIdx > 1) {
                sb.append(freeSeq, 0, tempIdx - 1);
            }
            int i = tempIdx - 1;
            tempIdx = varPtmMap.lastKey()+1;
            while (i < freeSeq.length()) {
                boolean hasMod = false;
                if (tempIdx > i) {
                    for (Integer co : varPtmMap.keySet()) {
                        if (co == i) {
                            sb.append(String.format(Locale.US, "%c(%.3f)", freeSeq.charAt(i), varPtmMap.get(co)));
                            hasMod = true;
                            ++i;
                            break;
                        }
                    }
                    if (!hasMod) {
                        sb.append(freeSeq.charAt(i));
                        ++i;
                    }
                } else {
                    break;
                }
            }
            if (tempIdx < freeSeq.length()) {
                sb.append(freeSeq.substring(tempIdx));
            }
            varPtmContainingSeq = sb.toString();
        } else {
            varPtmContainingSeq = freeSeq;
        }
        return varPtmContainingSeq;
    }

    public String getNormalizedPeptideString() {
        return normalizedPeptideString;
    }

    public boolean isDecoy() {
        return isDecoy;
    }

    public double getTheoMass() {
        if (theoMass < 0) {
            varPtmContainingSeq = getVarPtmContainingSeq();
            ionMatrix = massTool.buildIonArray(varPtmContainingSeq);
            theoMass = massTool.calResidueMass(varPtmContainingSeq) + massTool.H2O;
//            chargeOneBIonArray = ionMatrix[0];
        }
        return theoMass;
    }

    @Override
    public boolean equals(Object other) {
        if (!(other instanceof Peptide)) {
            return false;
        }

        Peptide otherPeptide = (Peptide) other;
        return this.hashCode == otherPeptide.hashCode;
    }

    public Peptide clone() throws CloneNotSupportedException {
        super.clone();//???
        Peptide other = new Peptide(freeSeq, isDecoy, massTool);
        if (varPtmMap != null) {
            other.setVarPTM(varPtmMap.clone());
            other.setScore(score);
            other.setMatchedHighestIntensityFrac(matchedHighestIntensityFrac);
            other.setExplainedAaFrac(explainedAaFrac);
            other.setIonFrac(ionFrac);
            other.setaScore(aScore);
            other.setQValue(qValue);
            if (this.posVarPtmResMap != null){
                other.posVarPtmResMap = new TreeMap<>(this.posVarPtmResMap);
            }
        }
        other.nDeltaMass = this.nDeltaMass;
        other.cDeltaMass = this.cDeltaMass;
//        other.isTarget = this.isTarget;
        other.finderTag = null;
        if (this.finderTag != null){
            other.finderTag = new ExpTag(this.finderTag.expAaList);
        }
        other.tagPosInPep = this.tagPosInPep;
        return other;
    }

    public int hashCode() {
        return hashCode;
    }

    public int length() {
        return freeSeq.length();
    }

    public void setVarPTM(PosMassMap ptmMap) {
        this.varPtmMap = ptmMap;
        if (ptmMap != null) {
            // reset these fields to make them being regenerated again.
            theoMass = -1;
            ionMatrix = null;
//            chargeOneBIonArray = null;
            varPtmContainingSeq = null;
            ptmContainingSeq = null;

            String toString = freeSeq + "." + ptmMap.toString();
            hashCode = toString.hashCode();
        }
    }

    private int getTotalPriority(){//todo
        int res = 0;
//        for (int )

        return res;
    }

    public String toString() {
        if (varPtmMap == null) {
            return freeSeq;
        }
        return freeSeq + "." + varPtmMap.toString();
    }

    public boolean hasVarPTM() {
        return varPtmMap != null;
    }

    private int getVarPTMNum() {
        if (hasVarPTM()) {
            return varPtmMap.size();
        } else {
            return 0;
        }
    }

    public String getFreeSeq() {
        return freeSeq;
    }

    public PosMassMap getVarPTMs() {
        return varPtmMap;
    }

    private String getVarPtmContainingSeq() {
        if (varPtmContainingSeq == null) {
            if (varPtmMap != null) {
                StringBuilder sb = new StringBuilder(freeSeq.length() * 5);

                int tempIdx = 0;
                try {
                    tempIdx = varPtmMap.firstKey()+1;

                } catch (Exception ex){
                    System.out.println("error");
                }
                if (tempIdx > 1) {
                    sb.append(freeSeq.substring(0, tempIdx - 1));
                }
                int i = tempIdx - 1;
                tempIdx = varPtmMap.lastKey()+1;
                while (i < freeSeq.length()) {
                    boolean hasMod = false;
                    if (tempIdx > i) {
                        for (Integer co : varPtmMap.keySet()) {
                            if (co == i) {
                                sb.append(String.format(Locale.US, "%c(%.3f)", freeSeq.charAt(i), varPtmMap.get(co)));
                                hasMod = true;
                                ++i;
                                break;
                            }
                        }
                        if (!hasMod) {
                            sb.append(freeSeq.charAt(i));
                            ++i;
                        }
                    } else {
                        break;
                    }
                }
                if (tempIdx < freeSeq.length()) {
                    sb.append(freeSeq.substring(tempIdx));
                }
                varPtmContainingSeq = sb.toString();
            } else {
                varPtmContainingSeq = freeSeq;
            }
        }

        return varPtmContainingSeq;
    }

    public String getPtmContainingSeq(Map<Character, Double> fixModMap) { // caution: containing fix modification. Calculating ion masses based on it is incorrect.
        if (ptmContainingSeq == null) {
            ptmContainingSeq = getVarPtmContainingSeq();
            for (char aa : fixModMap.keySet()) {
                if (Math.abs(fixModMap.get(aa)) > 0.01) {
                    ptmContainingSeq = ptmContainingSeq.replaceAll(String.valueOf(aa), String.format(Locale.US, "%c(%.3f)", aa, fixModMap.get(aa)));
                }
            }
        }

        return ptmContainingSeq;
    }

//    public String getPtmContainingSeq(Map<Character, Double> fixModMap) {
//        if (posVarPtmResMap != null) {
//            StringBuilder sb = new StringBuilder(freeSeq);
//            for (int pos : posVarPtmResMap.descendingKeySet()) {
//                if (pos == 0 && (posVarPtmResMap.get(pos).position == 0 || posVarPtmResMap.get(pos).position == 2)) {
//                    sb.replace(pos, pos+1, String.format(Locale.US, "(n%.3f)%c", posVarPtmResMap.get(pos).mass, freeSeq.charAt(pos)));
//                } else if (pos == freeSeq.length()-1 && (posVarPtmResMap.get(pos).position == 1 || posVarPtmResMap.get(pos).position == 3)) {
//                    sb.replace(pos, pos+1, String.format(Locale.US, "%c(c%.3f)", freeSeq.charAt(pos), posVarPtmResMap.get(pos).mass));
//                } else {
//                    sb.replace(pos, pos+1, String.format(Locale.US, "%c(%.3f)", freeSeq.charAt(pos), posVarPtmResMap.get(pos).mass));
//                }
//            }
//            ptmContainingSeq = sb.toString();
//        } else {
//            ptmContainingSeq = freeSeq;
//        }
//
//        // this below for-loop should be commented when developping to hide C(57.021) for better looking
////        for (char aa : fixModMap.keySet()) {
////            if (Math.abs(fixModMap.get(aa)) > 0.01) {
////                ptmContainingSeq = ptmContainingSeq.replaceAll(String.valueOf(aa), String.format(Locale.US, "%c(%.3f)", aa, fixModMap.get(aa)));
////            }
////        }
//        return ptmContainingSeq;
//    }

    public void setScore(double score) {
        this.score = score;
    }

    public void setMatchedPeakNum(int matchedPeakNum) {
        this.matchedPeakNum = matchedPeakNum;
    }

    public void setIonFrac(double ionFrac) {
        this.ionFrac = ionFrac;
    }

    public void setMatchedHighestIntensityFrac(double matchedHighestIntensityFrac) {
        this.matchedHighestIntensityFrac = matchedHighestIntensityFrac;
    }

    public void setExplainedAaFrac(double explainedAaFrac) {
        this.explainedAaFrac = explainedAaFrac;
    }

    public void setaScore(String aScore) {
        this.aScore = aScore;
    }

    private void setQValue(double qValue) {
        this.qValue = qValue;
    }

    public double getScore() {
        return score;
    }

    public int getMatchedPeakNum() {
        return matchedPeakNum;
    }

    public double getIonFrac() {
        return ionFrac;
    }

    public double getMatchedHighestIntensityFrac() {
        return matchedHighestIntensityFrac;
    }

    public double getExplainedAaFrac() {
        return explainedAaFrac;
    }

    public String getaScore() {
        return aScore;
    }

    public double getQValue() {
        return qValue;
    }

//    public int compareTo(Peptide peptide) {
//        if (score > peptide.getScore()) {
//            return 1;
//        } else if (score < peptide.getScore()) {
//            return -1;
//        } else {
//            if (matchedPeakNum > peptide.getMatchedPeakNum()) {
//                return 1;
//            } else if (matchedPeakNum < peptide.getMatchedPeakNum()) {
//                return -1;
//            } else {
//                if (explainedAaFrac > peptide.getExplainedAaFrac()) {
//                    return 1;
//                } else if (explainedAaFrac < peptide.getExplainedAaFrac()) {
//                    return -1;
//                } else {
//                    if (getVarPTMNum() < peptide.getVarPTMNum()) {
//                        return 1;
//                    } else if (getVarPTMNum() > peptide.getVarPTMNum()) {
//                        return -1;
//                    } else if (tagVecScore > peptide.getTagVecScore()) {
//                        return 1;
//                    } else if (tagVecScore < peptide.getTagVecScore()) {
//                        return -1;
//                    } else {
//                        if (!isDecoy && peptide.isDecoy()) {
//                            return 1;
//                        } else if (isDecoy && !peptide.isDecoy()) {
//                            return -1;
//                        } else{
//                            if (absDeltaMass < peptide.absDeltaMass){
//                                return 1;
//                            } else if (absDeltaMass > peptide.absDeltaMass){
//                                return -1;
//                            } else {
//                                return 0;
//                            }
//                        }
//                    }
//                }
//            }
//        }
//    }

    public int compareTo(Peptide peptide) {
        if (hasVarPTM() && peptide.hasVarPTM()) { // only when two peptide are both modified, use priority to compare them in the very beginning
            if (score*this.getPriority() > peptide.score*peptide.getPriority()){
                return 1;
            } else if(score*this.getPriority() < peptide.score*peptide.getPriority()){
                return -1;
            } else {
                return this.subCompareTo(peptide);
            }
        } else {
            return this.subCompareTo(peptide);
        }
    }

    public int subCompareTo(Peptide peptide) {
        if (score > peptide.getScore()) {
            return 1;
        } else if (score < peptide.getScore()) {
            return -1;
        } else {
            if (matchedPeakNum > peptide.getMatchedPeakNum()) {
                return 1;
            } else if (matchedPeakNum < peptide.getMatchedPeakNum()) {
                return -1;
            } else {
                if (getPriority() > peptide.getPriority()) {
                    return 1;
                } else if (getPriority() < peptide.getPriority()) {
                    return -1;
                } else {
                    if (getVarPTMNum() < peptide.getVarPTMNum()) {
                        return 1;
                    } else if (getVarPTMNum() > peptide.getVarPTMNum()) {
                        return -1;
                    } else {
                        if (this.hashCode < peptide.hashCode){
                            return -1;
                        } else if (this.hashCode > peptide.hashCode) {
                            return 1;
                        } else {
                            return 0;
                        }
                    }
                }
            }
        }
    }
//    public int compareTo(Peptide peptide) {
//        if (score > peptide.getScore()) {
//            return 1;
//        } else if (score < peptide.getScore()) {
//            return -1;
//        } else {
//            if (matchedPeakNum > peptide.getMatchedPeakNum()) {
//                return 1;
//            } else if (matchedPeakNum < peptide.getMatchedPeakNum()) {
//                return -1;
//            } else {
//                if (getPriority() > peptide.getPriority()) {
//                    return 1;
//                } else if (getPriority() < peptide.getPriority()) {
//                    return -1;
//                } else {
//                    if (getVarPTMNum() < peptide.getVarPTMNum()) {
//                        return 1;
//                    } else if (getVarPTMNum() > peptide.getVarPTMNum()) {
//                        return -1;
//                    } else {
//                        if (finderTag != null && peptide.finderTag != null) {
//                            if (finderTag.getTotalIntensity() > peptide.finderTag.getTotalIntensity()) {
//                                return 1;
//                            } else if (finderTag.getTotalIntensity() < peptide.finderTag.getTotalIntensity()) {
//                                return -1;
//                            } else {
//                                if (finderTag.size() < peptide.finderTag.size()){
//                                    return 1;
//                                } else if (finderTag.size() > peptide.finderTag.size()){
//                                    return -1;
//                                } else {
//                                    return 0;
//                                }
//                            }
//                        } else {
//                            return 0;
//                        }
//                    }
//                }
//            }
//        }
//    }
}
