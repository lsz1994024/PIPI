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
    public double nDeltaMass = 0.0;
    public double cDeltaMass = 0.0;
    public boolean isTarget;
    public ExpTag tagFinder = null;
    public int posInPepSeq = -1;

    public int scanNum = 0;
    private final String freeSeq;
    public String ptmSeq;
    public boolean isDecoy;
    private final String normalizedPeptideString;
    private final MassTool massTool;
    private final int maxMs2Charge;
    private final int globalRank;
    private final double tagVecScore;
    public double absDeltaMass = 0d;
    public Peptide bestPep = null;
    private int hashCode;
    public boolean shouldPTM = false;

    // these fields need to be changed every time PTM changed.
    private PositionDeltaMassMap varPTMMap = null;
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

    public Peptide(String freeSeq, boolean isDecoy, MassTool massTool, int maxMs2Charge, double tagVecScore, int globalRank) {
        this.freeSeq = freeSeq;
        this.ptmSeq = freeSeq;
        this.isDecoy = isDecoy;
        this.normalizedPeptideString = InferSegment.normalizeSequence(freeSeq);
        this.tagVecScore = tagVecScore;
        this.massTool = massTool;
        this.maxMs2Charge = maxMs2Charge;
        this.globalRank = globalRank;

        hashCode = freeSeq.hashCode();
    }

    public void setScanNum(int scanNum) {this.scanNum = scanNum;}

    public int getGlobalRank() {
        return globalRank;
    }

    public double[][] getIonMatrix() {
        if (ionMatrix == null) {
            varPtmContainingSeq = getVarPtmContainingSeq();
            ionMatrix = massTool.buildIonArray(varPtmContainingSeq, maxMs2Charge);
            theoMass = massTool.calResidueMass(varPtmContainingSeq) + massTool.H2O;
        }
        return ionMatrix;
    }

    public double[][] getIonMatrixNow() {
        varPtmContainingSeq = getVarPtmContainingSeqNow();
        ionMatrix = massTool.buildIonArray(varPtmContainingSeq, maxMs2Charge);
        theoMass = massTool.calResidueMass(varPtmContainingSeq) + massTool.H2O;
        return ionMatrix;
    }

    public String getVarPtmContainingSeqNow() {
        if (varPTMMap != null) {
            StringBuilder sb = new StringBuilder(freeSeq.length() * 5);
            int tempIdx = varPTMMap.firstKey().y;
            if (tempIdx > 1) {
                sb.append(freeSeq.substring(0, tempIdx - 1));
            }
            int i = tempIdx - 1;
            tempIdx = varPTMMap.lastKey().y;
            while (i < freeSeq.length()) {
                boolean hasMod = false;
                if (tempIdx > i) {
                    for (Coordinate co : varPTMMap.keySet()) {
                        if (co.y - 1 == i) {
                            sb.append(String.format(Locale.US, "%c(%.3f)", freeSeq.charAt(i), varPTMMap.get(co)));
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

//    public double[][] getTheoIonMatrix() {
//        if (theoIonMatrix == null) {
//            varPtmContainingSeq = getVarPtmContainingSeq();
//            theoIonMatrix = massTool.buildTheoBYIonsArray(varPtmContainingSeq, 1);
//            theoMass = massTool.calResidueMass(varPtmContainingSeq) + massTool.H2O;
//        }
//        return theoIonMatrix;
//    }

    public String getNormalizedPeptideString() {
        return normalizedPeptideString;
    }

    public boolean isDecoy() {
        return isDecoy;
    }

    public double getTheoMass() {
        if (theoMass < 0) {
            varPtmContainingSeq = getVarPtmContainingSeq();
            ionMatrix = massTool.buildIonArray(varPtmContainingSeq, maxMs2Charge);
            theoMass = massTool.calResidueMass(varPtmContainingSeq) + massTool.H2O;
//            chargeOneBIonArray = ionMatrix[0];
        }
        return theoMass;
    }

//    public double[] getChargeOneBIonArray() {
//        if (chargeOneBIonArray == null) {
//            varPtmContainingSeq = getVarPtmContainingSeq();
//            ionMatrix = massTool.buildIonArray(varPtmContainingSeq, maxMs2Charge);
//            theoMass = massTool.calResidueMass(varPtmContainingSeq) + massTool.H2O;
//            chargeOneBIonArray = ionMatrix[0];
//        }
//        return chargeOneBIonArray;
//    }

    public boolean equals(Object other) {
        if (!(other instanceof Peptide)) {
            return false;
        }

        Peptide otherPeptide = (Peptide) other;
        return this.hashCode == otherPeptide.hashCode;
    }

    public Peptide clone() throws CloneNotSupportedException {
        super.clone();//???
        Peptide other = new Peptide(freeSeq, isDecoy, massTool, maxMs2Charge, tagVecScore, globalRank);
        if (varPTMMap != null) {
            other.setVarPTM(varPTMMap.clone());
            other.setScore(score);
            other.setMatchedHighestIntensityFrac(matchedHighestIntensityFrac);
            other.setExplainedAaFrac(explainedAaFrac);
            other.setIonFrac(ionFrac);
            other.setaScore(aScore);
            other.setQValue(qValue);
        }

        return other;
    }

    public int hashCode() {
        return hashCode;
    }

    public int length() {
        return freeSeq.length();
    }

    public void setVarPTM(PositionDeltaMassMap ptmMap) {
        this.varPTMMap = ptmMap;
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

    public String toString() {
        if (varPTMMap == null) {
            return freeSeq;
        }
        return freeSeq + "." + varPTMMap.toString();
    }

    public boolean hasVarPTM() {
        return varPTMMap != null;
    }

    private int getVarPTMNum() {
        if (hasVarPTM()) {
            return varPTMMap.size();
        } else {
            return 0;
        }
    }

    public String getPTMFreePeptide() {
        return freeSeq;
    }

    public PositionDeltaMassMap getVarPTMs() {
        return varPTMMap;
    }

    private String getVarPtmContainingSeq() {
        if (varPtmContainingSeq == null) {
            if (varPTMMap != null) {
                StringBuilder sb = new StringBuilder(freeSeq.length() * 5);

                int tempIdx = 0;
                try {
                    tempIdx = varPTMMap.firstKey().y;

                } catch (Exception ex){
                    System.out.println("error");
                }
                if (tempIdx > 1) {
                    sb.append(freeSeq.substring(0, tempIdx - 1));
                }
                int i = tempIdx - 1;
                tempIdx = varPTMMap.lastKey().y;
                while (i < freeSeq.length()) {
                    boolean hasMod = false;
                    if (tempIdx > i) {
                        for (Coordinate co : varPTMMap.keySet()) {
                            if (co.y - 1 == i) {
                                sb.append(String.format(Locale.US, "%c(%.3f)", freeSeq.charAt(i), varPTMMap.get(co)));
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

    public double getNormalizedCrossCorr() {
        return tagVecScore;
    }

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

    public int compareTo(Peptide peptide) {
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
                if (explainedAaFrac > peptide.getExplainedAaFrac()) {
                    return 1;
                } else if (explainedAaFrac < peptide.getExplainedAaFrac()) {
                    return -1;
                } else {
                    if (getVarPTMNum() < peptide.getVarPTMNum()) {
                        return 1;
                    } else if (getVarPTMNum() > peptide.getVarPTMNum()) {
                        return -1;
                    } else if (tagVecScore > peptide.getNormalizedCrossCorr()) {
                        return 1;
                    } else if (tagVecScore < peptide.getNormalizedCrossCorr()) {
                        return -1;
                    } else {
                        if (!isDecoy && peptide.isDecoy()) {
                            return 1;
                        } else if (isDecoy && !peptide.isDecoy()) {
                            return -1;
                        } else{
                            if (absDeltaMass < peptide.absDeltaMass){
                                return 1;
                            } else if (absDeltaMass > peptide.absDeltaMass){
                                return -1;
                            } else {
                                return 0;
                            }
                        }
                    }
                }
            }
        }
    }
}
