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


import org.apache.commons.math3.analysis.function.Exp;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ThreeExpAA implements Comparable<ThreeExpAA> {

    private final ExpAA[] threeExpAa;
    private int hashCode;
    private final double totalIntensity;
    private final String ptmFreeAAString;
    private int regionIdx;
    public Map<Integer, Double> bAlignPosMassMap = new HashMap<>();
    public Map<Integer, Double> yAlignPosMassMap = new HashMap<>();
    public List<ExpAA> expAAList = null;
    public double[] intensityArray;
    public ThreeExpAA(ExpAA aa1, ExpAA aa2, ExpAA aa3) {
        threeExpAa = new ExpAA[]{aa1, aa2, aa3};
        String toString = threeExpAa[0].toString() + "-" + threeExpAa[1].toString() + "-" + threeExpAa[2].toString();
        hashCode = toString.hashCode();

        StringBuilder sb = new StringBuilder(5);
        for (ExpAA aa : threeExpAa) {
            sb.append(aa.getPtmFreeAA());
        }
        ptmFreeAAString = sb.toString();

        double intensity = threeExpAa[0].getHeadIntensity();
        for (ExpAA aa : threeExpAa) {
            intensity += aa.getTailIntensity();
        }
        totalIntensity = intensity;
    }
    public ThreeExpAA(ExpAA aa1, ExpAA aa2, ExpAA aa3, ExpAA aa4) {
        threeExpAa = new ExpAA[]{aa1, aa2, aa3, aa4};
        String toString = threeExpAa[0].toString() + "-" + threeExpAa[1].toString() + "-" + threeExpAa[2].toString()+ "-" + threeExpAa[3].toString();
        hashCode = toString.hashCode();

        StringBuilder sb = new StringBuilder(7);
        for (ExpAA aa : threeExpAa) {
            sb.append(aa.getPtmFreeAA());
        }
        ptmFreeAAString = sb.toString();

        double intensity = threeExpAa[0].getHeadIntensity();
        for (ExpAA aa : threeExpAa) {
            intensity += aa.getTailIntensity();
        }
        totalIntensity = intensity;
    }

    public ThreeExpAA(ExpAA aa1, ExpAA aa2, ExpAA aa3, ExpAA aa4, ExpAA aa5) {
        threeExpAa = new ExpAA[]{aa1, aa2, aa3, aa4, aa5};
        String toString = threeExpAa[0].toString() + "-" + threeExpAa[1].toString() + "-" + threeExpAa[2].toString()+ "-" + threeExpAa[3].toString()+ "-" + threeExpAa[4].toString();
        hashCode = toString.hashCode();

        StringBuilder sb = new StringBuilder(9);
        for (ExpAA aa : threeExpAa) {
            sb.append(aa.getPtmFreeAA());
        }
        ptmFreeAAString = sb.toString();

        double intensity = threeExpAa[0].getHeadIntensity();
        for (ExpAA aa : threeExpAa) {
            intensity += aa.getTailIntensity();
        }
        totalIntensity = intensity;
    }
    public ThreeExpAA(List<ExpAA> expAAList) {
        intensityArray = new double[expAAList.size()+1];
        this.expAAList = expAAList;
        intensityArray[0] = expAAList.get(0).getHeadIntensity();
        for (int i = 0; i < expAAList.size(); i++) {
            intensityArray[i+1] = expAAList.get(i).getTailIntensity();
        }
        threeExpAa = expAAList.toArray(new ExpAA[0]);
        String toString = threeExpAa[0].toString();
        for (int i = 1; i < threeExpAa.length; i++){
            toString += "-" + threeExpAa[i].toString();
        }
        hashCode = toString.hashCode();
        StringBuilder sb = new StringBuilder(2*threeExpAa.length - 1);
        for (ExpAA aa : threeExpAa) {
            sb.append(aa.getPtmFreeAA());
        }
        ptmFreeAAString = sb.toString();

        double intensity = threeExpAa[0].getHeadIntensity();
        for (ExpAA aa : threeExpAa) {
            intensity += aa.getTailIntensity();
        }
        totalIntensity = intensity;
    }

    public int hashCode() {
        return hashCode;
    }

    public boolean equals(Object other) {
        return (other instanceof ThreeExpAA) && (this.hashCode() == other.hashCode());
    }

    public boolean approximateEquals(ThreeExpAA other, double tolerance) {
        for (int i = 0; i < this.size(); ++i) {
            if (!this.get(i).approximateEquals(other.get(i), tolerance)) {
                return false;
            }
        }
        return true;
    }

    public int compareTo(ThreeExpAA other) {
        return Double.compare(threeExpAa[0].getHeadLocation(), other.threeExpAa[0].getHeadLocation());
    }

    public ExpAA[] getExpAAs() {
        return threeExpAa;
    }

    public String getPtmFreeAAString() {
        return ptmFreeAAString;
    }

    public String toString() {
        return ptmFreeAAString;
    }

    public double getTotalIntensity() {
        return totalIntensity;
    }

    public double getHeadLocation() {
        return threeExpAa[0].getHeadLocation();
    }

    public double getTailLocation() {
        return threeExpAa[threeExpAa.length - 1].getTailLocation();
    }

    public ThreeExpAA clone() throws CloneNotSupportedException {
        super.clone();
        if (threeExpAa.length == 3 ){
            return new ThreeExpAA(threeExpAa[0].clone(), threeExpAa[1].clone(), threeExpAa[2].clone());
        }else {
            return new ThreeExpAA(threeExpAa[0].clone(), threeExpAa[1].clone(), threeExpAa[2].clone(),threeExpAa[3].clone());
        }
    }

    public int size() {
        return threeExpAa.length;
    }

    public ExpAA get(int i) {
        return threeExpAa[i];
    }

    public void setRegionIdx(int regionIdx) {
        this.regionIdx = regionIdx;
    }

    public int getRegionIdx() {
        return regionIdx;
    }

    public ThreeExpAA revTag(double totalMass){
        List<ExpAA> revExpAAList = new ArrayList<>(threeExpAa.length);
        for (int i = threeExpAa.length-1; i>=0; i--) {
            revExpAAList.add(threeExpAa[i].revAA(totalMass));
        }
        return new ThreeExpAA(revExpAAList);
    }

    public double getIntesFrom(int aaId1, int aaId2){
        double intes = 0;
        for (int i = aaId1; i <= aaId2+1; i++){
            intes += intensityArray[i];
        }
        return intes;
    }
}
