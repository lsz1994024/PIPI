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


import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;

public class ExpTag implements Comparable<ExpTag> {
    public int isNorC = 0; //N -1, none 0, C 1;
    private int hashCode;
    private double totalIntensity;
    private String freeAaString;
    private String ptmAaString;
    private int regionIdx;
    public List<ExpAa> expAaList;
    public List<Double> intensityList;

    public ExpTag(List<ExpAa> expAaList) {
        isNorC = expAaList.get(0).isNorC;
        intensityList = new ArrayList<>(expAaList.size()+1);
        this.expAaList = expAaList;
        intensityList.add( expAaList.get(0).getHeadIntensity() );
        for (int i = 0; i < expAaList.size(); i++) {
            intensityList.add( expAaList.get(i).getTailIntensity());
        }
        String toString = expAaList.get(0).toString();
        for (int i = 1; i < expAaList.size(); i++){
            toString += "-" + expAaList.get(i).toString();
        }
        hashCode = toString.hashCode();
        StringBuilder sbFreeSeq = new StringBuilder(expAaList.size());
        StringBuilder sbPtmSeq = new StringBuilder(2* expAaList.size() - 1);
        for (ExpAa aa : expAaList) {
            sbFreeSeq.append(aa.getPtmFreeAA());
            sbPtmSeq.append(aa.getAA());
        }
        freeAaString = sbFreeSeq.toString();
        ptmAaString = sbPtmSeq.toString();

        double intensity = expAaList.get(0).getHeadIntensity();
        for (ExpAa aa : expAaList) {
            intensity += aa.getTailIntensity();
        }
        totalIntensity = intensity;
    }

    public int hashCode() {
        return hashCode;
    }

    public boolean equals(Object other) {
        return (other instanceof ExpTag) && (this.hashCode() == other.hashCode());
    }

    public boolean approximateEquals(ExpTag other, double tolerance) {
        for (int i = 0; i < this.size(); ++i) {
            if (!this.expAaList.get(i).approximateEquals(other.expAaList.get(i), tolerance)) {
                return false;
            }
        }
        return true;
    }

    public int compareTo(ExpTag other) {
        return Double.compare(this.expAaList.get(0).getHeadLocation(), other.expAaList.get(0).getHeadLocation());
    }

//    public ExpAa[] getExpAAs() {
//        return expAaArray;
//    }

    public String getFreeAaString() {
        return freeAaString;
    }
    public String getPtmAaString() {
        return ptmAaString;
    }

    public String toString() {
        return isNorC+","+ ptmAaString ;
    }

    public double getTotalIntensity() {
        return totalIntensity;
    }

    public double getHeadLocation() {
        return expAaList.get(0).getHeadLocation();
    }

    public double getTailLocation() {
        return expAaList.get(expAaList.size()-1).getTailLocation();
    }

//    public ExpTag clone() throws CloneNotSupportedException {
//        super.clone();
//        if (expAaArray.length == 3 ){
//            return new ExpTag(expAaArray[0].clone(), expAaArray[1].clone(), expAaArray[2].clone());
//        }else {
//            return new ExpTag(expAaArray[0].clone(), expAaArray[1].clone(), expAaArray[2].clone(), expAaArray[3].clone());
//        }
//    }

    public int size() {
        return expAaList.size();
    }


    public ExpTag revTag(double totalMass){
        List<ExpAa> revExpAaList = new ArrayList<>(expAaList.size());
        for (int i = expAaList.size()-1; i>=0; i--) {
            revExpAaList.add(expAaList.get(i).revAA(totalMass));
        }
        ExpTag revedTag = new ExpTag(revExpAaList);
        revedTag.isNorC = isNorC;
        return revedTag;
    }

    public void appendTag(int pos, ExpTag newTag){
        int diffLen = this.size() - pos;
        if (newTag.size() <= diffLen) return;

        for (int i = diffLen; i < newTag.size(); i++) {
            this.intensityList.add(newTag.intensityList.get(i+1));
            this.expAaList.add(newTag.expAaList.get(i));
        }

        String toString = expAaList.get(0).toString();
        for (int i = 1; i < expAaList.size(); i++){
            toString += "-" + expAaList.get(i).toString();
        }
        hashCode = toString.hashCode();
        StringBuilder sbFreeSeq = new StringBuilder(expAaList.size());
        StringBuilder sbPtmSeq = new StringBuilder(2* expAaList.size() - 1);
        for (ExpAa aa : expAaList) {
            sbFreeSeq.append(aa.getPtmFreeAA());
            sbPtmSeq.append(aa.getAA());
        }
        freeAaString = sbFreeSeq.toString();
        ptmAaString = sbPtmSeq.toString();

        double intensity = expAaList.get(0).getHeadIntensity();
        for (ExpAa aa : expAaList) {
            intensity += aa.getTailIntensity();
        }
        totalIntensity = intensity;
    }

    public ExpTag subTag(int start, int end){ // start included  end not included
        List<ExpAa> subExpAaList = new ArrayList<>(end-start);
        for (int i = start; i < end; i++) {
            subExpAaList.add(expAaList.get(i));
        }
        ExpTag subTag = new ExpTag(subExpAaList);
//        if (isNorC == -1 )
        subTag.isNorC = isNorC;
        return subTag;
    }

}
