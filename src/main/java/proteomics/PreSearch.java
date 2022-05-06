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
import ProteomicsLibrary.PrepareSpectrum;
import ProteomicsLibrary.Types.SparseVector;
import proteomics.Index.BuildIndex;
import proteomics.PTM.InferPTM;
import proteomics.Search.Search;
import proteomics.Segment.InferSegment;
import proteomics.Spectrum.PreSpectra;
import proteomics.Types.*;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;

import java.util.*;
import java.util.concurrent.Callable;
import java.util.concurrent.locks.ReentrantLock;

public class PreSearch implements Callable<PreSearch.Entry> {
    private static final int candisNum = 20;
    private final BuildIndex buildIndex;
    private final MassTool massTool;
    private final double ms1Tolerance;
    private final double leftInverseMs1Tolerance;
    private final double rightInverseMs1Tolerance;
    private final int ms1ToleranceUnit;
    private final double ms2Tolerance;
    private final double minPtmMass;
    private final double maxPtmMass;
    private final int localMaxMs2Charge;
    private final Map<String, Peptide0> peptide0Map;
    private final JMzReader spectraParser;
    private final double minClear;
    private final double maxClear;
    private final ReentrantLock lock;
    private final String scanId;
    private final int precursorCharge;
    private final double precursorMass;
    private final InferPTM inferPTM;
    private final PrepareSpectrum preSpectrum;
    private final String sqlPath;
    private final int scanNum;
    private final int precursorScanNo;



    public PreSearch(int scanNum, BuildIndex buildIndex, MassTool massTool, double ms1Tolerance, double leftInverseMs1Tolerance, double rightInverseMs1Tolerance
            , int ms1ToleranceUnit, double ms2Tolerance, double minPtmMass, double maxPtmMass, int localMaxMs2Charge
            , JMzReader spectraParser, double minClear, double maxClear, ReentrantLock lock, String scanId, int precursorCharge, double precursorMass
            , InferPTM inferPTM, PrepareSpectrum preSpectrum, String sqlPath, int precursorScanNo) {

        this.buildIndex = buildIndex;
        this.massTool = massTool;
        this.ms1Tolerance = ms1Tolerance;
        this.leftInverseMs1Tolerance = leftInverseMs1Tolerance;
        this.rightInverseMs1Tolerance = rightInverseMs1Tolerance;
        this.ms1ToleranceUnit = ms1ToleranceUnit;
        this.ms2Tolerance = ms2Tolerance;
        this.minPtmMass = minPtmMass;
        this.maxPtmMass = maxPtmMass;
        this.localMaxMs2Charge = localMaxMs2Charge;
        this.spectraParser = spectraParser;
        this.minClear = minClear;
        this.maxClear = maxClear;
        this.lock = lock;
        this.scanId = scanId;
        this.precursorCharge = precursorCharge;
        this.precursorMass = precursorMass;
        this.inferPTM = inferPTM;
        this.preSpectrum = preSpectrum;
        this.sqlPath = sqlPath;
        peptide0Map = buildIndex.getPeptide0Map();
        this.scanNum = scanNum;
        this.precursorScanNo = precursorScanNo;
    }

    @Override
    public Entry call() throws Exception {
        Map<Double, Double> rawPLMap;
        try {
            lock.lock();
            // Reading peak list.
            rawPLMap = spectraParser.getSpectrumById(scanId).getPeakList();
        } finally {
            lock.unlock();
        }
        // preprocess peak list
        TreeMap<Double, Double> plMap = preSpectrum.preSpectrumTopNStyle(rawPLMap, precursorMass, precursorCharge, minClear, maxClear, PreSpectra.topN);

        if (plMap.isEmpty()) {
            return null;
        }

        // Coding
        InferSegment inferSegment = buildIndex.getInferSegment();
        List<ThreeExpAA> expAaLists = inferSegment.inferSegmentLocationFromSpectrum(precursorMass, plMap, scanNum);

        if (!expAaLists.isEmpty()) {
            SparseVector scanCode = inferSegment.generateSegmentIntensityVector(expAaLists);

            Entry entry = new Entry();
            Search search = new Search(entry, scanNum, buildIndex, precursorMass, scanCode, massTool, ms1Tolerance, leftInverseMs1Tolerance, rightInverseMs1Tolerance, ms1ToleranceUnit, minPtmMass, maxPtmMass, localMaxMs2Charge);
//            entry.candidateList = search.candidatesList;
            //move search to here as presearch.

            return entry;
        } else {
            return null;
        }
    }


    public class Entry {

        public int scanNum = PreSearch.this.scanNum;
        public double precursorMass = PreSearch.this.precursorMass;
        public List<Peptide> ptmOnlyList = new ArrayList<>();
        public List<Peptide> ptmFreeList = new ArrayList<>();
        public String candidateList = null;
        Entry() {
//            this.scanNum = scanNum;
        }
//        Entry(List<Peptide> ptmOnlyList, List<Peptide> ptmFreeList) {
//
//            this.ptmOnlyList = ptmOnlyList;
//            this.ptmFreeList = ptmFreeList;
//        }
    }
}
