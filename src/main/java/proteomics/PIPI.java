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

import ProteomicsLibrary.Types.SparseVector;
import org.apache.commons.math3.util.Pair;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.PTM.InferPTM;
import ProteomicsLibrary.Binomial;
import ProteomicsLibrary.SpecProcessor;
import proteomics.Types.*;
import proteomics.Index.BuildIndex;
import proteomics.Parameter.Parameter;
import proteomics.Spectrum.DatasetReader;
import ProteomicsLibrary.MassTool;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.mgf_parser.MgfFile;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLFile;

import java.io.*;
import java.net.InetAddress;
import java.net.UnknownHostException;
import java.sql.*;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.locks.ReentrantLock;

public class PIPI {

    private static final Logger logger = LoggerFactory.getLogger(PIPI.class);
    public static final String versionStr = "1.4.7";
    static final boolean useXcorr = false;

    public static final int[] debugScanNumArray = new int[]{};

    public static void main(String[] args) {
        long startTime = System.nanoTime();

        // Process inputs
        if (args.length != 2) {
            help();
        }

        // Set parameters
        String parameterPath = args[0].trim();
        String spectraPath = args[1].trim();

        logger.info("Running PIPI version {}.", versionStr);

        String dbName = null;
        String hostName = "unknown-host";
        try {
            hostName = InetAddress.getLocalHost().getHostName();
            logger.info("Computer: {}.", hostName);
        } catch (UnknownHostException ex) {
            logger.warn("Cannot get the computer's name.");
        }

        try {
            logger.info("Spectra: {}, parameter: {}.", spectraPath, parameterPath);

            dbName = String.format(Locale.US, "PIPI.%s.%s.temp.db", hostName, new SimpleDateFormat("yyyy-MM-dd-HH-mm-ss").format(Calendar.getInstance().getTime()));
            new PIPI(parameterPath, spectraPath, dbName);



        } catch (Exception ex) {
            ex.printStackTrace();
            logger.error(ex.toString());
        } finally {
            if (dbName != null) {
                (new File(dbName)).delete();
                (new File(dbName + "-wal")).delete();
                (new File(dbName + "-shm")).delete();
            }
        }

        double totalHour = (double) (System.nanoTime() - startTime) * 1e-9 / 3600;
        logger.info("Running time: {} hours.", totalHour);
        logger.info("Done!");
    }

    private PIPI(String parameterPath, String spectraPath, String dbName) throws Exception {
        // Get the parameter map
        Parameter parameter = new Parameter(parameterPath);
        Map<String, String> parameterMap = parameter.returnParameterMap();
        double ms2Tolerance = Double.valueOf(parameterMap.get("ms2_tolerance"));
        double ms1Tolerance = Double.valueOf(parameterMap.get("ms1_tolerance"));
        double leftInverseMs1Tolerance = 1 / (1 + ms1Tolerance * 1e-6);
        double rightInverseMs1Tolerance = 1 / (1 - ms1Tolerance * 1e-6);
        int ms1ToleranceUnit = Integer.valueOf(parameterMap.get("ms1_tolerance_unit"));
        double minClear = Double.valueOf(parameterMap.get("min_clear_mz"));
        double maxClear = Double.valueOf(parameterMap.get("max_clear_mz"));
        String percolatorPath = parameterMap.get("percolator_path");
        boolean outputPercolatorInput = (Integer.valueOf(parameterMap.get("output_percolator_input")) == 1);

        // Check if Percolator can be executed.
        if (!(new File(percolatorPath)).exists()) {
            throw new NullPointerException(String.format(Locale.US, "Cannot find Percolator from %s.", percolatorPath));
        }

        if (!(new File(percolatorPath)).canExecute()) {
            throw new Exception(String.format(Locale.US, "Percolator (%s) exits but cannot be executed.", percolatorPath));
        }

        String[] tempArray = parameterMap.get("ms_level").split(",");
        Set<Integer> msLevelSet = new HashSet<>(tempArray.length + 1, 1);
        for (String temp : tempArray) {
            msLevelSet.add(Integer.valueOf(temp));
        }

        String labelling = "N14";
        if (parameterMap.get("15N").trim().contentEquals("1")) {
            labelling = "N15";
        }

        logger.info("Loading parameters and database...");
        BuildIndex buildIndex = new BuildIndex(parameterMap);
        MassTool massTool = buildIndex.returnMassTool();
        System.out.println("lsz db length "+ buildIndex.getPeptide0Map().size());
        InferPTM inferPTM = buildIndex.getInferPTM();

        Set<String> tagPepSize = new HashSet<>();
//        for(Set<String> pepset : buildIndex.getInferSegment().tagPepMap.values()) {
//            tagPepSize.addAll(pepset);
//        }
//        System.out.println("lsz tagMap length "+ tagPepSize.size());
//        System.out.println(buildIndex.getInferSegment().tagPepMap.keySet());


        logger.info("Reading spectra...");
        File spectraFile = new File(spectraPath);
        DatasetReader datasetReader;
        JMzReader[] spectraParserArray;
        String sqlPath = "jdbc:sqlite:" + dbName;
        Class.forName("org.sqlite.JDBC").newInstance();
        Map<Integer, String> fileIdNameMap = new HashMap<>();
        if ((!spectraFile.exists())) {
            throw new FileNotFoundException("The spectra file not found.");
        }

        if ( ! spectraFile.isDirectory()) {
            spectraParserArray = new JMzReader[1];
            JMzReader spectraParser;
            String ext = spectraPath.substring(spectraPath.lastIndexOf(".")+1);
            if (ext.contentEquals("mzXML")) {
                spectraParser = new MzXMLFile(spectraFile);
            } else if (ext.toLowerCase().contentEquals("mgf")) {
                spectraParser = new MgfFile(spectraFile);
            } else {
                throw new Exception(String.format(Locale.US, "Unsupported file format %s. Currently, PIPI only support mzXML and MGF.", ext));
            }
            spectraParserArray[0] = spectraParser;
            fileIdNameMap.put(0, spectraPath.substring(spectraPath.lastIndexOf("/")+1));
            datasetReader = new DatasetReader(spectraParserArray, ms1Tolerance, ms1ToleranceUnit, massTool, ext, msLevelSet, sqlPath);
        } else {
            String[] fileList = spectraFile.list(new FilenameFilter() {
                @Override
                public boolean accept(File dir, String name) {
                    return name.endsWith(".mgf");
                }
            });
            spectraParserArray = new JMzReader[fileList.length];
            for (int i = 0; i < fileList.length; i++){
                spectraParserArray[i] = new MgfFile(new File(spectraPath + fileList[i]));
                fileIdNameMap.put(i, fileList[i]);
            }

            String ext = fileList[0].substring(fileList[0].lastIndexOf(".")+1);
            datasetReader = new DatasetReader(spectraParserArray, ms1Tolerance, ms1ToleranceUnit, massTool, ext, msLevelSet, sqlPath);
        }
//        BufferedReader parameterReader = new BufferedReader(new FileReader("/home/slaiad/Data/PXD004732/pool121/truth.txt"));
        BufferedReader parameterReader = new BufferedReader(new FileReader("/home/slaiad/Code/PIPI/src/main/resources/ChickOpenTruth.txt"));

        Map<Integer, String> pepTruth = new HashMap<>();
        String line;
        while ((line = parameterReader.readLine()) != null) {
            if (line.isEmpty()) continue;
            line = line.trim();
            String[] splitRes = line.split(",");
            if (splitRes.length == 1) {
                int a = 1;
            }
            String pepWithMod = splitRes[1].replace('L','#').replace('I', '#');
            pepTruth.put(Integer.valueOf(splitRes[0]), pepWithMod);
        }

        BufferedReader parameterReader1 = new BufferedReader(new FileReader("/home/slaiad/Data/PXD001468/omProtUnion.txt"));
        Set<String> omProts = new HashSet<>();
        String line1;
        while ((line1 = parameterReader1.readLine()) != null) {
            if (line1.isEmpty()) continue;
            line1 = line1.trim();
            omProts.add(line1);
        }
        SpecProcessor specProcessor = new SpecProcessor(massTool);
        Map<String, Integer> precursorChargeMap = new HashMap<>();
        Map<String, Double> precursorMassMap = new HashMap<>();
        TreeMap<Double, Set<String>> pcMassScanNameMap = new TreeMap<>();

        logger.info("Scan decoding and protein db reducing...");
        int threadNum1 = Integer.valueOf(parameterMap.get("thread_num"));
        if (threadNum1 == 0) {
            threadNum1 = 3 + Runtime.getRuntime().availableProcessors();
        }
        if (java.lang.management.ManagementFactory.getRuntimeMXBean().getInputArguments().toString().indexOf("jdwp") >= 0){
//            threadNum1 = 1;
        }
        System.out.println("thread NUM "+ threadNum1);

        Map<String, String> mgfTitleMap = new HashMap<>();
        Map<String, Integer> isotopeCorrectionNumMap = new HashMap<>();
        Map<String, Double> ms1PearsonCorrelationCoefficientMap = new HashMap<>();

        ExecutorService threadPoolSpecCoder = Executors.newFixedThreadPool(threadNum1);
        ArrayList<Future<SpecCoder.Entry>> taskListSpecCoder = new ArrayList<>(datasetReader.getUsefulSpectraNum() + 10);
        Connection sqlConSpecCoder = DriverManager.getConnection(sqlPath);
        Statement sqlStatementSpecCoder = sqlConSpecCoder.createStatement();
        ResultSet sqlResSetSpecCoder = sqlStatementSpecCoder.executeQuery("SELECT scanName, scanNum, precursorCharge" +
                ", precursorMass, precursorScanNo, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient FROM spectraTable");

        ReentrantLock lockSpecCoder = new ReentrantLock();

        int submitNumSpecCoder = 0;

        while (sqlResSetSpecCoder.next()) {
            String scanName = sqlResSetSpecCoder.getString("scanName");
            int scanNum = sqlResSetSpecCoder.getInt("scanNum");
            int precursorCharge = sqlResSetSpecCoder.getInt("precursorCharge");
            double precursorMass = sqlResSetSpecCoder.getDouble("precursorMass");
            mgfTitleMap.put(scanName, sqlResSetSpecCoder.getString("mgfTitle"));
            isotopeCorrectionNumMap.put(scanName, sqlResSetSpecCoder.getInt("isotopeCorrectionNum"));
            ms1PearsonCorrelationCoefficientMap.put(scanName, sqlResSetSpecCoder.getDouble("ms1PearsonCorrelationCoefficient"));

//            if (scanNum > 2500) {  //22459
//                continue;
//            }l
//            if (scanName.split("\\.")[1].contentEquals("25") ) {  //22459
//                continue;
//            }
            int fileId = Integer.valueOf( scanName.split("\\.")[0] );

            precursorChargeMap.put(scanName, precursorCharge);
            precursorMassMap.put(scanName, precursorMass);

            taskListSpecCoder.add(threadPoolSpecCoder.submit(new SpecCoder(scanNum, buildIndex, massTool, spectraParserArray[fileId], minClear, maxClear, lockSpecCoder, scanName, precursorCharge
                    , precursorMass, specProcessor, pepTruth.get(1886))));
        }
        System.out.println("totalSubmit in SpecCoder, "+ submitNumSpecCoder);
        sqlResSetSpecCoder.close();
        sqlStatementSpecCoder.close();
//        sqlConSpecCoder.close();

//        Map<String, SparseVector> specCodeMap = new HashMap<>();
        Set<String> validScanSet = new HashSet<>();
        int lastProgressSpecCoder = 0;
//        int resultCountSpecCoder = 0;
        int totalCountSpecCoder = taskListSpecCoder.size();
        int countSpecCoder = 0;
        Map<String, Map<String, Double>> protTagScoreMapMap = new HashMap<>();
        while (countSpecCoder < totalCountSpecCoder) {
            // record search results and delete finished ones.
            List<Future<SpecCoder.Entry>> toBeDeleteTaskList = new ArrayList<>(totalCountSpecCoder - countSpecCoder);
            for (Future<SpecCoder.Entry> task : taskListSpecCoder) {
                if (task.isDone()) {
                    if (task.get() != null ) {
                        SpecCoder.Entry entry = task.get();
                        validScanSet.add(entry.scanName);

                        int coId = Integer.valueOf(entry.scanName.split("\\.")[3]);
                        if (coId > 0) {
                            toBeDeleteTaskList.add(task);
                            continue;
                        }


                        for (String prot : entry.protTagScoreMap.keySet()){
                            List<Pair<String, Double>> tagScoreList = entry.protTagScoreMap.get(prot);
                            if (protTagScoreMapMap.containsKey(prot)) {
                                Map<String, Double> tagScoreMap =  protTagScoreMapMap.get(prot);
                                for (Pair<String, Double> tagScore : tagScoreList){
                                    String tag = tagScore.getFirst();
                                    if (tagScoreMap.containsKey(tag)){
                                        tagScoreMap.put(tag, Math.min(tag.length()*2, tagScore.getSecond() + tagScoreMap.get(tag)));
                                    } else {
                                        tagScoreMap.put(tag, tagScore.getSecond());
                                    }
                                }
                            } else {
                                Map<String, Double> tagScoreMap =  new HashMap<>();
                                for (Pair<String, Double> tagScore : tagScoreList){
                                    String tag = tagScore.getFirst();
                                    if (tagScoreMap.containsKey(tag)){
                                        tagScoreMap.put(tag, Math.min(tag.length()*2, tagScore.getSecond() + tagScoreMap.get(tag)));
                                    } else {
                                        tagScoreMap.put(tag, tagScore.getSecond());
                                    }
                                }
                                protTagScoreMapMap.put(prot, tagScoreMap);
                            }
                        }

                    }

                    toBeDeleteTaskList.add(task);
                } else if (task.isCancelled()) {
                    toBeDeleteTaskList.add(task);
                }
            }
            countSpecCoder += toBeDeleteTaskList.size();
            taskListSpecCoder.removeAll(toBeDeleteTaskList);
            taskListSpecCoder.trimToSize();

            int progress = countSpecCoder * 20 / totalCountSpecCoder;
            if (progress != lastProgressSpecCoder) {
                logger.info("Spectra Coding {}%...", progress * 5);
                lastProgressSpecCoder = progress;
            }

            if (countSpecCoder == totalCountSpecCoder) {
                break;
            }
            Thread.sleep(6000);
        }
        // shutdown threads.
        threadPoolSpecCoder.shutdown();
        if (!threadPoolSpecCoder.awaitTermination(60, TimeUnit.SECONDS)) {
            threadPoolSpecCoder.shutdownNow();
            if (!threadPoolSpecCoder.awaitTermination(60, TimeUnit.SECONDS))
                throw new Exception("Pool did not terminate");
        }

        if (lockSpecCoder.isLocked()) {
            lockSpecCoder.unlock();
        }

        Map<String, Integer> protLengthMap = buildIndex.protLengthMap;
        Map<String, Double> protScoreFinalMap = new HashMap<>();
        Map<String, Integer> tagNumMap = new HashMap<>();
        for (String prot : protTagScoreMapMap.keySet()){
            double score = 0;
            double totalLen = 0;
            Map<String , Double> tagScoreMap = protTagScoreMapMap.get(prot);
            if (tagScoreMap.isEmpty()) {
                continue;
            }
            for  (String tag :tagScoreMap.keySet()) {
                score += tagScoreMap.get(tag) * tag.length()* tag.length();
                totalLen+= tag.length();
            }
            protScoreFinalMap.put(prot, score*totalLen/Math.sqrt(protLengthMap.get(prot)));

            for (String tag :tagScoreMap.keySet()) {
                if (tagNumMap.containsKey(tag)) {
                    tagNumMap.put(tag, tagNumMap.get(tag)+1);
                } else {
                    tagNumMap.put(tag, 1);
                }
            }
        }

        List<Pair<String, Double>> protScoreLongList = new ArrayList<>();
        for (String prot : protScoreFinalMap.keySet()){
            protScoreLongList.add(new Pair<>(prot, protScoreFinalMap.get(prot)));
        }
        protScoreLongList.sort(Comparator.comparingDouble(Pair::getSecond));
        int num = 0;
        int correct = 0;
        int totalCorrect = 0;
        for (int i = protScoreLongList.size()-1; i >= 0; i--){
//            System.out.print(protScoreLongList.get(i).getFirst()+","+(omProts.contains(protScoreLongList.get(i).getFirst()) ? 1:0)+","+protScoreLongList.get(i).getSecond() +","
//                    + protLengthMap.get(protScoreLongList.get(i).getFirst())+","+protTagNumMap.get(protScoreLongList.get(i).getFirst()));
//            System.out.print("\n");
            num++;
            if (num < 5300 && omProts.contains(protScoreLongList.get(i).getFirst()) ) {
                correct++;
            }
            if (omProts.contains(protScoreLongList.get(i).getFirst())) {
                totalCorrect ++;

            }
            System.out.println(num+ "," +protScoreLongList.get(i).getSecond() + ","+protLengthMap.get(protScoreLongList.get(i).getFirst() )
                    + ","+protScoreLongList.get(i).getFirst() +"," + (omProts.contains(protScoreLongList.get(i).getFirst()) ? 1:0));
            for (String tag : protTagScoreMapMap.get(protScoreLongList.get(i).getFirst()).keySet()){
                System.out.println(tag+ "," + protTagScoreMapMap.get(protScoreLongList.get(i).getFirst()).get(tag) + "," +tagNumMap.get(tag));
            }
            System.out.println("==============");
        }
        System.out.println("5300 contains,"+ correct);
        System.out.println("total,"+ num + ","+totalCorrect);
        Set<String> protHardSet = new HashSet<>();
        for (Pair<String, Double> pair : protScoreLongList){
            if (pair.getSecond() < 59) continue;
            protHardSet.add(pair.getFirst());
        }
//        Set<String> protHardSet = omProts;

        Map<String, Peptide0> peptide0Map = buildIndex.getPeptide0Map();
        Iterator<Map.Entry<String, Peptide0>> iter = peptide0Map.entrySet().iterator();
        while (iter.hasNext()) {
            boolean shouldKeep = false;
            for (String prot : iter.next().getValue().proteins) {
                if (prot.contains("DECOY_")) prot = prot.substring(6);
                if (protHardSet.contains(prot)){
                    shouldKeep = true;
                    break;
                }
            }
            if (!shouldKeep) {
                iter.remove();
            }
        }
        System.out.println("reduced massPepMap size," + peptide0Map.size());

        TreeMap<Double, Set<String>> massPeptideMap = buildIndex.getMassPeptideMap();
        for (double mass : massPeptideMap.keySet()) {
            Set<String> pepSet = massPeptideMap.get(mass);
            Iterator<String> iterPep = pepSet.iterator();
            while (iterPep.hasNext()) {
                if (!peptide0Map.containsKey(iterPep.next())) {
                    iterPep.remove();
                }
            }
        }
        System.out.println("reduced massPepMap size," + massPeptideMap.size());
        int numTargetPep = 0, numDecoyPep = 0;
        for (Peptide0 pep : peptide0Map.values()){
            if (pep.isTarget) {
                numTargetPep++;
            }else {
                numDecoyPep++;
            }
        }
        System.out.println("Reduced db size "+(numDecoyPep+ numTargetPep)+",targer : decoy = "+numTargetPep+" : "+numDecoyPep);

        System.out.println("prot Db," + protScoreLongList.size()+","+protHardSet.size());

        logger.info("Pre searching...");
        int threadNum = Integer.valueOf(parameterMap.get("thread_num"));
        if (threadNum == 0) {
            threadNum = 3 + Runtime.getRuntime().availableProcessors();
        }
        if (java.lang.management.ManagementFactory.getRuntimeMXBean().getInputArguments().toString().indexOf("jdwp") >= 0){
//            threadNum = 1;
        }
        System.out.println("thread NUM "+ threadNum);
        ExecutorService threadPoolBone = Executors.newFixedThreadPool(threadNum);
        ArrayList<Future<PreSearch.Entry>> taskListBone = new ArrayList<>(validScanSet.size() + 10);
        ReentrantLock lockBone = new ReentrantLock();
        int submitNumBone = 0;
        for (String scanName : validScanSet) {
            String[] scanNameStr = scanName.split("\\.");
            int precursorCharge = precursorChargeMap.get(scanName);

            double precursorMass = precursorMassMap.get(scanName);
            submitNumBone++;
            taskListBone.add(threadPoolBone.submit(new PreSearch(Integer.valueOf(scanNameStr[2]), buildIndex, massTool, ms1Tolerance, leftInverseMs1Tolerance, rightInverseMs1Tolerance
                    , ms1ToleranceUnit, inferPTM.getMinPtmMass(), inferPTM.getMaxPtmMass(), Math.min(precursorCharge > 1 ? precursorCharge - 1 : 1, 3)
                    , spectraParserArray[Integer.valueOf(scanNameStr[0])], minClear, maxClear, lockBone, scanName, precursorCharge, precursorMass, specProcessor, pepTruth.get(1886))));
        }
        System.out.println("totalSubmit in Bone, "+ submitNumBone);

        Map<String, List<Peptide>> ptmOnlyCandiMap = new HashMap<>();
        Map<String, List<Peptide>> ptmFreeCandiMap = new HashMap<>();

        Map<Integer, Map<Integer, Set<Integer>>> coGroupsMap = new HashMap<>();
        Map<Integer, Map<Integer, Map<Integer, Integer>>> fileIdScanNumChargeTimeMap = new HashMap<>();
        for (int fileId = 0; fileId < fileIdNameMap.size(); fileId++) {
            fileIdScanNumChargeTimeMap.put(fileId, new HashMap<>());
        }
        Map<String, Map<Integer, Set<String>>> fileIdScanNumToChargeScanNameMap = new HashMap<>(); //key: fileId.scanNum, value: fileId.scanId.scanNum.coId
        int lastProgressBone = 0;
        int resultCountBone = 0;
        int totalCountBone = taskListBone.size();
        int countBone = 0;
        while (countBone < totalCountBone) {
            // record search results and delete finished ones.
            List<Future<PreSearch.Entry>> toBeDeleteTaskList = new ArrayList<>(totalCountBone - countBone);
            for (Future<PreSearch.Entry> task : taskListBone) {
                if (task.isDone()) {
                    if (task.get() != null) {
                        PreSearch.Entry entry = task.get();

                        ptmOnlyCandiMap.put(entry.scanName, entry.ptmOnlyList);
                        ptmFreeCandiMap.put(entry.scanName, entry.ptmFreeList);
                        if (pcMassScanNameMap.containsKey(entry.precursorMass)) {
                            pcMassScanNameMap.get(entry.precursorMass).add(entry.scanName);
                        }else {
                            Set<String> scanNumSet = new HashSet<>();
                            scanNumSet.add(entry.scanName);
                            pcMassScanNameMap.put(entry.precursorMass, scanNumSet);
                        }

                        int charge = precursorChargeMap.get(entry.scanName);
                        String[] scanNameStr = entry.scanName.split("\\.");
                        String fileIdScanNum = scanNameStr[0] + "." + scanNameStr[2];
                        if (fileIdScanNumToChargeScanNameMap.containsKey(fileIdScanNum)) {
                            Map<Integer, Set<String>> chargeScanNameMap = fileIdScanNumToChargeScanNameMap.get(fileIdScanNum);
                            if (chargeScanNameMap.containsKey(charge)){
                                chargeScanNameMap.get(charge).add(entry.scanName);
                            } else {
                                Set<String> scanNameSet = new HashSet<>();
                                scanNameSet.add(entry.scanName);
                                chargeScanNameMap.put(charge,scanNameSet);
                            }
                        } else {
                            Map<Integer, Set<String>> chargeScanNameMap = new HashMap<>();
                            if (chargeScanNameMap.containsKey(charge)){
                                chargeScanNameMap.get(charge).add(entry.scanName);
                            } else {
                                Set<String> scanNameSet = new HashSet<>();
                                scanNameSet.add(entry.scanName);
                                chargeScanNameMap.put(charge,scanNameSet);
                            }
                            fileIdScanNumToChargeScanNameMap.put(fileIdScanNum, chargeScanNameMap);
                        }

                        ++resultCountBone;
                    }

                    toBeDeleteTaskList.add(task);
                } else if (task.isCancelled()) {
                    toBeDeleteTaskList.add(task);
                }
            }
            countBone += toBeDeleteTaskList.size();
            taskListBone.removeAll(toBeDeleteTaskList);
            taskListBone.trimToSize();

            int progress = countBone * 20 / totalCountBone;
            if (progress != lastProgressBone) {
                logger.info("PreSearching {}%...", progress * 5);
                lastProgressBone = progress;
            }

            if (countBone == totalCountBone) {
                break;
            }
            Thread.sleep(6000);
        }

        // shutdown threads.
        threadPoolBone.shutdown();
        if (!threadPoolBone.awaitTermination(60, TimeUnit.SECONDS)) {
            threadPoolBone.shutdownNow();
            if (!threadPoolBone.awaitTermination(60, TimeUnit.SECONDS))
                throw new Exception("Pool did not terminate");
        }

        if (lockBone.isLocked()) {
            lockBone.unlock();
        }

        System.out.println("resultCount," +resultCountBone);


        Iterator<Map.Entry<String, Map<Integer, Set<String>>>> iter1 = fileIdScanNumToChargeScanNameMap.entrySet().iterator();
        while (iter1.hasNext()) {
            Map<Integer, Set<String>> chargeScanNameMap = iter1.next().getValue();
            Iterator<Map.Entry<Integer, Set<String>>> iter2 = chargeScanNameMap.entrySet().iterator();

            while (iter2.hasNext()) {
                if (iter2.next().getValue().size() <= 1) {
                    iter2.remove();
                }
            }
            if (chargeScanNameMap.isEmpty()) {
                iter1.remove();
            }
        }

        int numScanNamesToCo = 0;
        for (String fileIdScanNum : fileIdScanNumToChargeScanNameMap.keySet()){
            Map<Integer, Set<String>> chargeScanNameMap = fileIdScanNumToChargeScanNameMap.get(fileIdScanNum);
            for (int charge : chargeScanNameMap.keySet()) {
                numScanNamesToCo += chargeScanNameMap.get(charge).size();
            }
        }
        System.out.println("numScanNamesToCo +" + numScanNamesToCo);


        System.out.println("lsz +" +","+ptmOnlyCandiMap.size()+","+pcMassScanNameMap.size() + "," + ptmOnlyCandiMap.keySet().size());
        logger.info("Start searching...");
        if (threadNum == 0) {
            threadNum = 3 + Runtime.getRuntime().availableProcessors();
        }
        if (java.lang.management.ManagementFactory.getRuntimeMXBean().getInputArguments().toString().indexOf("jdwp") >= 0){
//            threadNum = 1;
        }
        ExecutorService threadPoolPtm = Executors.newFixedThreadPool(threadNum);
        ArrayList<Future<PtmSearch.Entry>> taskListPTM = new ArrayList<>(ptmOnlyCandiMap.keySet().size() + 10);
        ReentrantLock lockPtm = new ReentrantLock();
        Binomial binomial = new Binomial(Integer.valueOf(parameterMap.get("max_peptide_length")) * 2);
        int submitTimePtm = 0;


        Map<String, Peptide0> pep0Map = buildIndex.getPeptide0Map();
        for (double mass : pcMassScanNameMap.keySet()) {
            for (String thisScanName : pcMassScanNameMap.get(mass)){
                int thisScanNum = Integer.valueOf( thisScanName.split("\\.")[2] );
                int thisFileId = Integer.valueOf( thisScanName.split("\\.")[0] );
                Set<Peptide> realPtmOnlyList = new HashSet<>();
                Set<Peptide> realPtmFreeList = new HashSet<>();
                for (Peptide pep : ptmOnlyCandiMap.get(thisScanName)) {
//                    for (String prot : pep0Map.get(pep.getPTMFreePeptide()).proteins) {
//                        if (prot.contains("DECOY_")) prot = prot.substring(6);
//
//                        if (protHardSet.contains(prot)) {
                            realPtmOnlyList.add(pep.clone());
//                            break;
//                        }
//                    }
                }
                for (Peptide pep : ptmFreeCandiMap.get(thisScanName)) {
//                    for (String prot : pep0Map.get(pep.getPTMFreePeptide()).proteins) {
//                        if (prot.contains("DECOY_")) prot = prot.substring(6);
//                        if (protHardSet.contains(prot)) {
                            realPtmFreeList.add(pep.clone());
//                            break;
//                        }
//                    }
                }
//                System.out.println(thisScanNum+","+ptmOnlyCandiMap.get(thisScanNum).size()+","+ptmFreeCandiMap.get(thisScanNum).size()+"," + realPtmOnlyList.size()+","+realPtmFreeList.size());
                for (int i = -3; i <= 3; i++) {
//                    if (i == 0) continue;

                    for (Set<String> otherScanNameSet : pcMassScanNameMap.subMap(mass + i*MassTool.PROTON - 0.02, true, mass + i*MassTool.PROTON + 0.02, true).values()) {
                        for (String otherScanName : otherScanNameSet) {
                            int otherScanNum = Integer.valueOf( otherScanName.split("\\.")[2] );
                            int otherFileId = Integer.valueOf( otherScanName.split("\\.")[0] );
                            if (otherFileId != thisFileId) continue;

                            if (otherScanNum < thisScanNum+2000 && otherScanNum > thisScanNum-2000 && otherScanNum != thisScanNum) {
                                for (Peptide pep : ptmOnlyCandiMap.get(otherScanName)) {
//                                    for (String prot : pep0Map.get(pep.getPTMFreePeptide()).proteins) {
//                                        if (prot.contains("DECOY_")) prot = prot.substring(6);
//                                        if (protHardSet.contains(prot)) {
                                            realPtmOnlyList.add(pep.clone());
//                                            break;
//                                        }
//                                    }
                                }
                                for (Peptide pep : ptmFreeCandiMap.get(otherScanName)) {
//                                    for (String prot : pep0Map.get(pep.getPTMFreePeptide()).proteins) {
//                                        if (prot.contains("DECOY_")) prot = prot.substring(6);
//                                        if (protHardSet.contains(prot)) {
                                            realPtmFreeList.add(pep.clone());
//                                            break;
//                                        }
//                                    }
                                }
                            }
                        }
                    }
                }

                int precursorCharge = precursorChargeMap.get(thisScanName);
                int precursorScanNo = 0;
                double precursorMass = precursorMassMap.get(thisScanName);
                submitTimePtm++;
                boolean hasExtra = false;
                if (fileIdScanNumToChargeScanNameMap.containsKey(thisFileId + "."+ thisScanNum) && fileIdScanNumToChargeScanNameMap.get(thisFileId + "."+ thisScanNum).containsKey(precursorCharge)) {
                    hasExtra = true;
                }
                taskListPTM.add(threadPoolPtm.submit(new PtmSearch(thisScanNum, buildIndex, massTool, ms1Tolerance
                        , leftInverseMs1Tolerance, rightInverseMs1Tolerance, ms1ToleranceUnit, ms2Tolerance
                        , inferPTM.getMinPtmMass(), inferPTM.getMaxPtmMass(), Math.min(precursorCharge > 1 ? precursorCharge - 1 : 1, 3)
                        , spectraParserArray[thisFileId], minClear, maxClear, lockPtm, thisScanName, precursorCharge
                        , precursorMass, inferPTM, specProcessor, sqlPath, binomial, precursorScanNo, realPtmOnlyList, realPtmFreeList, hasExtra)));
            }
        }
        System.out.println("submit times in Ptm" + submitTimePtm);
        // check progress every minute, record results,and delete finished tasks.
        PreparedStatement sqlPreparedStatement = sqlConSpecCoder.prepareStatement("REPLACE INTO spectraTable (scanNum, scanName,  precursorCharge, precursorMass, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient, labelling, peptide, theoMass, isDecoy, globalRank, normalizedCorrelationCoefficient, score, deltaLCn, deltaCn, matchedPeakNum, ionFrac, matchedHighestIntensityFrac, explainedAaFrac, otherPtmPatterns, aScore, candidates, peptideSet, whereIsTopCand, shouldPtm, hasPTM, ptmNum, isSettled) VALUES (?, ?, ?, ?,?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)");
        sqlConSpecCoder.setAutoCommit(false);
        int lastProgressPtm = 0;
        int resultCountPtm = 0;
        int totalCountPtm = taskListPTM.size();
        int countPtm = 0;

        Map<String, List<ExRes>> exResForScanNameMap = new HashMap<>();
//        System.out.println("lsz totalCount, " + totalCount);
        while (countPtm < totalCountPtm) {
            // record search results and delete finished ones.
            List<Future<PtmSearch.Entry>> toBeDeleteTaskList = new ArrayList<>(totalCountPtm - countPtm);
            for (Future<PtmSearch.Entry> task : taskListPTM) {
                if (task.isDone()) {
                    if (task.get() != null) {
                        PtmSearch.Entry entry = task.get();
                        sqlPreparedStatement.setInt(1, entry.scanNum);
                        sqlPreparedStatement.setString(2, entry.scanName);
                        sqlPreparedStatement.setInt(3, entry.precursorCharge);
                        sqlPreparedStatement.setDouble(4, entry.precursorMass);
                        sqlPreparedStatement.setString(5, mgfTitleMap.get(entry.scanName));
                        sqlPreparedStatement.setInt(6, isotopeCorrectionNumMap.get(entry.scanName));
                        sqlPreparedStatement.setDouble(7, ms1PearsonCorrelationCoefficientMap.get(entry.scanName));
                        sqlPreparedStatement.setString(8, entry.labelling);
                        sqlPreparedStatement.setString(9, entry.peptide);
                        sqlPreparedStatement.setDouble(10, entry.theoMass);
                        sqlPreparedStatement.setInt(11, entry.isDecoy);
                        sqlPreparedStatement.setInt(12, entry.globalRank);
                        sqlPreparedStatement.setDouble(13, entry.normalizedCorrelationCoefficient);
                        sqlPreparedStatement.setDouble(14, entry.score);
                        sqlPreparedStatement.setDouble(15, entry.deltaLCn);
                        sqlPreparedStatement.setDouble(16, entry.deltaCn);
                        sqlPreparedStatement.setInt(17, entry.matchedPeakNum);
                        sqlPreparedStatement.setDouble(18, entry.ionFrac);
                        sqlPreparedStatement.setDouble(19, entry.matchedHighestIntensityFrac);
                        sqlPreparedStatement.setDouble(20, entry.explainedAaFrac);
                        sqlPreparedStatement.setString(21, entry.otherPtmPatterns);
                        sqlPreparedStatement.setString(22, entry.aScore);
                        sqlPreparedStatement.setString(23, entry.candidates);
                        sqlPreparedStatement.setString(24, entry.peptideSet);
                        sqlPreparedStatement.setInt(25, entry.whereIsTopCand);
                        sqlPreparedStatement.setInt(26, entry.shouldPtm);
                        sqlPreparedStatement.setInt(27, entry.hasPTM);
                        sqlPreparedStatement.setInt(28, entry.ptmNum);
                        sqlPreparedStatement.setInt(29, entry.isSettled);

                        if (!entry.extraEntryList.isEmpty()) {
                            List<ExRes> exResList = new ArrayList<>();
                            for (PtmSearch.ExtraEntry extraEntry : entry.extraEntryList) {
                                ExRes res = new ExRes();
                                res.scanNum = extraEntry.scanNum;
                                res.scanName = extraEntry.scanName;
                                res.precursorCharge = extraEntry.precursorCharge;
                                res.precursorMass = extraEntry.precursorMass;
                                res.mgfTitle = mgfTitleMap.get(entry.scanName);
                                res.isotopeCorrectionNum = isotopeCorrectionNumMap.get(entry.scanName);
                                res.ms1PearsonCorrelationCoefficient = ms1PearsonCorrelationCoefficientMap.get(entry.scanName);
                                res.labelling = extraEntry.labelling;
                                res.peptide = extraEntry.peptide;
                                res.theoMass = extraEntry.theoMass;
                                res.isDecoy = extraEntry.isDecoy;
                                res.globalRank = extraEntry.globalRank;
                                res.normalizedCorrelationCoefficient = extraEntry.normalizedCorrelationCoefficient;
                                res.score = extraEntry.score;
                                res.deltaLCn = extraEntry.deltaLCn;
                                res.deltaCn = extraEntry.deltaCn;
                                res.matchedPeakNum = extraEntry.matchedPeakNum;
                                res.ionFrac = extraEntry.ionFrac;
                                res.matchedHighestIntensityFrac = extraEntry.matchedHighestIntensityFrac;
                                res.explainedAaFrac = extraEntry.explainedAaFrac;
                                res.otherPtmPatterns = extraEntry.otherPtmPatterns;
                                res.aScore = extraEntry.aScore;
                                res.candidates = extraEntry.candidates;
                                res.peptideSet = extraEntry.peptideSet;
                                res.whereIsTopCand = extraEntry.whereIsTopCand;
                                res.shouldPtm = extraEntry.shouldPtm;
                                res.hasPTM = extraEntry.hasPTM;
                                res.ptmNum = extraEntry.ptmNum;
                                res.isSettled = extraEntry.isSettled;
                                exResList.add(res);
                            }
                            exResForScanNameMap.put(entry.scanName, exResList);
                        }


                        sqlPreparedStatement.executeUpdate();
                        ++resultCountPtm;
                    }

                    toBeDeleteTaskList.add(task);
                } else if (task.isCancelled()) {
                    toBeDeleteTaskList.add(task);
                }
            }
            countPtm += toBeDeleteTaskList.size();
            taskListPTM.removeAll(toBeDeleteTaskList);
            taskListPTM.trimToSize();

            sqlConSpecCoder.commit();

            int progress = countPtm * 20 / totalCountPtm;
            if (progress != lastProgressPtm) {
                logger.info("Searching {}%...", progress * 5);
//                System.out.println(toBeDeleteTaskList.size()+ ","+ taskListPTM.size() +"," + count2 + ", " + totalCount + "," + lastProgress);
                lastProgressPtm = progress;
            }

            if (countPtm == totalCountPtm) {
                break;
            }
            Thread.sleep(6000);
        }
        System.out.println("exResForScanNameMap, "+exResForScanNameMap.size());
        // shutdown threads.
        threadPoolPtm.shutdown();
        if (!threadPoolPtm.awaitTermination(60, TimeUnit.SECONDS)) {
            threadPoolPtm.shutdownNow();
            if (!threadPoolPtm.awaitTermination(60, TimeUnit.SECONDS))
                throw new Exception("Pool did not terminate");
        }

        sqlConSpecCoder.commit();
        sqlConSpecCoder.setAutoCommit(true);
        sqlConSpecCoder.close();
        if (lockPtm.isLocked()) {
            lockPtm.unlock();
        }

        if (resultCountPtm == 0) {
            throw new Exception("There is no useful results.");
        }

        String percolatorInputFileName = spectraPath + "." + labelling + ".input.temp";
        writePercolator(percolatorInputFileName, buildIndex.getPeptide0Map(), sqlPath, fileIdScanNumToChargeScanNameMap, exResForScanNameMap);
        Map<String, PercolatorEntry> percolatorResultMap = null;

        if (parameterMap.get("add_decoy").contentEquals("0")) {
            logger.warn("add_decoy = 0. Don't estimate FDR.");
        } else {
            logger.info("Estimating FDR...");
            String percolatorOutputFileName = spectraPath + "." + labelling + ".output.temp";
            String percolatorProteinOutputFileName = spectraPath + "." + labelling + ".protein.tsv";
            percolatorResultMap = runPercolator(percolatorPath, percolatorInputFileName, percolatorOutputFileName, percolatorProteinOutputFileName, parameterMap.get("db") + ".TD.fasta", parameterMap.get("enzyme_name_1"));
            if (percolatorResultMap.isEmpty()) {
                throw new Exception(String.format(Locale.US, "Percolator failed to estimate FDR. Please check if Percolator is installed and the percolator_path in %s is correct.", parameterPath));
            }
            if (!outputPercolatorInput) {
                (new File(percolatorOutputFileName)).delete();
            }
        }

        if (!outputPercolatorInput) {
            (new File(percolatorInputFileName)).delete();
        }

        logger.info("Saving results...");
        writeTop20Candidates(sqlPath, spectraPath);

        writeFinalResult(percolatorResultMap, spectraPath + "." + labelling + ".pipi.csv", buildIndex.getPeptide0Map(), sqlPath);
//        new WritePepXml(spectraPath + "." + labelling + ".pipi.pep.xml", spectraPath, parameterMap, massTool.getMassTable(), percolatorResultMap, buildIndex.getPeptide0Map(), buildIndex.returnFixModMap(), sqlPath);
    }

    class ExRes {
        ExRes (){};
        public int scanNum = 0;
        public String scanName = "";
        public int precursorCharge = 0;
        public double precursorMass = 0d;
        public String mgfTitle = "";
        public int isotopeCorrectionNum = 0;
        public double ms1PearsonCorrelationCoefficient = 0d;
        public String labelling = "";
        public String peptide = "";
        public double theoMass = 0d;
        public int isDecoy = 0;
        public int globalRank = 0;
        public double normalizedCorrelationCoefficient = 0d;
        public double score = 0d;
        public double deltaLCn = 0d;
        public double deltaCn = 0d;
        public int matchedPeakNum = 0;
        public double ionFrac = 0d;
        public double matchedHighestIntensityFrac = 0d;
        public double explainedAaFrac = 0d;
        public String otherPtmPatterns = "";
        public String aScore = "";
        public String candidates = "";
        public String peptideSet = "";
        public int whereIsTopCand = 0;
        public int shouldPtm = 0;
        public int hasPTM = 0;
        public int ptmNum = 0;
        public int isSettled = 0;
    }
    private static void help() {
        String helpStr = "PIPI version " + versionStr + "\r\n"
                + "A tool identifying peptides with unlimited PTM.\r\n"
                + "Author: Fengchao Yu\r\n"
                + "Email: fyuab@connect.ust.hk\r\n"
                + "PIPI usage: java -Xmx25g -jar /path/to/PIPI.jar <parameter_file> <data_file>\r\n"
                + "\t<parameter_file>: parameter file. Can be download along with PIPI.\r\n"
                + "\t<data_file>: spectra data file (mzXML)\r\n"
                + "\texample: java -Xmx32g -jar PIPI.jar parameter.def data.mzxml\r\n";
        System.out.print(helpStr);
        System.exit(1);
    }
    private void writePercolator(String resultPath, Map<String, Peptide0> peptide0Map, String sqlPath, Map<String, Map<Integer, Set<String>>> fileIdScanNumToChargeScanNameMap, Map<String, List<ExRes>> exResForScanNameMap) throws IOException, SQLException {
        BufferedWriter writer = new BufferedWriter(new FileWriter(resultPath));
        writer.write("id\tlabel\tscannr\texpmass\tcalcmass\tscore\tdelta_c_n\tdelta_L_c_n\tnormalized_cross_corr\tglobal_search_rank\tabs_ppm\tion_frac\tmatched_high_peak_frac\tcharge1\tcharge2\tcharge3\tcharge4\tcharge5\tcharge6\texplained_aa_frac\tpeptide\tprotein\n");
        Connection sqlConnection = DriverManager.getConnection(sqlPath);
        Statement sqlStatement = sqlConnection.createStatement();
        ResultSet sqlResultSet = sqlStatement.executeQuery("SELECT scanNum, scanName, precursorCharge, precursorMass, peptide, theoMass, isDecoy, globalRank, normalizedCorrelationCoefficient, score, deltaLCn, deltaCn, matchedPeakNum, ionFrac, matchedHighestIntensityFrac, explainedAaFrac FROM spectraTable");

        while (sqlResultSet.next()) {
            String peptide = sqlResultSet.getString("peptide");
            if (!sqlResultSet.wasNull()) {
                int charge = sqlResultSet.getInt("precursorCharge");
                double theoMass = sqlResultSet.getDouble("theoMass");
                double expMass = sqlResultSet.getDouble("precursorMass");
                double massDiff = getMassDiff(expMass, theoMass, MassTool.C13_DIFF);

                Peptide0 peptide0 = peptide0Map.get(peptide.replaceAll("[^ncA-Z]+", ""));
                TreeSet<String> proteinIdSet = new TreeSet<>();
                for (String protein : peptide0.proteins) {
                    proteinIdSet.add(protein.trim());
                }

                StringBuilder sb = new StringBuilder(20);
                for (int i = 0; i < 6; ++i) {
                    if (i == charge - 1) {
                        sb.append(1);
                    } else {
                        sb.append(0);
                    }
                    sb.append("\t");
                }

                String scanName = sqlResultSet.getString("scanName");
                int scanNum = sqlResultSet.getInt("scanNum");
                int isDecoy = sqlResultSet.getInt("isDecoy");
                int globalRank = sqlResultSet.getInt("globalRank");
                double normalizedCorrelationCoefficient = sqlResultSet.getDouble("normalizedCorrelationCoefficient");
                double score = sqlResultSet.getDouble("score");
                double deltaLCn = sqlResultSet.getDouble("deltaLCn");
                double deltaCn = sqlResultSet.getDouble("deltaCn");
                double ionFrac = sqlResultSet.getDouble("ionFrac");
                double matchedHighestIntensityFrac = sqlResultSet.getDouble("matchedHighestIntensityFrac");
                double explainedAaFrac = sqlResultSet.getDouble("explainedAaFrac");

                String[] scanNameStr = scanName.split("\\.");
                int fileId = Integer.valueOf(scanNameStr[0]);
                boolean hasCo = fileIdScanNumToChargeScanNameMap.containsKey(fileId+"."+scanNum) && fileIdScanNumToChargeScanNameMap.get(fileId+"."+scanNum).containsKey(charge);
                if (!hasCo) {
                    if (isDecoy == 1) {
                        writer.write(scanName + "\t-1\t" + scanNum + "\t" + expMass + "\t" + theoMass + "\t" + score + "\t" + deltaCn + "\t" + deltaLCn + "\t" + normalizedCorrelationCoefficient + "\t" + globalRank + "\t" + Math.abs(massDiff * 1e6 / theoMass) + "\t" + ionFrac + "\t" + matchedHighestIntensityFrac + "\t" + sb.toString() + explainedAaFrac + "\t" + peptide0.leftFlank + "." + peptide.replaceAll("\\(", "[").replaceAll("\\)", "]") + "." + peptide0.rightFlank + "\t" + String.join("\t", proteinIdSet) + "\n"); // Percolator only recognize "[]".
                    } else {
                        writer.write(scanName + "\t1\t" + scanNum + "\t" + expMass + "\t" + theoMass + "\t" + score + "\t" + deltaCn + "\t" + deltaLCn + "\t" + normalizedCorrelationCoefficient + "\t" + globalRank + "\t" + Math.abs(massDiff * 1e6 / theoMass) + "\t" + ionFrac + "\t" + matchedHighestIntensityFrac + "\t" + sb.toString() + explainedAaFrac + "\t" + peptide0.leftFlank + "." + peptide.replaceAll("\\(", "[").replaceAll("\\)", "]") + "." + peptide0.rightFlank + "\t" + String.join("\t", proteinIdSet) + "\n"); // Percolator only recognize "[]".
                    }
                }
            }
        }
        sqlResultSet.close();
        sqlStatement.close();

        PreparedStatement sqlUpdateCo = sqlConnection.prepareStatement("REPLACE INTO spectraTable (scanNum, scanName,  precursorCharge, precursorMass, mgfTitle" +
                ", isotopeCorrectionNum, ms1PearsonCorrelationCoefficient, labelling, peptide, theoMass" +
                ", isDecoy, globalRank, normalizedCorrelationCoefficient, score, deltaLCn" +
                ", deltaCn, matchedPeakNum, ionFrac, matchedHighestIntensityFrac, explainedAaFrac" +
                ", otherPtmPatterns, aScore, candidates, peptideSet, whereIsTopCand" +
                ", shouldPtm, hasPTM, ptmNum, isSettled) VALUES (?, ?, ?, ?,?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)");
        sqlConnection.setAutoCommit(false);
        for (String fileIdScanNum : fileIdScanNumToChargeScanNameMap.keySet()) {
            if (fileIdScanNum.contentEquals("0.1972")) {
                int a = 1;
            }
            Map<Integer, Set<String>> chargeScanNameMap = fileIdScanNumToChargeScanNameMap.get(fileIdScanNum);
            for (int charge : chargeScanNameMap.keySet()) {
                Map<String, Double> seqScoreMap = new HashMap<>();
                Set<String> scanNameSet = chargeScanNameMap.get(charge);
                Iterator<String> scanNameIter = scanNameSet.iterator();
                while (scanNameIter.hasNext()) {
                    if (!exResForScanNameMap.containsKey(scanNameIter.next())) {
                        scanNameIter.remove();
                    }
                }

                Set<String> selectedPep = new HashSet<>();
                while (0 < scanNameSet.size()) {
                    String bestPep = "";
                    String bestScanName = "";
                    double bestScore = 0d;

                    for (String scanName : scanNameSet) {
                        for (ExRes res : exResForScanNameMap.get(scanName)) {
                            if (res.score > bestScore && !selectedPep.contains(res.peptide)) {
                                bestScore = res.score;
                                bestPep = res.peptide;
                                bestScanName = scanName;
                            }
                        }
                    }

                    if (exResForScanNameMap.get(bestScanName) == null) {
                        int a = 1;
                    }
                    if (bestScanName.contentEquals("")) {
                        break;
                    }
                    ExRes bestRes = exResForScanNameMap.get(bestScanName).get(0);
                    for (ExRes res : exResForScanNameMap.get(bestScanName)) {
                        if (res.peptide.contentEquals(bestPep)) {
                            bestRes = res;
                        }
                    }
                    StringBuilder sb = new StringBuilder(20);
                    for (int i = 0; i < 6; ++i) {
                        if (i == bestRes.precursorCharge - 1) {
                            sb.append(1);
                        } else {
                            sb.append(0);
                        }
                        sb.append("\t");
                    }
                    double massDiff = getMassDiff(bestRes.precursorMass, bestRes.theoMass, MassTool.C13_DIFF);
                    Peptide0 peptide0 = peptide0Map.get(bestRes.peptide.replaceAll("[^ncA-Z]+", ""));
                    TreeSet<String> proteinIdSet = new TreeSet<>();
                    for (String protein : peptide0.proteins) {
                        proteinIdSet.add(protein.trim());
                    }
                    if (bestRes.isDecoy == 1) {
                        writer.write(bestScanName + "\t-1\t" + bestRes.scanNum + "\t" + bestRes.precursorMass + "\t" + bestRes.theoMass + "\t" + bestRes.score + "\t" + bestRes.deltaCn + "\t" + bestRes.deltaLCn + "\t" + bestRes.normalizedCorrelationCoefficient + "\t" + bestRes.globalRank + "\t" + Math.abs(massDiff * 1e6 / bestRes.theoMass) + "\t" + bestRes.ionFrac + "\t" + bestRes.matchedHighestIntensityFrac + "\t" + sb.toString() + bestRes.explainedAaFrac + "\t" + peptide0.leftFlank + "." + bestRes.peptide.replaceAll("\\(", "[").replaceAll("\\)", "]") + "." + peptide0.rightFlank + "\t" + String.join("\t", proteinIdSet) + "\n"); // Percolator only recognize "[]".
                    } else {
                        writer.write(bestScanName + "\t1\t" + bestRes.scanNum + "\t" + bestRes.precursorMass + "\t" + bestRes.theoMass + "\t" + bestRes.score + "\t" + bestRes.deltaCn + "\t" + bestRes.deltaLCn + "\t" + bestRes.normalizedCorrelationCoefficient + "\t" + bestRes.globalRank + "\t" + Math.abs(massDiff * 1e6 / bestRes.theoMass) + "\t" + bestRes.ionFrac + "\t" + bestRes.matchedHighestIntensityFrac + "\t" + sb.toString() + bestRes.explainedAaFrac + "\t" + peptide0.leftFlank + "." + bestRes.peptide.replaceAll("\\(", "[").replaceAll("\\)", "]") + "." + peptide0.rightFlank + "\t" + String.join("\t", proteinIdSet) + "\n"); // Percolator only recognize "[]".
                    }
                    scanNameSet.remove(bestScanName);
                    selectedPep.add(bestPep);
                    sqlUpdateCo.setInt(1, bestRes.scanNum);
                    sqlUpdateCo.setString(2, bestRes.scanName);
                    sqlUpdateCo.setInt(3, bestRes.precursorCharge);
                    sqlUpdateCo.setDouble(4, bestRes.precursorMass);
                    sqlUpdateCo.setString(5, "fakeMgfTitle");
                    sqlUpdateCo.setInt(6, -99);
                    sqlUpdateCo.setDouble(7, -99.9);
                    sqlUpdateCo.setString(8, bestRes.labelling);
                    sqlUpdateCo.setString(9, bestRes.peptide);
                    sqlUpdateCo.setDouble(10, bestRes.theoMass);
                    sqlUpdateCo.setInt(11, bestRes.isDecoy);
                    sqlUpdateCo.setInt(12, bestRes.globalRank);
                    sqlUpdateCo.setDouble(13, bestRes.normalizedCorrelationCoefficient);
                    sqlUpdateCo.setDouble(14, bestRes.score);
                    sqlUpdateCo.setDouble(15, bestRes.deltaLCn);
                    sqlUpdateCo.setDouble(16, bestRes.deltaCn);
                    sqlUpdateCo.setInt(17, bestRes.matchedPeakNum);
                    sqlUpdateCo.setDouble(18, bestRes.ionFrac);
                    sqlUpdateCo.setDouble(19, bestRes.matchedHighestIntensityFrac);
                    sqlUpdateCo.setDouble(20, bestRes.explainedAaFrac);
                    sqlUpdateCo.setString(21, bestRes.otherPtmPatterns);
                    sqlUpdateCo.setString(22, bestRes.aScore);
                    sqlUpdateCo.setString(23, bestRes.candidates);
                    sqlUpdateCo.setString(24, bestRes.peptideSet);
                    sqlUpdateCo.setInt(25, bestRes.whereIsTopCand);
                    sqlUpdateCo.setInt(26, bestRes.shouldPtm);
                    sqlUpdateCo.setInt(27, bestRes.hasPTM);
                    sqlUpdateCo.setInt(28, bestRes.ptmNum);
                    sqlUpdateCo.setInt(29, bestRes.isSettled);

                    sqlUpdateCo.executeUpdate();
                }
            }
        }
        sqlConnection.commit();
        sqlConnection.setAutoCommit(true);
        writer.close();
        sqlConnection.close();
    }

    private static Map<String, PercolatorEntry> runPercolator(String percolatorPath, String percolatorInputFileName, String percolatorOutputFileName, String percolatorProteinOutputFileName, String tdFastaPath, String enzymeName) throws Exception {
        Map<String, PercolatorEntry> percolatorResultMap = new HashMap<>();
        if ((new File(percolatorInputFileName)).exists()) {
            String percolatorEnzyme;
            switch (enzymeName) {
                case "Trypsin": percolatorEnzyme = "trypsin";
                    break;
                case "Trypsin/P": percolatorEnzyme = "trypsinp";
                    break;
                case "TrypsinR": percolatorEnzyme = "trypsin";
                    break;
                case "LysC": percolatorEnzyme = "lys-c";
                    break;
                case "ArgC": percolatorEnzyme = "arg-c";
                    break;
                case "Chymotrypsin": percolatorEnzyme = "chymotrypsin";
                    break;
                case "GluC": percolatorEnzyme = "glu-c";
                    break;
                case "LysN": percolatorEnzyme = "lys-n";
                    break;
                case "AspN": percolatorEnzyme = "asp-n";
                    break;
                default: percolatorEnzyme = "trypsin";
                    break;
            }

            String[] commands = new String[]{percolatorPath, "--only-psms", "--verbose", "0", "--no-terminate", "--search-input", "concatenated", "--protein-decoy-pattern", "DECOY_", "--picked-protein", tdFastaPath, "--protein-enzyme", percolatorEnzyme, "--protein-report-fragments", "--protein-report-duplicates", "--results-proteins", percolatorProteinOutputFileName, "--results-psms", percolatorOutputFileName, percolatorInputFileName};

            Process ps = Runtime.getRuntime().exec(commands);
            ps.waitFor();

            if (!(new File(percolatorOutputFileName).exists()) || ps.exitValue() != 0) {
                BufferedReader reader = new BufferedReader(new InputStreamReader(ps.getInputStream()));
                String line;
                while ((line = reader.readLine()) != null) {
                    logger.info("[Percolator info]: {}", line.trim());
                }
                reader.close();
                reader = new BufferedReader(new InputStreamReader(ps.getErrorStream()));
                while ((line = reader.readLine()) != null) {
                    logger.info("[Percolator error]: {}", line.trim());
                }
                reader.close();
                throw new NullPointerException("Percolator didn't exit normally.");
            } else {
                BufferedReader reader = new BufferedReader(new InputStreamReader(ps.getInputStream()));
                String line;
                while ((line = reader.readLine()) != null) {
                    logger.info("[Percolator info]: {}", line.trim());
                }
                reader.close();
                reader = new BufferedReader(new InputStreamReader(ps.getErrorStream()));
                while ((line = reader.readLine()) != null) {
                    logger.info("[Percolator info]: {}", line.trim());
                }
                reader.close();
            }

            logger.info("Percolator finished normally");

            BufferedReader reader = new BufferedReader(new FileReader(percolatorOutputFileName));
            String line;
            while ((line = reader.readLine()) != null) {
                line = line.trim();
                if (!line.startsWith("PSMId")) {
                    String[] parts = line.split("\t");
                    percolatorResultMap.put(parts[0], new PercolatorEntry(Double.valueOf(parts[1]), parts[2], parts[3]));
                }
            }
            reader.close();
        } else {
            logger.error("Cannot find Percolator input file (from {}) for estimating Percolator Q-Value.", percolatorInputFileName);
            return percolatorResultMap;
        }

        return percolatorResultMap;
    }

    private void writeTop20Candidates(String sqlPath, String spectraPath) throws IOException, SQLException {
        System.out.println("PFM starts");
        Connection sqlConnection = DriverManager.getConnection(sqlPath);
        Statement sqlStatement = sqlConnection.createStatement();
        ResultSet sqlResultSet = sqlStatement.executeQuery("SELECT scanNum, candidates, whereIsTopCand FROM spectraTable");

        String resultPath = spectraPath + "Chick9268.csv";
//        String resultPath = spectraPath + "So70S.csv";
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(resultPath))) {
            writer.write("scanNo,whereIsTopCand,pep1,s1,b1,m1,p1,pep2,s2,b2,m2,p2,pep3,s3,b3,m3,p3,pep4,s4,b4,m4,p4,pep5,s5,b5,m5,p5,pep6,s6,b6,m6,p6,pep7,s7,b7,m7,p7,pep8,s8,b8,m8,p8,pep9,s9,b9,m9,p9,pep10,s10,b10,m10,p10,pep11,s11,b11,m11,p11,pep12,s12,b12,m12,p12,pep13,s13,b13,m13,p13,pep14,s14,b14,m14,p14,pep15,s15,b15,m15,p15,pep16,s16,b16,m16,p16,pep17,s17,b17,m17,p17,pep18,s18,b18,m18,p18,pep19,s19,b19,m19,p19,pep20,s20,b20,m20,p20\n");
            while (sqlResultSet.next()) {
                int scanNum = sqlResultSet.getInt("scanNum");
                int whereIsTopCand = sqlResultSet.getInt("whereIsTopCand");
                String candidatesList = sqlResultSet.getString("candidates");
                if (!sqlResultSet.wasNull()) {
                    writer.write(scanNum + "," + whereIsTopCand + "," + candidatesList + "\n");
                }
            }
        } catch (IOException ex) {
            ex.printStackTrace();
            logger.error(ex.getMessage());
            System.exit(1);
        }
        sqlResultSet.close();
        sqlStatement.close();
        sqlConnection.close();

        Connection sqlConnection2 = DriverManager.getConnection(sqlPath);
        Statement sqlStatement2 = sqlConnection2.createStatement();
        ResultSet sqlResultSet2 = sqlStatement2.executeQuery("SELECT scanNum, peptide, isDecoy, score FROM spectraTable");
        String noPercoPath = spectraPath + "Chick9268TopNoPerco.csv";
//        String noPercoPath = spectraPath + "So70STopNoPerco.csv";
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(noPercoPath))) {
            writer.write("scanNo,pep1,s1,b1,m1,p1\n");
            while (sqlResultSet2.next()) {
                int scanNum = sqlResultSet2.getInt("scanNum");
                String peptide = sqlResultSet2.getString("peptide");
                double score = sqlResultSet2.getDouble("score");
                int isDecoy = sqlResultSet2.getInt("isDecoy");
                if (!sqlResultSet2.wasNull()) {
                    writer.write(scanNum + "," + peptide + "," + score +"," + isDecoy +"\n");
                }
            }
        } catch (IOException ex) {
            ex.printStackTrace();
            logger.error(ex.getMessage());
            System.exit(1);
        }

        sqlResultSet2.close();
        sqlStatement2.close();
        sqlConnection2.close();

        Connection sqlConnection3 = DriverManager.getConnection(sqlPath);
        Statement sqlStatement3 = sqlConnection3.createStatement();
        ResultSet sqlResultSet3 = sqlStatement3.executeQuery("SELECT scanNum, peptideSet FROM spectraTable");
        String candisWithPTMXcorrPath = spectraPath + "Chick9268Candis20WithPTMXcorr.csv";
//        String candisWithPTMXcorrPath = spectraPath + "So70SCandis20WithPTMXcorr.csv";
        try (BufferedWriter writer = new BufferedWriter(new FileWriter(candisWithPTMXcorrPath))) {
            writer.write("scanNo,pep1,s1,b1,m1,p1,pep2,s2,b2,m2,p2,pep3,s3,b3,m3,p3,pep4,s4,b4,m4,p4,pep5,s5,b5,m5,p5,pep6,s6,b6,m6,p6,pep7,s7,b7,m7,p7,pep8,s8,b8,m8,p8,pep9,s9,b9,m9,p9,pep10,s10,b10,m10,p10,pep11,s11,b11,m11,p11,pep12,s12,b12,m12,p12,pep13,s13,b13,m13,p13,pep14,s14,b14,m14,p14,pep15,s15,b15,m15,p15,pep16,s16,b16,m16,p16,pep17,s17,b17,m17,p17,pep18,s18,b18,m18,p18,pep19,s19,b19,m19,p19,pep20,s20,b20,m20,p20\n");
            while (sqlResultSet3.next()) {
                int scanNum = sqlResultSet3.getInt("scanNum");
                String peptideSet = sqlResultSet3.getString("peptideSet");
                if (!sqlResultSet3.wasNull()) {
                    writer.write(scanNum + "," + peptideSet +"\n");
                }
            }
        } catch (IOException ex) {
            ex.printStackTrace();
            logger.error(ex.getMessage());
            System.exit(1);
        }

        sqlResultSet3.close();
        sqlStatement3.close();
        sqlConnection3.close();
    }
    private void writeFinalResult(Map<String, PercolatorEntry> percolatorResultMap, String outputPath, Map<String, Peptide0> peptide0Map, String sqlPath) throws IOException, SQLException {
        TreeMap<Double, List<String>> tempMap = new TreeMap<>();

        BufferedWriter writer = new BufferedWriter(new FileWriter(outputPath));
        if (percolatorResultMap == null) {
            writer.write("scanName,scan_num,peptide,isDecoy, shouldPtm, hasPTM, ptmNum, isSettled,charge,theo_mass,exp_mass,abs_ppm,A_score,protein_ID,score,delta_C_n,deltaLCn,globalRank, " +
                    "cosScore,ionFrac,matchedHighestIntensityFrac,explainedAaFrac, other_PTM_patterns,MGF_title,labelling,isotope_correction,MS1_pearson_correlation_coefficient\n");
        } else {
            writer.write("scanName,scan_num,peptide, isDecoy, shouldPtm, hasPTM, ptmNum, isSettled,charge,theo_mass,exp_mass,abs_ppm,A_score,protein_ID,score,delta_C_n,deltaLCn,globalRank, " +
                    "cosScore,ionFrac,matchedHighestIntensityFrac,explainedAaFrac,percolator_score,posterior_error_prob,q_value,other_PTM_patterns,MGF_title,labelling,isotope_correction,MS1_pearson_correlation_coefficient\n");
        }

        Connection sqlConnection = DriverManager.getConnection(sqlPath);
        Statement sqlStatement = sqlConnection.createStatement();
        ResultSet sqlResultSet = sqlStatement.executeQuery("SELECT scanName,scanNum, shouldPtm, hasPTM, ptmNum, isSettled, precursorCharge, precursorMass, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient, " +
                "labelling, peptide, theoMass, isDecoy, score, otherPtmPatterns, aScore, deltaCn,deltaLCn, isDecoy, globalRank, normalizedCorrelationCoefficient, matchedPeakNum, ionFrac, matchedHighestIntensityFrac, explainedAaFrac FROM spectraTable");
        int numScansQ001 = 0;

        while (sqlResultSet.next()) {
            int isDecoy = sqlResultSet.getInt("isDecoy");
            if (!sqlResultSet.wasNull()) {
                if (isDecoy == 0) {
                    String scanName = sqlResultSet.getString("scanName");
                    int scanNum = sqlResultSet.getInt("scanNum");
                    double expMass = sqlResultSet.getDouble("precursorMass");
                    String peptide = sqlResultSet.getString("peptide");
                    double theoMass = sqlResultSet.getDouble("theoMass");
                    double massDiff = getMassDiff(expMass, theoMass, MassTool.C13_DIFF);
                    double ppm = Math.abs(massDiff * 1e6 / theoMass);
                    int globalRank = sqlResultSet.getInt("globalRank");
                    double cosScore = sqlResultSet.getDouble("normalizedCorrelationCoefficient");
                    double score = sqlResultSet.getDouble("score");
                    double deltaLCn = sqlResultSet.getDouble("deltaLCn");
                    double ionFrac = sqlResultSet.getDouble("ionFrac");
                    double matchedHighestIntensityFrac = sqlResultSet.getDouble("matchedHighestIntensityFrac");
                    double explainedAaFrac = sqlResultSet.getDouble("explainedAaFrac");

                    Peptide0 peptide0 = peptide0Map.get(peptide.replaceAll("[^ncA-Z]+", ""));
                    TreeSet<String> proteinIdSet = new TreeSet<>();
                    for (String protein : peptide0.proteins) {
                        proteinIdSet.add(protein.trim());
                    }

                    String aScore = sqlResultSet.getString("aScore");

//                    if (percolatorResultMap == null) {
//                        writer.write("scan_num,peptide,isDecoy, shouldPtm, hasPTM, ptmNum, isSettled,charge,theo_mass,exp_mass,abs_ppm,A_score,protein_ID,score,delta_C_n,deltaLCn,globalRank, " +
//                                "cosScore,matchedPeakNum,ionFrac,matchedHighestIntensityFrac,explainedAaFrac, other_PTM_patterns,MGF_title,labelling,isotope_correction,MS1_pearson_correlation_coefficient\n");
//                    } else {
//                        writer.write("scan_num,peptide, isDecoy, shouldPtm, hasPTM, ptmNum, isSettled,charge,theo_mass,exp_mass,abs_ppm,A_score,protein_ID,score,delta_C_n,deltaLCn,globalRank, " +
//                                "cosScore, matchedPeakNum,ionFrac,matchedHighestIntensityFrac,explainedAaFrac,percolator_score,posterior_error_prob,q_value,other_PTM_patterns,MGF_title,labelling,isotope_correction,MS1_pearson_correlation_coefficient\n");
//                    }
                    if (percolatorResultMap == null) {
//                        double score = sqlResultSet.getDouble("score");
                        String str = String.format(Locale.US, "%s,%d,%s,%d,%d,%d,%d,%d,%d,%f,%f,%f,%s,%s,%f,%f,%f,%d,%f,%f,%f,%f,%s,\"%s\",%s,%d,%f\n"
                                , scanName,scanNum, peptide,isDecoy, sqlResultSet.getInt("shouldPtm"),sqlResultSet.getInt("hasPTM"),sqlResultSet.getInt("ptmNum"),sqlResultSet.getInt("isSettled")
                                , sqlResultSet.getInt("precursorCharge"), theoMass, expMass, ppm, aScore, String.join(";", proteinIdSet).replaceAll(",", "~"), score
                                , sqlResultSet.getDouble("deltaCn"),deltaLCn,globalRank,cosScore, ionFrac, matchedHighestIntensityFrac,explainedAaFrac
                                , sqlResultSet.getString("otherPtmPatterns"), sqlResultSet.getString("mgfTitle")
                                , sqlResultSet.getString("labelling"), sqlResultSet.getInt("isotopeCorrectionNum"), sqlResultSet.getDouble("ms1PearsonCorrelationCoefficient"));

                        if (tempMap.containsKey(score)) {
                            tempMap.get(score).add(str);
                        } else {
                            List<String> tempList = new LinkedList<>();
                            tempList.add(str);
                            tempMap.put(score, tempList);
                        }
                    } else {
                        if (!percolatorResultMap.containsKey(scanName)) {
                            continue; // this is temp way to avoid that the isDecoy of a scanName can change and be different from sql because of losing competition against other co-spectra
                        }
                        PercolatorEntry percolatorEntry = percolatorResultMap.get(scanName);

                        if (Double.valueOf(percolatorEntry.qValue) <= 0.01) numScansQ001++;
                        String str = String.format(Locale.US, "%s, %d,%s,%d,%d,%d,%d,%d,%d,%f,%f,%f,%s,%s,%f,%f,%f,%d,%f,%f,%f,%f,%f,%s,%s,%s,\"%s\",%s,%d,%f\n"
                                , scanName, scanNum, peptide,isDecoy, sqlResultSet.getInt("shouldPtm"),sqlResultSet.getInt("hasPTM"),sqlResultSet.getInt("ptmNum"),sqlResultSet.getInt("isSettled")
                                , sqlResultSet.getInt("precursorCharge"), theoMass, expMass, ppm, aScore, String.join(";", proteinIdSet).replaceAll(",", "~"), score
                                , sqlResultSet.getDouble("deltaCn"),deltaLCn,globalRank,cosScore, ionFrac, matchedHighestIntensityFrac,explainedAaFrac, percolatorEntry.percolatorScore, percolatorEntry.PEP, percolatorEntry.qValue
                                , sqlResultSet.getString("otherPtmPatterns"), sqlResultSet.getString("mgfTitle")
                                , sqlResultSet.getString("labelling"), sqlResultSet.getInt("isotopeCorrectionNum"), sqlResultSet.getDouble("ms1PearsonCorrelationCoefficient"));

                        if (tempMap.containsKey(percolatorResultMap.get(scanName).percolatorScore)) {
                            tempMap.get(percolatorResultMap.get(scanName).percolatorScore).add(str);
                        } else {
                            List<String> tempList = new LinkedList<>();
                            tempList.add(str);
                            tempMap.put(percolatorResultMap.get(scanName).percolatorScore, tempList);
                        }
                    }
                }
            }
        }

        sqlResultSet.close();
        sqlStatement.close();
        sqlConnection.close();

        Double[] tempArray = tempMap.keySet().toArray(new Double[0]); //this is for writing the final results in descending order of percolator score
        for (int i = tempArray.length - 1; i >= 0; --i) {
            List<String> tempList = tempMap.get(tempArray[i]);
            for (String tempStr : tempList) {
                writer.write(tempStr);
            }
        }
        writer.close();
        logger.info("Number of PSMs with q < 0.01: {}", numScansQ001);

// ======== down is for comparison ======
        BufferedReader parameterReader = new BufferedReader(new FileReader("/home/slaiad/Code/PIPI/src/main/resources/modPlusRes.txt"));
        Map<Integer, String> pepTruth = new HashMap<>();
        Map<Integer, String> backboneTruth = new HashMap<>();
        String line;
        while ((line = parameterReader.readLine()) != null) {
            line = line.trim();
            String[] splitRes = line.split(",");
            String pepWithMod = splitRes[1].substring(2, splitRes[1].length()-2).replace('L','I');
            String backbone = pepWithMod.replaceAll("[^a-zA-Z]", "");
            pepTruth.put(Integer.valueOf(splitRes[0]), pepWithMod);
            backboneTruth.put(Integer.valueOf(splitRes[0]), backbone);
        }

        Set<Integer> numNoSuchScan = new HashSet<>(pepTruth.keySet());
        Set<Integer>  numNotContainsTruth = new HashSet<>();
        Set<Integer>  numContains_WrongPep = new HashSet<>();
        Set<Integer>  numContains_CorrectPep_badQ = new HashSet<>();
        Set<Integer>  numContains_CorrectPep_goodQ_wrongPtm = new HashSet<>();
        Set<Integer>  numContains_CorrectPep_goodQ_goodPtm = new HashSet<>();
//        Set<Integer> scanInPipiAndModp = new HashSet<>();
        Connection sqlConnection3 = DriverManager.getConnection(sqlPath);
        Statement sqlStatement3 = sqlConnection3.createStatement();
        ResultSet sqlResultSet3 = sqlStatement3.executeQuery("SELECT scanNum, peptideSet, peptide , candidates FROM spectraTable");
        while (sqlResultSet3.next()) {
            int scanNum = sqlResultSet3.getInt("scanNum");
            if (!pepTruth.containsKey(scanNum)) continue;
            numNoSuchScan.remove(scanNum);
            String modpPep = pepTruth.get(scanNum).replace('L','I');
            String modpBackbone = backboneTruth.get(scanNum).replace('L','I');
            String peptide = sqlResultSet3.getString("peptide");
//            String peptideSetString = sqlResultSet3.getString("peptideSet");
//            if (scanNum == 2050 ){
//                System.out.println(peptideSetString);
//            }
            String peptideSetString = sqlResultSet3.getString("candidates");
            if (!sqlResultSet3.wasNull()) {
                String[] tempPeptideSet = peptideSetString.split(",");
                Set<String> peptideSet = new HashSet<>();
                for (String str : tempPeptideSet) {
                    if (str.startsWith("n") && str.endsWith("c")){
                        peptideSet.add(str.substring(1, str.length()-1));
                    }
                }
                if (!peptideSet.contains(modpBackbone)) {
                    numNotContainsTruth.add(scanNum);
                    continue;
                }

                String myPep = peptide.substring(1, peptide.length()-1).replace('L','I');
                String myBackbone = myPep.replaceAll("\\(-", "-");
                myBackbone = myBackbone.replaceAll("[\\(]", "+");
                myBackbone = myBackbone.replaceAll("[\\)]", "");

                if (!myBackbone.contentEquals(modpBackbone)) {
                    numContains_WrongPep.add(scanNum);
                    continue;
                }
                if (Double.valueOf(percolatorResultMap.get(scanNum).qValue) > 0.01 ){
                    numContains_CorrectPep_badQ.add(scanNum);
                    continue;
                }
                if (!myPep.contentEquals(modpPep)) {
                    numContains_CorrectPep_goodQ_wrongPtm.add(scanNum);
                    continue;
                }
                numContains_CorrectPep_goodQ_goodPtm.add(scanNum);
            }
        }
//        System.out.println("numNoSuchScan, "+numNoSuchScan.size());
//        System.out.println("numNotContainsTruth, "+numNotContainsTruth.size());
//        System.out.println("numContains_WrongPep, "+numContains_WrongPep.size());
//        System.out.println("numContains_CorrectPep_badQ, "+numContains_CorrectPep_badQ.size());
//        System.out.println("numContains_CorrectPep_goodQ_wrongPtm, "+numContains_CorrectPep_goodQ_wrongPtm.size());
//        System.out.println("numContains_CorrectPep_goodQ_goodPtm, "+numContains_CorrectPep_goodQ_goodPtm.size());
//        System.out.println("=======");
//        System.out.println("numNoSuchScan, "+numNoSuchScan);
////        System.out.println("numNotContainsTruth, "+numNotContainsTruth);
//        System.out.println("numContains_WrongPep, "+numContains_WrongPep);
//        System.out.println("numContains_CorrectPep_badQ, "+numContains_CorrectPep_badQ);
//        System.out.println("numContains_CorrectPep_goodQ_wrongPtm, "+numContains_CorrectPep_goodQ_wrongPtm);
//        System.out.println("numContains_CorrectPep_goodQ_goodPtm, "+numContains_CorrectPep_goodQ_goodPtm);


        sqlResultSet3.close();
        sqlStatement3.close();
        sqlConnection3.close();

    }

    public static double getMassDiff(double expMass, double theoMass, double C13Diff) {
        double massDiff1 = expMass - theoMass;
        double massDiff2 = expMass - theoMass - C13Diff;
        double massDiff3 = expMass - theoMass - 2 * C13Diff;
        double absMassDiff1 = Math.abs(massDiff1);
        double absMassDiff2 = Math.abs(massDiff2);
        double absMassDiff3 = Math.abs(massDiff3);

        if ((absMassDiff1 <= absMassDiff2) && (absMassDiff1 <= absMassDiff2)) {
            return massDiff1;
        } else if ((absMassDiff2 <= absMassDiff1) && (absMassDiff2 <= absMassDiff3)) {
            return massDiff2;
        } else {
            return massDiff3;
        }
    }
}
