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

import ProteomicsLibrary.*;
import org.apache.commons.math3.util.Pair;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;
import proteomics.FM.FMIndex;
import proteomics.PTM.InferPTM;
import proteomics.Types.*;
import proteomics.Index.BuildIndex;
import proteomics.Parameter.Parameter;
import proteomics.Spectrum.DatasetReader;
import uk.ac.ebi.pride.tools.jmzreader.JMzReader;
import uk.ac.ebi.pride.tools.mgf_parser.MgfFile;
import uk.ac.ebi.pride.tools.mzxml_parser.MzXMLFile;


import java.io.*;
import java.net.InetAddress;
import java.net.UnknownHostException;
import java.sql.*;
import java.text.DecimalFormat;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.locks.ReentrantLock;

public class PIPI {
    private static final Logger logger = LoggerFactory.getLogger(PIPI.class);
    public static final String versionStr = "3";
    static final boolean useXcorr = false;
    public static final int minTagLenToReduceProtDb = 5;
    public static final boolean isDebugMode = java.lang.management.ManagementFactory.getRuntimeMXBean().getInputArguments().toString().indexOf("jdwp") >= 0;

    //// normal
//    static final int minTagLenToExtract = 4;  //normal //todo
//    static final boolean nTermSpecific = false; //normal //todo
//    public    static final double MIN_PEAK_SUM_INFER_AA = 0.4;
//    static final double proteinCovThres = 0.02;//0.02 is good for normal and DL dataset.0.1 is good for synthetic
//    static final int  maxNumVarPtmConsidered = 5;

    // dimethyl
//    static final int minTagLenToExtract = 4;  //normal //todo
//    static final boolean nTermSpecific = false; //normal //todo
//    public    static final double MIN_PEAK_SUM_INFER_AA = 0.4;
//    static final double proteinCovThres = 0.02;//0.02 is good for normal and DL dataset.0.1 is good for synthetic
//    static final int  maxNumVarPtmConsidered = 2;

    ///synthetic
//    static final int minTagLenToExtract = 3;  //normal //todo
//    static final boolean nTermSpecific = true; //normal //todo
//    public static final double MIN_PEAK_SUM_INFER_AA = 0.0;
//    static final double proteinCovThres = 0.1;//0.02 is good for normal and DL dataset.0.1 is good for synthetic
//    static final int  maxNumVarPtmConsidered = 1;

    //    / DL simu
    static final int minTagLenToExtract = 3;  //normal //todo
    public static final double MIN_PEAK_SUM_INFER_AA = 0.0;
    static final int  maxNumVarPtmConsidered = 18;
//    /debuging parameters
    public static HashSet<Integer> lszDebugScanNum = new HashSet<>(Arrays.asList(82001));//178,179,180,181,183,184,192

    public static void main(String[] args) {
        long startTime = System.nanoTime();
        if (args.length != 1) {
            help();
        }
        // Set parameters
        String parameterPath = args[0].trim();
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
            dbName = String.format(Locale.US, "PIPI.temp.db");
            new PIPI(parameterPath, dbName, hostName);
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

    private PIPI(String parameterPath, String dbName, String hostName) throws Exception {
//        if (isDebugMode) {
            ch.qos.logback.classic.Logger root = (ch.qos.logback.classic.Logger) org.slf4j.LoggerFactory.getLogger(ch.qos.logback.classic.Logger.ROOT_LOGGER_NAME);
            root.setLevel(ch.qos.logback.classic.Level.DEBUG);
//        }
        // Get the parameter map
        Parameter parameter = new Parameter(parameterPath);
        Map<String, String> parameterMap = parameter.returnParameterMap();
        double ms2Tol = Double.valueOf(parameterMap.get("ms2_tolerance"));
        double ms1Tol = Double.valueOf(parameterMap.get("ms1_tolerance"));
//        double leftInverseMs1Tolerance = 1 / (1 + ms1Tol * 1e-6);
//        double rightInverseMs1Tolerance = 1 / (1 - ms1Tol * 1e-6);
        boolean addContaminants = Boolean.valueOf(parameterMap.get("add_contaminant"));
        double proteinCovThres = Double.valueOf(parameterMap.get("min_prot_coverage"));
        String spectraPath = parameterMap.get("spectra_path");
        String outputDir = parameterMap.get("output_dir");



        logger.info("Spectra: {}, parameter: {}.", parameterPath);

        Set<Integer> msLevelSet = new HashSet<>();
        msLevelSet.add(2);

        logger.info("Loading parameters and build fmIndex...");
        BuildIndex buildIndex = new BuildIndex(parameterMap);
        MassTool massTool = buildIndex.returnMassTool();
        InferPTM inferPTM = buildIndex.getInferPTM();

        logger.info("Reading spectra...");
        File spectraFile = new File(spectraPath);
        DatasetReader datasetReader;
        JMzReader[] spectraParserArray;
        String sqlPath = "jdbc:sqlite:" + dbName;
        Class.forName("org.sqlite.JDBC").newInstance();
        Map<Integer, String> fileIdNameMap = new HashMap<>();
        Map<String, Integer> fileNameIdMap = new HashMap<>();
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
            fileIdNameMap.put(0, spectraPath.substring(spectraPath.lastIndexOf("/")+1).split("\\.")[0].replaceAll("\\.","_"));
            fileNameIdMap.put(spectraPath.substring(spectraPath.lastIndexOf("/")+1).split("\\.")[0].replaceAll("\\.","_"), 0);
            datasetReader = new DatasetReader(spectraParserArray, ms1Tol, massTool, ext, msLevelSet, sqlPath, fileIdNameMap);
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
                fileIdNameMap.put(i, fileList[i].split("\\.")[0].replaceAll("\\.","_"));
                fileNameIdMap.put(fileList[i].split("\\.")[0].replaceAll("\\.","_"), i);
            }
            String ext = fileList[0].substring(fileList[0].lastIndexOf(".")+1);
            datasetReader = new DatasetReader(spectraParserArray, ms1Tol, massTool, ext, msLevelSet, sqlPath, fileIdNameMap);
        }

        SpecProcessor specProcessor = new SpecProcessor(massTool);
        Map<String, Integer> precursorChargeMap = new HashMap<>();
        Map<String, Double> precursorMassMap = new HashMap<>();

        logger.info("Preprocessing protein database...");
        int p_thread_num = Runtime.getRuntime().availableProcessors(); // for non main search use all threads

//        if (isDebugMode) p_thread_num = 1;
        logger.debug("thread NUM "+ p_thread_num);

        Map<String, String> mgfTitleMap = new HashMap<>();
        Map<String, Integer> isotopeCorrectionNumMap = new HashMap<>();
        Map<String, Integer> precursorScanNoMap = new HashMap<>();
        Map<String, Double> ms1PearsonCorrelationCoefficientMap = new HashMap<>();

        //////////==================================================
        //   Get Long Tags
        ExecutorService threadPoolGetLongTag = Executors.newFixedThreadPool(p_thread_num);
        ArrayList<Future<GetLongTag.Entry>> taskListGetLongTag = new ArrayList<>(datasetReader.getUsefulSpectraNum() + 10);
        Connection sqlConSpecCoderX = DriverManager.getConnection(sqlPath);
        Statement sqlStateGetLongTag = sqlConSpecCoderX.createStatement();
        ResultSet sqlResSetGetLongTag = sqlStateGetLongTag.executeQuery("SELECT scanName, scanNum, precursorCharge" +
                ", precursorMass, precursorScanNo, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient FROM spectraTable");

        ReentrantLock lockGetLongTag = new ReentrantLock();
        int submitNum = 0;
        Set<String> validScanSet = new HashSet<>();
        while (sqlResSetGetLongTag.next()) {
            String scanName = sqlResSetGetLongTag.getString("scanName");
            int scanNum = sqlResSetGetLongTag.getInt("scanNum");
            int precursorCharge = sqlResSetGetLongTag.getInt("precursorCharge");
            double precursorMass = sqlResSetGetLongTag.getDouble("precursorMass");
            mgfTitleMap.put(scanName, sqlResSetGetLongTag.getString("mgfTitle"));
            isotopeCorrectionNumMap.put(scanName, sqlResSetGetLongTag.getInt("isotopeCorrectionNum"));
            precursorScanNoMap.put(scanName, sqlResSetGetLongTag.getInt("precursorScanNo"));
            ms1PearsonCorrelationCoefficientMap.put(scanName, sqlResSetGetLongTag.getDouble("ms1PearsonCorrelationCoefficient"));
            if (isDebugMode){
                boolean shouldRun = false;
                for (int debugScanNum : lszDebugScanNum) {
                    if (Math.abs(scanNum-debugScanNum) < 2) {
                        shouldRun = true;
                    }
                }
                if (!shouldRun) continue;// comment this line, if dont skip
            }
            int fileId = fileNameIdMap.get( scanName.split("\\.")[0] );
            precursorChargeMap.put(scanName, precursorCharge);
            precursorMassMap.put(scanName, precursorMass);
            submitNum++;
            validScanSet.add(scanName);
            taskListGetLongTag.add(threadPoolGetLongTag.submit(new GetLongTag(scanNum, buildIndex, massTool, spectraParserArray[fileId], lockGetLongTag, scanName, precursorCharge
                    , precursorMass, specProcessor, ms2Tol)));
        }
        logger.debug("totalSubmit in GetLongTag, "+ submitNum);
        sqlResSetGetLongTag.close();
        sqlStateGetLongTag.close();

        int lastProgressGetLongTag = 0;
        int totalCountGetLongTag = taskListGetLongTag.size();
        int countGetLongTag = 0;
        Map<String, List<GetLongTag.TagRes>> prot_TagResList_Map = new HashMap<>();
        while (countGetLongTag < totalCountGetLongTag) {
            // record search results and delete finished ones.
            List<Future<GetLongTag.Entry>> toBeDeleteTaskList = new ArrayList<>(totalCountGetLongTag - countGetLongTag);
            for (Future<GetLongTag.Entry> task : taskListGetLongTag) {
                if (task.isDone()) {
                    if (task.get() != null ) {
                        GetLongTag.Entry entry = task.get();
                        for (String prot : entry.prot_TagResList_Map.keySet()){
                            List<GetLongTag.TagRes> tagResList = entry.prot_TagResList_Map.get(prot);
                            if  (prot_TagResList_Map.containsKey(prot)) {
                                prot_TagResList_Map.get(prot).addAll(tagResList);
                            } else {
                                List<GetLongTag.TagRes> tmpTagResList = new LinkedList<>();
                                tmpTagResList.addAll(tagResList);
                                prot_TagResList_Map.put(prot, tagResList);
                            }
                        }
                    }

                    toBeDeleteTaskList.add(task);
                } else if (task.isCancelled()) {
                    toBeDeleteTaskList.add(task);
                }
            }
            countGetLongTag += toBeDeleteTaskList.size();
            taskListGetLongTag.removeAll(toBeDeleteTaskList);
            taskListGetLongTag.trimToSize();

            int progress = countGetLongTag * 20 / totalCountGetLongTag;
            if (progress != lastProgressGetLongTag) {
                logger.info("Getting long tags for prot {}%...", progress * 5);
                lastProgressGetLongTag = progress;
            }

            if (countGetLongTag == totalCountGetLongTag) {
                break;
            }
            Thread.sleep(6000);
        }
        // shutdown threads.
        threadPoolGetLongTag.shutdown();
        if (!threadPoolGetLongTag.awaitTermination(60, TimeUnit.SECONDS)) {
            threadPoolGetLongTag.shutdownNow();
            if (!threadPoolGetLongTag.awaitTermination(60, TimeUnit.SECONDS))
                throw new Exception("Pool did not terminate");
        }
        if (lockGetLongTag.isLocked()) {
            lockGetLongTag.unlock();
        }
        //   ==========END=============Get Long Tags

        //////////==================================================
        //   Score for Proteins
        Map<String, Integer> protLengthMap = buildIndex.protLengthMap;
        Map<String, Double> protCoverageMap = new HashMap<>();
        for (String prot : prot_TagResList_Map.keySet()){
            int protLen = buildIndex.protSeqMap.get(prot).length();
            List<GetLongTag.TagRes> tagResList = prot_TagResList_Map.get(prot);
            Map<Integer, Double> posProductMap = new HashMap<>(protLen);
            for (int pos = 0; pos < protLen; pos++) {
                posProductMap.put(pos, 1.0); //intial, all pos is with 1
            }
            for (GetLongTag.TagRes tagRes : tagResList) {
                int relPos = tagRes.relPos;
                List<Double> normedIaaList = tagRes.normedIaaList;
                for (int aaPos = relPos; aaPos < relPos+tagRes.tagSeq.length(); aaPos++) {
                    posProductMap.put(aaPos, posProductMap.get(aaPos)*(1- normedIaaList.get(aaPos-relPos)));
                }
            } // finish aa cal
            double tempCoverage = 0;
            for (int pos : posProductMap.keySet()) {
                tempCoverage += 1-posProductMap.get(pos);
            }
            protCoverageMap.put(prot, tempCoverage/protLen);
        }

        List<Pair<String, Double>> protScoreLongList = new ArrayList<>();
        for (String prot : protCoverageMap.keySet()){
            protScoreLongList.add(new Pair<>(prot, protCoverageMap.get(prot)));
        }
        Collections.sort(protScoreLongList, Comparator.comparing(o -> o.getSecond(), Comparator.reverseOrder()));

        Set<String> reducedProtIdSet = new HashSet<>();

        for (Pair<String, Double> pair : protScoreLongList){
            if (pair.getSecond() < proteinCovThres) break;
            reducedProtIdSet.add(pair.getFirst());
        }

        if (isDebugMode) reducedProtIdSet = protLengthMap.keySet();  //dont reduce for debug
//        reducedProtIdSet = protLengthMap.keySet();  //dont reduce

        logger.debug("Protein size: "+protLengthMap.size()+" reduced to " + reducedProtIdSet.size());
//        for (String protId : reducedProtIdSet) {
//            System.out.println(protId);
//        }
        //  END Score for Proteins

        //////////==================================================
        //   Get Tags and Get Peptide Candidates
        logger.info("Building decoy...");
        String dbPath = parameterMap.get("database");
        String decoyFullDBPath = dbPath + ".decoy.fasta";
        File decoyFullFile = new File(decoyFullDBPath);
        if (decoyFullFile.exists()) {
            DbTool decoyDbTool = new DbTool(decoyFullDBPath, "others");
            Map<String, String> decoyProtSeqMap = decoyDbTool.getProtSeqMap();
            for (String decoyProtId : decoyProtSeqMap.keySet()) {
                String targetProtId = decoyProtId.split("DECOY_")[1];
                if (reducedProtIdSet.contains(targetProtId)) {
                    buildIndex.protSeqMap.put(decoyProtId, decoyProtSeqMap.get(decoyProtId)); // using the target sequence to replace contaminant sequence if there is conflict.
                }
            }
        } else {
            BufferedWriter decoyWriter = new BufferedWriter(new FileWriter(decoyFullDBPath));
//            int threadNum_1 = Integer.valueOf(parameterMap.get("thread_num"));
//            if (threadNum_1 == 0)  threadNum_1 = 3 + Runtime.getRuntime().availableProcessors();
            ExecutorService threadPoolBuildDecoyProts = Executors.newFixedThreadPool(p_thread_num);
            ArrayList<Future<BuildDecoyProts.Entry>> taskListBuildDecoyProts = new ArrayList<>(reducedProtIdSet.size() + 10);
            ReentrantLock lockSpecCoder = new ReentrantLock();
            for (String protId : protLengthMap.keySet()) {    //All protSeq are recorded but some decoy version are not because their targets are not in the reducedProtIdSet
                taskListBuildDecoyProts.add(threadPoolBuildDecoyProts.submit(new BuildDecoyProts(parameterMap, buildIndex, protId)));
            }
            int lastProgressBuildDecoyProts = 0;
            int totalCountBuildDecoyProts = taskListBuildDecoyProts.size();
            int countBuildDecoyProts = 0;
            while (countBuildDecoyProts < totalCountBuildDecoyProts) {
                List<Future<BuildDecoyProts.Entry>> toBeDeleteTaskList = new ArrayList<>(totalCountBuildDecoyProts - countBuildDecoyProts);
                for (Future<BuildDecoyProts.Entry> task : taskListBuildDecoyProts) {
                    if (task.isDone()) {
                        if (task.get() != null ) {
                            BuildDecoyProts.Entry entry = task.get();
                            String protId = entry.protId;
                            String decoyProtSeq = entry.decoyProtSeq;
                            String decoyProtId = "DECOY_"+protId;
                            if (reducedProtIdSet.contains(protId)) buildIndex.protSeqMap.put(decoyProtId, decoyProtSeq);
                            decoyWriter.write(String.format(Locale.US, ">%s\n", decoyProtId));
                            decoyWriter.write(decoyProtSeq + "\n");
                        }
                        toBeDeleteTaskList.add(task);
                    } else if (task.isCancelled()) {
                        toBeDeleteTaskList.add(task);
                    }
                }
                countBuildDecoyProts += toBeDeleteTaskList.size();
                taskListBuildDecoyProts.removeAll(toBeDeleteTaskList);
                taskListBuildDecoyProts.trimToSize();

                int progress = countBuildDecoyProts * 20 / totalCountBuildDecoyProts;
                if (progress != lastProgressBuildDecoyProts) {
                    logger.info("Build Decoy Prots {}%...", progress * 5);
                    lastProgressBuildDecoyProts = progress;
                }

                if (countBuildDecoyProts == totalCountBuildDecoyProts) {
                    break;
                }
                Thread.sleep(6000);
            }
            // shutdown threads.
            threadPoolBuildDecoyProts.shutdown();
            if (!threadPoolBuildDecoyProts.awaitTermination(60, TimeUnit.SECONDS)) {
                threadPoolBuildDecoyProts.shutdownNow();
                if (!threadPoolBuildDecoyProts.awaitTermination(60, TimeUnit.SECONDS))
                    throw new Exception("Pool did not terminate");
            }
            if (lockSpecCoder.isLocked()) {
                lockSpecCoder.unlock();
            }
            decoyWriter.close();
        }

        // reduce proteins from buildIndex if it is not contained in reducedProtIdSet
//        Set<String> decoyProtIdToRemove = new HashSet<>();
        Iterator<String> iter = buildIndex.protSeqMap.keySet().iterator();
        while (iter.hasNext()) {
//            String protId = iter.next();
            String pureId = iter.next().replace("DECOY_","");
            if (!reducedProtIdSet.contains(pureId)) {
//                decoyProtIdToRemove.add("DECOY_"+protId);
                iter.remove();
            }
        }
        logger.debug("protTD size "+ buildIndex.protSeqMap.size());
//        for (String protId : buildIndex.protSeqMap.keySet()) {
//            System.out.println(protId);
//        }
        logger.info("Generating reduced FM index...");
        buildBiDirectionFMIndex(buildIndex);
        buildIndex.minPeptideMass = 0;
        buildIndex.maxPeptideMass = 9999;
        // writer concatenated fasta
        Map<String, String> proteinAnnotationMap;

        DbTool dbTool = buildIndex.dbTool;
        DbTool contaminantsDb;
        if (addContaminants) { //addContaminants = true
            contaminantsDb = new DbTool(null, "contaminants");
            proteinAnnotationMap = contaminantsDb.getProteinAnnotateMap();
            proteinAnnotationMap.putAll(dbTool.getProteinAnnotateMap()); // using the target annotation to replace contaminant sequence if there is conflict.
        } else {
            proteinAnnotationMap = dbTool.getProteinAnnotateMap();
        }
        BufferedWriter writer = new BufferedWriter(new FileWriter(dbPath + ".TD.fasta"));
        for (String protId : buildIndex.protSeqMap.keySet()) {
            writer.write(String.format(Locale.US, ">%s %s\n", protId, proteinAnnotationMap.getOrDefault(protId, "")));
            writer.write(buildIndex.protSeqMap.get(protId) + "\n");
        }
        writer.close();

        logger.info("Main searching...");

        p_thread_num = Integer.valueOf(parameterMap.get("PIPI_thread_num"));
        if (Integer.valueOf(parameterMap.get("PIPI_thread_num")) == 0) {
            p_thread_num = (int) (Runtime.getRuntime().availableProcessors()/2);
        }
        if (isDebugMode) p_thread_num = 1;
        logger.debug("availabel processors for Main search"+ Runtime.getRuntime().availableProcessors());
        logger.debug("thread num for Main search"+ p_thread_num);
        ExecutorService threadPoolBone = Executors.newFixedThreadPool(p_thread_num);
        ArrayList<Future<MainSearch.Entry>> taskListBone = new ArrayList<>(validScanSet.size() + 10);
        ReentrantLock lockBone = new ReentrantLock();
        int submitNumBone = 0;
        for (String scanName : validScanSet) {
            String[] scanNameStr = scanName.split("\\.");
            int precursorCharge = precursorChargeMap.get(scanName);

            double precursorMass = precursorMassMap.get(scanName);
            int scanNum = Integer.valueOf(scanNameStr[2]);
            boolean shouldRun = false;
            if (isDebugMode) {
                for (int debugScanNum : lszDebugScanNum) {
                    if (Math.abs(scanNum-debugScanNum) < 1) {
                        shouldRun = true;
                    }
                }
                if (!shouldRun) continue;// skip scans not in the debug list
            }
            submitNumBone++;
            int fileId = fileNameIdMap.get( scanNameStr[0] );
            taskListBone.add(threadPoolBone.submit(new MainSearch(scanNum, buildIndex, massTool, ms2Tol, ms1Tol, inferPTM.getMinPtmMass(), inferPTM.getMaxPtmMass()
                    , spectraParserArray[fileId], lockBone, scanName, precursorCharge, precursorMass, specProcessor )));
        }
        logger.info(submitNumBone+ " MS2 submitted to Main Search.");

        Map<String, PeptideInfo> allPeptideInfoMap = new HashMap<>();
        int lastProgressBone = 0;
        int resultCountBone = 0;
        int totalCountBone = taskListBone.size();
        int countBone = 0;

        Map<Integer, TreeMap<Double, Set<String>>> fileId_pcMassScanNameMap = new HashMap<>();
        for (int fileId : fileIdNameMap.keySet()) {
            fileId_pcMassScanNameMap.put(fileId, new TreeMap<>());
        }
        Map<String, Peptide> scanName_TopPeptide_Map = new HashMap<>();
        Map<String, String> scanName_PepString_Map = new HashMap<>();
        Connection sqlConSpecCoder = DriverManager.getConnection(sqlPath);
        PreparedStatement sqlPreparedStatement = sqlConSpecCoder.prepareStatement("REPLACE INTO spectraTable (scanNum, scanName,  precursorCharge, precursorMass, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient,  precursorScanNo, labelling, peptide, theoMass, score,peptideSet) VALUES (?, ?, ?,?, ?, ?, ?, ?, ?, ?, ?, ?, ?)");
        sqlConSpecCoder.setAutoCommit(false);
        Map<VarPtm, Integer> varPtmCountMap = new HashMap<>();
        while (countBone < totalCountBone) {
            // record search results and delete finished ones.
            List<Future<MainSearch.Entry>> toBeDeleteTaskList = new ArrayList<>(totalCountBone - countBone);
            for (Future<MainSearch.Entry> task : taskListBone) {
                if (task.isDone()) {
                    if (task.get() != null) {
                        MainSearch.Entry entry = task.get();
                        for (VarPtm varPtm : entry.varPtmList){
                            if (varPtmCountMap.containsKey(varPtm)) {
                                varPtmCountMap.put(varPtm, varPtmCountMap.get(varPtm)+1);
                            }else{
                                varPtmCountMap.put(varPtm, 1);
                            }
                        }
                        allPeptideInfoMap.putAll(entry.peptideInfoMapForRef);
                        sqlPreparedStatement.setInt(1, entry.scanNum);
                        sqlPreparedStatement.setString(2, entry.scanName);
                        sqlPreparedStatement.setInt(3, entry.precursorCharge);
                        sqlPreparedStatement.setDouble(4, entry.precursorMass);
                        sqlPreparedStatement.setString(5, mgfTitleMap.get(entry.scanName));
                        sqlPreparedStatement.setInt(6, isotopeCorrectionNumMap.get(entry.scanName));
                        sqlPreparedStatement.setDouble(7, ms1PearsonCorrelationCoefficientMap.get(entry.scanName));
                        sqlPreparedStatement.setInt(8, precursorScanNoMap.get(entry.scanName));
                        sqlPreparedStatement.setString(9, entry.labelling);
                        sqlPreparedStatement.setString(10, entry.peptide);
                        sqlPreparedStatement.setDouble(11, entry.theoMass);
                        sqlPreparedStatement.setDouble(12, entry.score);
                        sqlPreparedStatement.setString(13, entry.peptideSet);

                        sqlPreparedStatement.executeUpdate();
                        int fileId = fileNameIdMap.get( entry.scanName.split("\\.")[0] );
                        TreeMap<Double,Set<String>> local_pcMassScanNameMap = fileId_pcMassScanNameMap.get(fileId);
                        if (local_pcMassScanNameMap.containsKey(entry.precursorMass)) {
                            local_pcMassScanNameMap.get(entry.precursorMass).add(entry.scanName);
                        }else {
                            Set<String> scanNumSet = new HashSet<>();
                            scanNumSet.add(entry.scanName);
                            local_pcMassScanNameMap.put(entry.precursorMass, scanNumSet);
                        }
                        scanName_TopPeptide_Map.put(entry.scanName, entry.topPeptide);
                        scanName_PepString_Map.put(entry.scanName, entry.peptideSet);

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

            sqlConSpecCoder.commit();
            int progress = countBone * 20 / totalCountBone;
            if (progress != lastProgressBone) {
                logger.info("Main Search {}%...", progress * 5);
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
        sqlConSpecCoder.commit();
        sqlConSpecCoder.setAutoCommit(true);
        sqlConSpecCoder.close();
        if (lockBone.isLocked()) {
            lockBone.unlock();
        }
        logger.info(resultCountBone + " MS2 finished in Main Search.");

        logger.info("Supplementary searching...");
        p_thread_num = 3 + Runtime.getRuntime().availableProcessors();
        ExecutorService threadPoolSupp = Executors.newFixedThreadPool(p_thread_num);
        ArrayList<Future<SuppSearch.Entry>> taskListSupp = new ArrayList<>(validScanSet.size() + 10);
        ReentrantLock lockSupp = new ReentrantLock();
        Binomial binomial = new Binomial(Integer.valueOf(parameterMap.get("max_peptide_length")) * 2);
        int submitTimeSupp = 0;

        int n_complementaryPeptides = 5;
        for (String thisScanName : validScanSet) {
            int thisFileId = fileNameIdMap.get( thisScanName.split("\\.")[0] );
            double mass = precursorMassMap.get(thisScanName);
            int thisScanNum = Integer.valueOf( thisScanName.split("\\.")[2] );

            String thisPtmSeq = "XFAKEPEP";
            if (scanName_TopPeptide_Map.containsKey(thisScanName)) {
                thisPtmSeq = scanName_TopPeptide_Map.get(thisScanName).getVarPtmContainingSeqNow();
            }
            Map<String, TreeMap<Integer, VarPtm>> ptmSeq_posVarPtmMap_Map = new HashMap<>();
            Map<String, PeptideInfo> local_PepSeq_PeptideInfo_Map = new HashMap<>();

            TreeMap<Double, Set<String>> local_other_pcMassScanNameMap = fileId_pcMassScanNameMap.get(thisFileId);
            for (int i = 0; i <= 0; i++) {
                for (Set<String> otherScanNameSet : local_other_pcMassScanNameMap.subMap(mass + i*MassTool.PROTON - mass*ms1Tol*1e-6, true, mass + i*MassTool.PROTON + mass*ms1Tol*1e-6, true).values()) {
                    for (String otherScanName : otherScanNameSet) {
                        int otherScanNum = Integer.valueOf( otherScanName.split("\\.")[2] );
                        if (otherScanNum < thisScanNum+2000 && otherScanNum > thisScanNum-2000 && otherScanNum != thisScanNum) {
                            String otherPtmSeq = scanName_TopPeptide_Map.get(otherScanName).getVarPtmContainingSeqNow();
                            String otherFreeSeq = scanName_TopPeptide_Map.get(otherScanName).getFreeSeq();
                            if (otherPtmSeq.contentEquals(thisPtmSeq)) continue;

                            if (!ptmSeq_posVarPtmMap_Map.containsKey(otherPtmSeq)) {
                                ptmSeq_posVarPtmMap_Map.put(otherPtmSeq, scanName_TopPeptide_Map.get(otherScanName).posVarPtmResMap);
                                local_PepSeq_PeptideInfo_Map.put(otherFreeSeq, allPeptideInfoMap.get(otherFreeSeq).clone());
                            }
                        }
                    }
                }
            }

            if (ptmSeq_posVarPtmMap_Map.isEmpty()) {
                continue;
            }
            int precursorCharge = precursorChargeMap.get(thisScanName);
            double precursorMass = precursorMassMap.get(thisScanName);
            submitTimeSupp++;
            taskListSupp.add(threadPoolSupp.submit(new SuppSearch(thisScanNum, buildIndex, massTool
                    , ms2Tol, spectraParserArray[thisFileId], lockSupp, thisScanName, precursorCharge
                    , precursorMass, specProcessor, binomial, ptmSeq_posVarPtmMap_Map, local_PepSeq_PeptideInfo_Map)));
        }
        logger.info(submitTimeSupp + " MS2 submitted to Supplementary Search" + submitTimeSupp);
        // check progress every minute, record results,and delete finished tasks.
        int lastProgressSupp = 0;
        int resultCountSupp = 0;
        int totalCountSupp = taskListSupp.size();
        int countSupp = 0;
        Connection sqlConSupp = DriverManager.getConnection(sqlPath);
        PreparedStatement sqlPreparedStatementPTM = sqlConSupp.prepareStatement("REPLACE INTO spectraTable (scanNum, scanName,  precursorCharge, precursorMass, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient,  precursorScanNo, labelling, peptide, theoMass, score, peptideSet) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)");
        sqlConSupp.setAutoCommit(false);
        while (countSupp < totalCountSupp) {
            List<Future<SuppSearch.Entry>> toBeDeleteTaskList = new ArrayList<>(totalCountSupp - countSupp);
            for (Future<SuppSearch.Entry> task : taskListSupp) {
                if (task.isDone()) {
                    if (task.get() != null) {
                        SuppSearch.Entry entry = task.get();
                        if (  scanName_TopPeptide_Map.containsKey(entry.scanName)
                                && entry.score < scanName_TopPeptide_Map.get(entry.scanName).getScore()) {
                            toBeDeleteTaskList.add(task);
                            continue; // if the complementary is not better dont change in sql
                        }

                        sqlPreparedStatementPTM.setInt(1, entry.scanNum);
                        sqlPreparedStatementPTM.setString(2, entry.scanName);
                        sqlPreparedStatementPTM.setInt(3, entry.precursorCharge);
                        sqlPreparedStatementPTM.setDouble(4, entry.precursorMass);
                        sqlPreparedStatementPTM.setString(5, mgfTitleMap.get(entry.scanName));
                        sqlPreparedStatementPTM.setInt(6, isotopeCorrectionNumMap.get(entry.scanName));
                        sqlPreparedStatementPTM.setDouble(7, ms1PearsonCorrelationCoefficientMap.get(entry.scanName));
                        sqlPreparedStatementPTM.setInt(8, precursorScanNoMap.get(entry.scanName));
                        sqlPreparedStatementPTM.setString(9, entry.labelling);
                        sqlPreparedStatementPTM.setString(10, entry.peptide);
                        sqlPreparedStatementPTM.setDouble(11, entry.theoMass);
                        sqlPreparedStatementPTM.setDouble(12, entry.score);
                        String peptideSetStr = null;
                        if (  scanName_PepString_Map.containsKey(entry.scanName) ){
                            peptideSetStr = entry.peptideSet + "," + scanName_PepString_Map.get(entry.scanName);
                        } else {
                            peptideSetStr = entry.peptideSet;
                        }
                        sqlPreparedStatementPTM.setString(13, peptideSetStr); //if a better complementary pep comes, need to insert it in the front

                        sqlPreparedStatementPTM.executeUpdate();
                        ++resultCountSupp;
                    }

                    toBeDeleteTaskList.add(task);
                } else if (task.isCancelled()) {
                    toBeDeleteTaskList.add(task);
                }
            }
            countSupp += toBeDeleteTaskList.size();
            taskListSupp.removeAll(toBeDeleteTaskList);
            taskListSupp.trimToSize();

            sqlConSupp.commit();

            int progress = countSupp * 20 / totalCountSupp;
            if (progress != lastProgressSupp) {
                logger.info("Supplementary Searching {}%...", progress * 5);
//                System.out.println(toBeDeleteTaskList.size()+ ","+ taskListSupp.size() +"," + count2 + ", " + totalCount + "," + lastProgress);
                lastProgressSupp = progress;
            }

            if (countSupp == totalCountSupp) {
                break;
            }
            Thread.sleep(6000);
        }
        // shutdown threads.
        threadPoolSupp.shutdown();
        if (!threadPoolSupp.awaitTermination(60, TimeUnit.SECONDS)) {
            threadPoolSupp.shutdownNow();
            if (!threadPoolSupp.awaitTermination(60, TimeUnit.SECONDS))
                throw new Exception("Pool did not terminate");
        }

        sqlConSupp.commit();
        sqlConSupp.setAutoCommit(true);
        sqlConSupp.close();
        if (lockSupp.isLocked()) {
            lockSupp.unlock();
        }

        postProcessing(outputDir, allPeptideInfoMap, sqlPath, buildIndex.protSeqMap, massTool, hostName, varPtmCountMap);
        logger.info("Saving results...");
    }

    private void buildBiDirectionFMIndex(BuildIndex buildIndex) throws IOException {
        BufferedWriter writerProtNormal  = new BufferedWriter(new FileWriter("catProtNormal.txt"));
        BufferedWriter writerProtReverse = new BufferedWriter(new FileWriter("catProtReverse.txt"));
        StringBuilder normalSb = new StringBuilder();
        int dotPos = 0;
        int dotNum = 0;
        buildIndex.dotPosArrNormal = new int[buildIndex.protSeqMap.keySet().size()];
        for (String protId : buildIndex.protSeqMap.keySet()) {
            buildIndex.dotPosArrNormal[dotNum] = dotPos;
            buildIndex.posProtMapNormal.put(dotNum, protId);
            String protSeq = buildIndex.protSeqMap.get(protId).replace('I', 'L');
            buildIndex.protSeqMap.put(protId, protSeq);
            normalSb.append(".");
            normalSb.append(protSeq.replace('I', 'L'));
            dotNum++;
            dotPos += protSeq.length()+1;
        }
        writerProtNormal.write(normalSb.toString());
        writerProtReverse.write(normalSb.reverse().toString());
        writerProtNormal.close();
        writerProtReverse.close();
        char[] textNormal = buildIndex.loadFile("catProtNormal.txt", true);
        char[] textReverse = buildIndex.loadFile("catProtReverse.txt", true);
        buildIndex.fmIndexNormal = new FMIndex(textNormal);
        buildIndex.fmIndexReverse = new FMIndex(textReverse);
        buildIndex.textNormalLen = textNormal.length-1; //remove the \r\n
    }
    private boolean isTarget(Collection<String> proteinIds) {
        for (String protein : proteinIds) {
            if (!protein.startsWith("DECOY_")) {
                return true;
            }
        }
        return false;
    }
    private static void help() {
        String helpStr = "PIPI version " + versionStr + "\r\n"
                + "Identification of Peptide with Multiple PTMs using Combinatorial Optimization.\r\n"
                + "Author: Shengzhi Lai\r\n"
                + "Email: slaiad@connect.ust.hk\r\n"
                + "java -Xmx8g -jar PIPI.jar parameter.def spectra_file output_directory\r\n"
                + "More details at: https://bioinformatics.hkust.edu.hk/";
        System.out.print(helpStr);
        System.exit(1);
    }

    class ScanRes {
        public int scanNum;
        public double qValue = -0.1;
        public List<CandiScore> peptideInfoScoreList;
        public double expMass;
        public int charge;
        public String scanName;
        ScanRes(String scanName, int scanNum, double expMass, List<CandiScore> peptideInfoScoreList, int charge){
            this.scanName = scanName;
            this.scanNum = scanNum;
            this.expMass = expMass;
            this.peptideInfoScoreList = peptideInfoScoreList;
            this.charge = charge;
        }
    }
    class CandiScore implements Comparable<CandiScore>{
        public String ptmContainingSeq;
        public PeptideInfo peptideInfo;
        public double pepScore;
        public double protScore = 0;
        public double varPtmTotalScore = 0;
        CandiScore(PeptideInfo peptideInfo, double pepScore, String ptmContainingSeq) {
            this.peptideInfo = peptideInfo;
            this.pepScore = pepScore;
            this.ptmContainingSeq = ptmContainingSeq;
        }

        public void setVarPtmTotalScore(Map<String, Double> varPtmStrScoreRefMap) { //varPtmScoreRefMap is like K(14.016)->1.2
            int startI = -1;
            double varPtmTotalScore = 0;
            int n_PtmOnPep = 0;
            for (int anyI = 0; anyI < ptmContainingSeq.length(); anyI++) {
                char thisChar = ptmContainingSeq.charAt(anyI);
                if (thisChar == '(') {
                    n_PtmOnPep++;
                    startI = anyI-1;
                } else if (thisChar == ')') {
                    String thisVarPtmStr = ptmContainingSeq.substring(startI, anyI + 1);
                    if (varPtmStrScoreRefMap.containsKey(thisVarPtmStr)) {
                        varPtmTotalScore += varPtmStrScoreRefMap.get(thisVarPtmStr);
                    }
                }
            }
            this.varPtmTotalScore = varPtmTotalScore/n_PtmOnPep;
        }
        public int compareTo(CandiScore o2) {
            if (this.protScore < o2.protScore) {
                return -1;
            } else if (this.protScore > o2.protScore) {
                return 1;
            } else {
                if (this.varPtmTotalScore < o2.varPtmTotalScore) {
                    return -1;
                } else if (this.varPtmTotalScore > o2.varPtmTotalScore) {
                    return 1;
                } else {
                    if (this.pepScore < o2.pepScore) {
                        return -1;
                    } else if (this.pepScore > o2.pepScore) {
                        return 1;
                    }
                }
            }
            return 0;
        }
    }
    private void postProcessing(String outputDir, Map<String, PeptideInfo> allPeptideInfoMap, String sqlPath, Map<String, String> protSeqMap, MassTool massTool, String hostName, Map<VarPtm, Integer> varPtmCountMap) throws IOException, SQLException {
        //preprocess varPtmCountMap
        Map<String, Double> varPtmRefScoreMap = new HashMap<>();
        for (VarPtm varPtm : varPtmCountMap.keySet()) {
            if (varPtmCountMap.get(varPtm) == 1) {
                continue; // no need to care for those ptm only appears once because their log10 is naturally 0
            }
            varPtmRefScoreMap.put(varPtm.site+"("+InferPTM.df3.format(varPtm.mass)+")", Math.sqrt(varPtmCountMap.get(varPtm)));
        }
        System.out.println(varPtmCountMap.size());
        System.out.println(varPtmRefScoreMap.size());
        List<Map.Entry<String, Double>> testList = new ArrayList<>(varPtmRefScoreMap.entrySet());
        Collections.sort(testList, Map.Entry.comparingByValue(Comparator.reverseOrder()));
        int j = 0;
        for (Map.Entry<String, Double> entry : testList) {
            String varPtmStr = entry.getKey();
//            System.out.println(varPtmStr+","+InferPTM.df3.format(entry.getValue()) + "," + Math.pow(entry.getValue(),2));
            if (j>maxNumVarPtmConsidered) {  // this 5 is too small for PT09449 where 17 modification is considered
//                System.out.println("removed, "+ varPtmStr);
                varPtmRefScoreMap.remove(varPtmStr);
            }
            j++;
        }

        for (String varPtmStr : varPtmRefScoreMap.keySet()) {
            System.out.println("after,"+varPtmStr + "," + InferPTM.df3.format(varPtmRefScoreMap.get(varPtmStr)));
        }
        //collect data
        Connection sqlConnection = DriverManager.getConnection(sqlPath);
        Statement sqlStatement = sqlConnection.createStatement();
        ResultSet sqlResultSet = sqlStatement.executeQuery("SELECT scanNum, scanName, precursorCharge, precursorMass, peptide, theoMass, score, peptideSet FROM spectraTable");
        Map<String, Map<String, Double>> protPepScoreMap = new HashMap<>();

        List<ScanRes> scanResList = new ArrayList<>();
        List<Double> topScoreList = new LinkedList<>();
        while (sqlResultSet.next()) {
            String topPeptide = sqlResultSet.getString("peptide"); //this is ptmContaining Seq without n or c
            if (!sqlResultSet.wasNull()) {
                int charge = sqlResultSet.getInt("precursorCharge");
//                double theoMass = sqlResultSet.getDouble("theoMass");
                double expMass = sqlResultSet.getDouble("precursorMass");
//                double massDiff = getMassDiff(expMass, theoMass, MassTool.C13_DIFF);
                double score = sqlResultSet.getDouble("score");
                String peptideSet = sqlResultSet.getString("peptideSet");
                int scanNum = sqlResultSet.getInt("scanNum");
                String scanName = sqlResultSet.getString("scanName");

                String[] candiSetStr = peptideSet.split(",");
                int numPep = candiSetStr.length/3;

                PeptideInfo pepInfo = allPeptideInfoMap.get(topPeptide.replaceAll("[^A-Z]+", ""));
                List<CandiScore> candiScoreList = new ArrayList<>();
                for (int i = 0; i < numPep; i++) {
                    //dont put in the impossible ones
                    String ptmContainingSeq = candiSetStr[3*i+0];
                    PeptideInfo candiPeptideInfo = allPeptideInfoMap.get(ptmContainingSeq.replaceAll("[^A-Z]+", ""));
                    double thisScore = Double.valueOf(candiSetStr[3*i+1]);
//                    if (lszDebugScanNum.contains(scanNum) ){
//                        System.out.println(ptmContainingSeq + ", "+ thisScore);
//                    }
                    CandiScore candiScore = new CandiScore(candiPeptideInfo, thisScore, ptmContainingSeq);
                    candiScore.setVarPtmTotalScore(varPtmRefScoreMap);
                    candiScoreList.add(candiScore); //peptideInfo and their score
                }
                Collections.sort(candiScoreList, Comparator.comparing(o -> o.pepScore, Comparator.reverseOrder()));
                double topPepScore = candiScoreList.get(0).pepScore;
//                System.out.println(scanNum + "," + topPepScore + "," + (candiScoreList.get(0).peptideInfo.isDecoy ? 1 : 0)); // print all top candidates scores
                topScoreList.add(topPepScore);
                Iterator<CandiScore> iter = candiScoreList.iterator();
                while (iter.hasNext()) {
                    CandiScore candiScore = iter.next();
                    if (candiScore.pepScore < 0.85 * topPepScore) { // candidates with score < 0.5* topscore are not consider in PFM
                        iter.remove();
                    }
                }
                scanResList.add(new ScanRes(scanName, scanNum, expMass, candiScoreList, charge));

                for (String protId : pepInfo.protIdSet){
                    if (protPepScoreMap.containsKey(protId)){
                        Map<String, Double> pepScoreMap = protPepScoreMap.get(protId);
                        if (pepScoreMap.containsKey(topPeptide)){
                            pepScoreMap.put(topPeptide, Math.max(pepScoreMap.get(topPeptide), score)); // use max score of repeated peptides for any protein
                        } else {
                            pepScoreMap.put(topPeptide, score);
                        }
                    } else {
                        Map<String, Double> pepScoreMap = new HashMap<>();
                        pepScoreMap.put(topPeptide, score);
                        protPepScoreMap.put(protId, pepScoreMap);
                    }
                }
            }
        }
        sqlResultSet.close();
        sqlStatement.close();
        sqlConnection.close();

        Collections.sort(topScoreList);
        int numScores = topScoreList.size();
        double medianScore = numScores%2 == 0 ? (topScoreList.get(numScores/2) + topScoreList.get(numScores/2-1))/2 : topScoreList.get((numScores-1)/2);
        medianScore*=1.5;
        //calculate prot score
        Map<String, Double> protScoreMap = new HashMap<>();
        for (String protId : protPepScoreMap.keySet()) {
            Map<String,Double> pepScoreMap = protPepScoreMap.get(protId);
            for (String pep : pepScoreMap.keySet()){
                if (pepScoreMap.get(pep) > medianScore) {     //this peptide score threshold is empirical
                    if (protScoreMap.containsKey(protId)){
                        protScoreMap.put(protId, protScoreMap.get(protId)+ pepScoreMap.get(pep));
                    } else {
                        protScoreMap.put(protId, pepScoreMap.get(pep));
                    }
                }
            }
        }
        //normalize prot score
        for (String protId : protScoreMap.keySet()){
            protScoreMap.put(protId, protScoreMap.get(protId) / Math.log(protSeqMap.get(protId).length()));
//            System.out.println(protId + "," + protScoreMap.get(protId) + "," + protSeqMap.get(protId).length()); // print all top candidates scores
        }
        //===============================

        //update prot score for candidates
        for (ScanRes scanRes : scanResList) {
            List<CandiScore> candiScoreList = scanRes.peptideInfoScoreList;
            for (CandiScore candiScore : candiScoreList) {
                double protScoreForCand = -1;
                for (String protId : candiScore.peptideInfo.protIdSet) {
                    if (protScoreMap.getOrDefault(protId, 0.0) > protScoreForCand){
                        protScoreForCand = protScoreMap.getOrDefault(protId, 0.0);
                    }
                }
                candiScore.protScore = protScoreForCand;
            }
        }

        DecimalFormat df= new  DecimalFormat( ".00000" );

        //calculate ori FDR
        for (ScanRes scanRes : scanResList) {
            Collections.sort(scanRes.peptideInfoScoreList, Comparator.comparing(o -> o.pepScore, Comparator.reverseOrder())); // rank candidates using peptide score
        }
        Collections.sort(scanResList, Comparator.comparing(o -> o.peptideInfoScoreList.get(0).pepScore, Comparator.reverseOrder()));
        List<Double> fdrList = new ArrayList<>(scanResList.size());
        int numTPlusD = 0;
        int numD = 0;
        for (ScanRes scanRes : scanResList) {  //scanResList is already ranked
            CandiScore candiScore = scanRes.peptideInfoScoreList.get(0);
            numTPlusD++;
            if (candiScore.peptideInfo.isDecoy) numD++;
            fdrList.add(2*(double)numD/numTPlusD);
//            System.out.println(scanRes.scanNum + "," + candiScore.pepScore + "," +(candiScore.peptideInfo.isDecoy ? 0 : 1));
        }
        int numQ001 = 0;
        double minQ = 1.0;
        boolean found = false;
        for (int i = fdrList.size()-1; i >= 0; i--) {
            minQ = Math.min(fdrList.get(i), minQ);
            scanResList.get(i).qValue = minQ;
            if (!found && minQ < 0.01) {
                numQ001 = i;
                found = true;
//                break;
            }
        }
        BufferedWriter oriWriter = new BufferedWriter(new FileWriter(outputDir+"Res."+hostName+".Peptides.csv"));
        oriWriter.write("scanNum,qValue,TorD,peptide,pepScore,proteins,protscore,peptide,pepScore,proteins,protscore,peptide,pepScore,proteins,protscore,peptide,pepScore,proteins,protscore,peptide,pepScore,proteins,protscore,peptide,pepScore,proteins,protscore,peptide,pepScore,proteins,protscore,peptide,pepScore,proteins,protscore,peptide,pepScore,proteins,protscore,peptide,pepScore,proteins,protscore,peptide,pepScore,proteins,protscore,peptide,pepScore,proteins,protscore,peptide,pepScore,proteins,protscore\n");
        for (ScanRes scanRes : scanResList) {
            List<CandiScore> candiScoreList = scanRes.peptideInfoScoreList;
            StringBuilder str = new StringBuilder();
            str.append(scanRes.scanNum+",").append(scanRes.qValue+",").append(candiScoreList.get(0).peptideInfo.isDecoy ? 0 : 1);
            for (CandiScore candiScore : candiScoreList){
                str.append(","+candiScore.peptideInfo.freeSeq).append(","+candiScore.pepScore).append(","+String.join(";", candiScore.peptideInfo.protIdSet)).append(","+candiScore.protScore);
            }
            str.append("\n");
            oriWriter.write(str.toString());
        }
        oriWriter.close();
        System.out.println("naive TDA original PSM number: " + numQ001);


        System.out.println("end ori start new");
//        for (ScanRes scanRes : scanResList) {
//            Collections.sort(scanRes.peptideInfoScoreList, Comparator.comparing(o -> o.protScore, Comparator.reverseOrder())); // rank candidates using peptide score
//        }
        for (ScanRes scanRes : scanResList) {
            try {
                Collections.sort(scanRes.peptideInfoScoreList, Comparator.reverseOrder()); // rank candidates using peptide score
            } catch (Exception e){
                System.out.println("lsz error +" +scanRes.scanNum);
            }
//            Collections.sort(scanRes.peptideInfoScoreList, Comparator.reverseOrder()); // rank candidates using peptide score
        }

        Collections.sort(scanResList, Comparator.comparing(o -> (o.peptideInfoScoreList.get(0).pepScore)*(o.peptideInfoScoreList.get(0).protScore+1), Comparator.reverseOrder())); // should still use peptideScore to do FDR
        //calculate new FDR
        fdrList = new ArrayList<>(scanResList.size());
        numTPlusD = 0;
        numD = 0;
        for (ScanRes scanRes : scanResList) {
            CandiScore candiScore = scanRes.peptideInfoScoreList.get(0);
            numTPlusD++;
            if (candiScore.peptideInfo.isDecoy) numD++;
            fdrList.add(2*(double)numD/numTPlusD);
//            System.out.println(scanRes.scanNum + "," + candiScore.protScore + "," + 2*(double)numD/numTPlusD + "," +(candiScore.peptideInfo.isTarget ? 1 : 0));
        }
        numQ001 = 0;
        minQ = 1.0;
        found = false;
        for (int i = fdrList.size()-1; i >= 0; i--) {
            minQ = Math.min(fdrList.get(i), minQ);
            scanResList.get(i).qValue = minQ;
            if (!found && minQ < 0.01) {
                numQ001 = i;
                found = true;
//                break;
            }
        }
        BufferedWriter newWriter = new BufferedWriter(new FileWriter(outputDir+"."+hostName+".Proteins.csv"));
        newWriter.write("scanNum,qValue,TorD,peptide,pepScore,proteins,protscore,varPtmScore,peptide,pepScore,proteins,protscore,varPtmScore,peptide,pepScore,proteins,protscore,varPtmScore,peptide,pepScore,proteins,protscore,varPtmScore,peptide,pepScore,proteins,protscore,varPtmScore,peptide,pepScore,proteins,protscore,varPtmScore,peptide,pepScore,proteins,protscore,varPtmScore,peptide,pepScore,proteins,protscore,varPtmScore,peptide,pepScore,proteins,protscore,varPtmScore,peptide,pepScore,proteins,protscore,varPtmScore,peptide,pepScore,proteins,protscore,varPtmScore,peptide,pepScore,proteins,protscore,varPtmScore,peptide,pepScore,proteins,protscore,varPtmScore\n");

//        Map<Integer, TempRes> scanNumFinalScoreMap = new HashMap<>();
        List<Pair<Double, String>> finalExcelList = new ArrayList<>(scanResList.size());
        for (ScanRes scanRes : scanResList) {
            if (lszDebugScanNum.contains(scanRes.scanNum)) {
                int a = 1;
            }
            List<CandiScore> candiScoreList = scanRes.peptideInfoScoreList;
            CandiScore topCandi = candiScoreList.get(0);
            StringBuilder str = new StringBuilder();
            str.append(scanRes.scanNum+",").append(scanRes.qValue+",").append(topCandi.peptideInfo.isDecoy ? 0 : 1);
            for (CandiScore candiScore : candiScoreList){
                str.append(","+candiScore.ptmContainingSeq).append(","+candiScore.pepScore).append(","+String.join(";", candiScore.peptideInfo.protIdSet)).append(","+candiScore.protScore).append(","+candiScore.varPtmTotalScore);
            }
            str.append("\n");
            newWriter.write(str.toString());

            double theoMass = massTool.calResidueMass(topCandi.ptmContainingSeq) + massTool.H2O;
            double massDiff = getMassDiff(scanRes.expMass, theoMass, MassTool.C13_DIFF);
            double ppm = Math.abs(massDiff * 1e6 / theoMass);
            //TempRes(double pepScore, double protScore, double qValue, boolean isDecoy, String ptmPepSeq)
//            scanNumFinalScoreMap.put(scanRes.scanNum, new TempRes(topCandi.pepScore, topCandi.protScore, scanRes.qValue, !topCandi.peptideInfo.isTarget, topCandi.peptideInfo.seq, 0)); //isdecoy
            double finalScore = topCandi.pepScore*(topCandi.protScore+1);
            String finalStr = String.format(Locale.US, "%s,%d,%f,%d,%s,%s,%s,%s,%s,%s,%f,%f,%f,%d\n"
                    , scanRes.scanName, scanRes.scanNum, scanRes.qValue, topCandi.peptideInfo.isDecoy ? 0 : 1, df.format(finalScore), topCandi.ptmContainingSeq, topCandi.peptideInfo.freeSeq,df.format(topCandi.pepScore)
                    , String.join(";",topCandi.peptideInfo.protIdSet), df.format(topCandi.protScore), ppm, theoMass, scanRes.expMass, scanRes.charge
            );

            finalExcelList.add(new Pair(finalScore, finalStr));
        }
        newWriter.close();
        System.out.println("naive TDA PFM PSM number: " + numQ001);

        // official output with pfm
        Collections.sort(finalExcelList, Comparator.comparing(o -> o.getFirst(), Comparator.reverseOrder()));
        BufferedWriter writer = new BufferedWriter(new FileWriter(outputDir+"PIPI3."+hostName+".csv"));
//        BufferedWriter writer = new BufferedWriter(new FileWriter(outputDir+"pipi7.csv"));

        writer.write("scanName,scanNum,qValue,TorD,finalScore,peptide,freeSeq,pepScore,proteins,protScore,ppm,theoMass,expMass,charge\n");
        for (Pair<Double, String> pair : finalExcelList) {
            writer.write(pair.getSecond());
        }
        writer.close();
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
