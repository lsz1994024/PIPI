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
import ProteomicsLibrary.Types.SparseBooleanVector;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
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
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.locks.ReentrantLock;

public class PIPI {
    private static final Logger logger = LoggerFactory.getLogger(PIPI.class);
    public static final String versionStr = "1.4.7";
    static final boolean useXcorr = false;
    static final int minTagLenToReduceProtDb = 5;

    //// normal
    public static final boolean isPtmSimuTest = false; //normal //todo
    static final boolean usePfmAndReduceDb = true;  //normal //todo
    static final int minTagLenToExtract = 4;  //normal //todo
    static final int maxTagLenToExtract = 99;  //normal //todo
    static final boolean nTermSpecific = false; //normal //todo
    ///// ptmTest
//    public static final boolean isPtmSimuTest = true; //simulation test //todo
//    static final boolean usePfmAndReduceDb = false;  //simulation test //todo
//    static final int minTagLenToExtract = 4;  //simulation test //todo
//    static final int maxTagLenToExtract = 4;  //simulation test //todo
//    static final boolean nTermSpecific = true; //simulation test //todo
    /////


    public static final int[] debugScanNumArray = new int[]{};
    public static final ArrayList<Integer> lszDebugScanNum = new ArrayList<>(Arrays.asList(55966));//35581
    public static void main(String[] args) {
        long startTime = System.nanoTime();

        // Process inputs
        if (args.length != 3) {
            help();
        }

        // Set parameters
        String parameterPath = args[0].trim();
        String spectraPath = args[1].trim();
        String outputDir = args[2].trim();

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
            new PIPI(parameterPath, spectraPath, dbName, hostName, outputDir);



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

    private PIPI(String parameterPath, String spectraPath, String dbName, String hostName, String outputDir) throws Exception {
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
//        String outputDir = parameterMap.get("output_dir");

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
            datasetReader = new DatasetReader(spectraParserArray, ms1Tolerance, ms1ToleranceUnit, massTool, ext, msLevelSet, sqlPath, fileIdNameMap);
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
            datasetReader = new DatasetReader(spectraParserArray, ms1Tolerance, ms1ToleranceUnit, massTool, ext, msLevelSet, sqlPath, fileIdNameMap);
        }
        BufferedReader parameterReader = new BufferedReader(new FileReader("/home/slaiad/Code/PIPI/src/main/resources/ChickOpenTruth.txt"));

        Map<Integer, String> pepTruth = new HashMap<>();
        String line;
        while ((line = parameterReader.readLine()) != null) {
            if (line.isEmpty()) continue;
            line = line.trim();
            String[] splitRes = line.split(",");
            String pepWithMod = splitRes[1].replace('I','L');
            pepTruth.put(Integer.valueOf(splitRes[0]), pepWithMod);
        }

        SpecProcessor specProcessor = new SpecProcessor(massTool);
        Map<String, Integer> precursorChargeMap = new HashMap<>();
        Map<String, Double> precursorMassMap = new HashMap<>();
        TreeMap<Double, Set<String>> pcMassScanNameMap = new TreeMap<>();
        Map<Integer, TreeMap<Double, Set<String>>> fileId_pcMassScanNameMap = new HashMap<>();
        for (int fileId : fileIdNameMap.keySet()) {
            fileId_pcMassScanNameMap.put(fileId, new TreeMap<>());
        }

        logger.info("Get long tags to reduce proteins...");
        int threadNum_0 = Integer.valueOf(parameterMap.get("thread_num"));
        if (threadNum_0 == 0) {
            threadNum_0 = 3 + Runtime.getRuntime().availableProcessors();
        }
        if (java.lang.management.ManagementFactory.getRuntimeMXBean().getInputArguments().toString().indexOf("jdwp") >= 0){
            //change thread 1
//            threadNum_0 = 1;
        }
        System.out.println("thread NUM "+ threadNum_0);

        Map<String, String> mgfTitleMap = new HashMap<>();
        Map<String, Integer> isotopeCorrectionNumMap = new HashMap<>();
        Map<String, Double> ms1PearsonCorrelationCoefficientMap = new HashMap<>();

        //////////==================================================
        //   Get Long Tags
        ExecutorService threadPoolGetLongTag = Executors.newFixedThreadPool(threadNum_0);
        ArrayList<Future<GetLongTag.Entry>> taskListGetLongTag = new ArrayList<>(datasetReader.getUsefulSpectraNum() + 10);
        Connection sqlConSpecCoderX = DriverManager.getConnection(sqlPath);
        Statement sqlStateGetLongTag = sqlConSpecCoderX.createStatement();
        ResultSet sqlResSetGetLongTag = sqlStateGetLongTag.executeQuery("SELECT scanName, scanNum, precursorCharge" +
                ", precursorMass, precursorScanNo, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient FROM spectraTable");

        ReentrantLock lockGetLongTag = new ReentrantLock();

        int submitNumSpecCoderX = 0;
        Set<String> validScanSet = new HashSet<>();

        while (sqlResSetGetLongTag.next()) {
            String scanName = sqlResSetGetLongTag.getString("scanName");
            int scanNum = sqlResSetGetLongTag.getInt("scanNum");
            int precursorCharge = sqlResSetGetLongTag.getInt("precursorCharge");
            double precursorMass = sqlResSetGetLongTag.getDouble("precursorMass");
            mgfTitleMap.put(scanName, sqlResSetGetLongTag.getString("mgfTitle"));
            isotopeCorrectionNumMap.put(scanName, sqlResSetGetLongTag.getInt("isotopeCorrectionNum"));
            ms1PearsonCorrelationCoefficientMap.put(scanName, sqlResSetGetLongTag.getDouble("ms1PearsonCorrelationCoefficient"));

            boolean shouldRun = false;
            for (int debugScanNum : lszDebugScanNum) {
                if (Math.abs(scanNum-debugScanNum) < 5) {
                    shouldRun = true;
                }
            }
            if (java.lang.management.ManagementFactory.getRuntimeMXBean().getInputArguments().toString().indexOf("jdwp") >= 0){
                if (!shouldRun) {
                    continue;//22459   comment this continue line, if run all scans
                }
            }

            int fileId = fileNameIdMap.get( scanName.split("\\.")[0] );

            precursorChargeMap.put(scanName, precursorCharge);
            precursorMassMap.put(scanName, precursorMass);
            submitNumSpecCoderX++;
            validScanSet.add(scanName);
            taskListGetLongTag.add(threadPoolGetLongTag.submit(new GetLongTag(scanNum, buildIndex, massTool, spectraParserArray[fileId], minClear, maxClear, lockGetLongTag, scanName, precursorCharge
                    , precursorMass, specProcessor, pepTruth.get(1886), ms2Tolerance)));
        }
        System.out.println("totalSubmit in SpecCoder, "+ submitNumSpecCoderX);
        sqlResSetGetLongTag.close();
        sqlStateGetLongTag.close();

        int lastProgressGetLongTag = 0;
        int totalCountGetLongTag = taskListGetLongTag.size();
        int countGetLongTag = 0;
        Map<String, Map<String, Double>> protTagScoreMapMap = new HashMap<>();
        while (countGetLongTag < totalCountGetLongTag) {
            // record search results and delete finished ones.
            List<Future<GetLongTag.Entry>> toBeDeleteTaskList = new ArrayList<>(totalCountGetLongTag - countGetLongTag);
            for (Future<GetLongTag.Entry> task : taskListGetLongTag) {
                if (task.isDone()) {
                    if (task.get() != null ) {
                        GetLongTag.Entry entry = task.get();
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


        //////////==================================================
        //   Score for Proteins
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
        Collections.sort(protScoreLongList, Comparator.comparing(o -> o.getSecond(), Comparator.reverseOrder()));

        Set<String> reducedProtIdSet = new HashSet<>();
        int ii = 0;
        for (Pair<String, Double> pair : protScoreLongList){
            ii++;
            if (ii > 35000) break;
            reducedProtIdSet.add(pair.getFirst());
        }

        if (! usePfmAndReduceDb)
        {
            reducedProtIdSet = protLengthMap.keySet();  //dont reduce for simulation dataset
        }

        if (java.lang.management.ManagementFactory.getRuntimeMXBean().getInputArguments().toString().indexOf("jdwp") >= 0) {
            reducedProtIdSet = protLengthMap.keySet();  //dont reduce for simulation dataset
        }
        // reduce proteins from buildIndex if it is not contained in reducedProtIdSet
        Iterator<String> iter = buildIndex.protSeqMap.keySet().iterator();
        while (iter.hasNext()) {
            if (!reducedProtIdSet.contains(iter.next())) {
                iter.remove();
            }
        }

        System.out.println("Protein size: "+protLengthMap.size()+" reduced to " + reducedProtIdSet.size());


        //  END Score for Proteins
        //////////==================================================


        //////////==================================================
        //   Get Tags and Get Peptide Candidates
        logger.info("Building decoy...");
        int threadNum_1 = Integer.valueOf(parameterMap.get("thread_num"));
        if (threadNum_1 == 0) {
            threadNum_1 = 3 + Runtime.getRuntime().availableProcessors();
        }
        if (java.lang.management.ManagementFactory.getRuntimeMXBean().getInputArguments().toString().indexOf("jdwp") >= 0){
            //change thread 1
//            threadNum_1 = 1;
        }
        System.out.println("thread NUM "+ threadNum_1);

        ExecutorService threadPoolBuildDecoyProts = Executors.newFixedThreadPool(threadNum_1);
        ArrayList<Future<BuildDecoyProts.Entry>> taskListBuildDecoyProts = new ArrayList<>(reducedProtIdSet.size() + 10);
        ReentrantLock lockSpecCoder = new ReentrantLock();

//        int submitNumGetPepCandis = 0;
        for (String protId : reducedProtIdSet) {    //All protSeq are recorded but some decoy version are not because their targets are not in the reducedProtIdSet
            taskListBuildDecoyProts.add(threadPoolBuildDecoyProts.submit(new BuildDecoyProts(parameterMap, buildIndex, protId)));
        }
//        System.out.println("totalSubmit in GetPepCandi, "+ submitNumGetPepCandis);

        int lastProgressBuildDecoyProts = 0;
        int totalCountBuildDecoyProts = taskListBuildDecoyProts.size();
        int countBuildDecoyProts = 0;
//        Multimap<String, String> pepProtsMap = HashMultimap.create();
//        Map<String, Double> pepMassMap = new HashMap<>(500000);
//        Map<String, SparseBooleanVector> pepCodeMap = new HashMap<>();
        double minPeptideMass = 9999;
        double maxPeptideMass = 0;

//        Set<String> targetPepSet = new HashSet<>();
//        Set<String> decoyPepSet = new HashSet<>();
        Map<String, Set<Pair<String, Integer>>> tagProtPosMap = new HashMap<>();
        while (countBuildDecoyProts < totalCountBuildDecoyProts) {
            // record search results and delete finished ones.
            List<Future<BuildDecoyProts.Entry>> toBeDeleteTaskList = new ArrayList<>(totalCountBuildDecoyProts - countBuildDecoyProts);
            for (Future<BuildDecoyProts.Entry> task : taskListBuildDecoyProts) {
                if (task.isDone()) {
                    if (task.get() != null ) {
                        BuildDecoyProts.Entry entry = task.get();
                        String protId = entry.protId;
                        String decoyProtSeq = entry.decoyProtSeq;
                        String decoyProtId = "DECOY_"+protId;
                        buildIndex.protSeqMap.put(decoyProtId, decoyProtSeq);

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

        // build reduced protein fm index
//        long t1=System.currentTimeMillis();
        logger.info("Generating reduced FM index...");

        BufferedWriter writerProt = new BufferedWriter(new FileWriter("catProtReduced.txt"));
        int dotPos = 0;
        int dotNum = 0;
        buildIndex.dotPosArrReduced = new int[buildIndex.protSeqMap.keySet().size()];

        for (String protId : buildIndex.protSeqMap.keySet()) {
            buildIndex.dotPosArrReduced[dotNum] = dotPos;
            buildIndex.posProtMapReduced.put(dotNum, protId);
            String protSeq = buildIndex.protSeqMap.get(protId).replace('I', 'L');
            buildIndex.protSeqMap.put(protId, protSeq);
            writerProt.write("." + protSeq.replace('I', 'L'));
            dotNum++;
            dotPos += protSeq.length()+1;
            int numOfTags = buildIndex.inferSegment.getLongTagNumForProt(protSeq);
            protLengthMap.put(protId, numOfTags);
        }
        writerProt.close();
//        long t2=System.currentTimeMillis();
        char[] text = buildIndex.loadFile("catProtReduced.txt", true);
//        long t3=System.currentTimeMillis();
        buildIndex.fmIndexReduced = new FMIndex(text);
//        long t4=System.currentTimeMillis();
//
//        System.out.println("time," + (t4-t3) + "," + (t3-t2)+ "," + (t2-t1));
        buildIndex.minPeptideMass = minPeptideMass;
        buildIndex.maxPeptideMass = maxPeptideMass;
        buildIndex.tagProtPosMap = tagProtPosMap;
        // writer concatenated fasta
        Map<String, String> proteinAnnotationMap;
        String dbPath = parameterMap.get("db");
        DbTool dbTool = buildIndex.dbTool;
        DbTool contaminantsDb = null;
        if (true) { //addContaminants = true
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

        logger.info("Pre searching...");
        int threadNum = Integer.valueOf(parameterMap.get("thread_num"));
        if (threadNum == 0) {
            threadNum = 3 + Runtime.getRuntime().availableProcessors();
        }
        if (java.lang.management.ManagementFactory.getRuntimeMXBean().getInputArguments().toString().indexOf("jdwp") >= 0){
            // change thread PreSearch
            threadNum = 1;
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
            int scanNum = Integer.valueOf(scanNameStr[2]);
            boolean shouldRun = false;
            for (int debugScanNum : lszDebugScanNum) {
                if (Math.abs(scanNum-debugScanNum) < 5) {
                    shouldRun = true;
                }
            }
            if (java.lang.management.ManagementFactory.getRuntimeMXBean().getInputArguments().toString().indexOf("jdwp") >= 0){
                if (!shouldRun) {
                    continue;//22459   comment this continue line, if run all scans
                }
            }
            int fileId = fileNameIdMap.get( scanNameStr[0] );
            taskListBone.add(threadPoolBone.submit(new PreSearch(scanNum, buildIndex, massTool, ms2Tolerance, ms1Tolerance, leftInverseMs1Tolerance, rightInverseMs1Tolerance
                    , ms1ToleranceUnit, inferPTM.getMinPtmMass(), inferPTM.getMaxPtmMass(), Math.min(precursorCharge > 1 ? precursorCharge - 1 : 1, 3)
                    , spectraParserArray[fileId], minClear, maxClear, lockBone, scanName, precursorCharge, precursorMass, specProcessor ,pepTruth.get(1886))));
        }
        System.out.println("totalSubmit in Bone, "+ submitNumBone);

        Map<String, PeptideInfo> allPeptideInfoMap = new HashMap<>();
        int lastProgressBone = 0;
        int resultCountBone = 0;
        int totalCountBone = taskListBone.size();
        int countBone = 0;

        Connection sqlConSpecCoder = DriverManager.getConnection(sqlPath);
        PreparedStatement sqlPreparedStatement = sqlConSpecCoder.prepareStatement("REPLACE INTO spectraTable (scanNum, scanName,  precursorCharge, precursorMass, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient, labelling, peptide, theoMass, isDecoy, globalRank, normalizedCorrelationCoefficient, score, deltaLCn, deltaCn, matchedPeakNum, ionFrac, matchedHighestIntensityFrac, explainedAaFrac, otherPtmPatterns, aScore, candidates, peptideSet, whereIsTopCand, shouldPtm, hasPTM, ptmNum, isSettled) VALUES (?, ?, ?, ?,?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)");
        sqlConSpecCoder.setAutoCommit(false);
        while (countBone < totalCountBone) {
            // record search results and delete finished ones.
            List<Future<PreSearch.Entry>> toBeDeleteTaskList = new ArrayList<>(totalCountBone - countBone);
            for (Future<PreSearch.Entry> task : taskListBone) {
                if (task.isDone()) {
                    if (task.get() != null) {
                        PreSearch.Entry entry = task.get();

                        allPeptideInfoMap.putAll(entry.peptideInfoMapForRef);
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
                        sqlPreparedStatement.setInt(25, -1);//whereIsTopCand dummy to be deleted
                        sqlPreparedStatement.setInt(26, entry.shouldPtm);
                        sqlPreparedStatement.setInt(27, entry.hasPTM);
                        sqlPreparedStatement.setInt(28, entry.ptmNum);
                        sqlPreparedStatement.setInt(29, entry.isSettled);

                        sqlPreparedStatement.executeUpdate();

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
        sqlConSpecCoder.commit();
        sqlConSpecCoder.setAutoCommit(true);
        sqlConSpecCoder.close();
        if (lockBone.isLocked()) {
            lockBone.unlock();
        }

        System.out.println("resultCount," +resultCountBone);

        System.out.println("lsz +" +","+allPeptideInfoMap.size());
//        logger.info("Start searching...");
//        if (threadNum == 0) {
//            threadNum = 3 + Runtime.getRuntime().availableProcessors();
//        }
//        if (java.lang.management.ManagementFactory.getRuntimeMXBean().getInputArguments().toString().indexOf("jdwp") >= 0){
////            change thread ptmSearch
//            threadNum = 1;
//        }
//        ExecutorService threadPoolPtm = Executors.newFixedThreadPool(threadNum);
//        ArrayList<Future<PtmSearch.Entry>> taskListPTM = new ArrayList<>(scanNamePeptideInfoMap.keySet().size() + 10);
//        ReentrantLock lockPtm = new ReentrantLock();
//        Binomial binomial = new Binomial(Integer.valueOf(parameterMap.get("max_peptide_length")) * 2);
//        int submitTimePtm = 0;
//
//        Map<String, PeptideInfo> allPeptideInfoMap = new HashMap<>();
//        Map<String, PeptideInfo> thisScanPeptideInfoMap;
//        Map<String, PeptideInfo> otherScanPeptideInfoMap;
//        for (int thisFileId : fileIdNameMap.keySet()){
//            TreeMap<Double, Set<String>> local_this_pcMassScanNameMap = fileId_pcMassScanNameMap.get(thisFileId);
//            for (double mass : local_this_pcMassScanNameMap.keySet()) {
//                for (String thisScanName : local_this_pcMassScanNameMap.get(mass)){
//                    thisScanPeptideInfoMap = scanNamePeptideInfoMap.get(thisScanName);
//                    int thisScanNum = Integer.valueOf( thisScanName.split("\\.")[2] );
//                    if (lszDebugScanNum.contains(thisScanNum)) {
//                        int a = 1;
//                    }
////                    int thisFileId = Integer.valueOf( thisScanName.split("\\.")[0] );
//                    Set<Peptide> realPtmOnlyList = new HashSet<>();
//                    Set<Peptide> realPtmFreeList = new HashSet<>();
//                    Map<String, PeptideInfo> realPeptideInfoMap = new HashMap<>();
//                    for (Peptide pep : ptmOnlyCandiMap.get(thisScanName)) {
//                        realPtmOnlyList.add(pep.clone());
//                        realPeptideInfoMap.put(pep.getFreeSeq(), thisScanPeptideInfoMap.get(pep.getFreeSeq()).clone());
//                        allPeptideInfoMap.put(pep.getFreeSeq(), thisScanPeptideInfoMap.get(pep.getFreeSeq()).clone());
//                    }
//                    for (Peptide pep : ptmFreeCandiMap.get(thisScanName)) {
//                        realPtmFreeList.add(pep.clone());
//                        realPeptideInfoMap.put(pep.getFreeSeq(), thisScanPeptideInfoMap.get(pep.getFreeSeq()).clone());
//                        allPeptideInfoMap.put(pep.getFreeSeq(), thisScanPeptideInfoMap.get(pep.getFreeSeq()).clone());
//                    }
//
//                    TreeMap<Double, Set<String>> local_other_pcMassScanNameMap = fileId_pcMassScanNameMap.get(thisFileId);
//                    for (int i = -3; i <= 3; i++) {
//                        for (Set<String> otherScanNameSet : local_other_pcMassScanNameMap.subMap(mass + i*MassTool.PROTON - 0.02, true, mass + i*MassTool.PROTON + 0.02, true).values()) {
//                            for (String otherScanName : otherScanNameSet) {
//                                int otherScanNum = Integer.valueOf( otherScanName.split("\\.")[2] );
////                                int otherFileId = Integer.valueOf( otherScanName.split("\\.")[0] );
////                                if (otherFileId != thisFileId) continue;
//                                otherScanPeptideInfoMap = scanNamePeptideInfoMap.get(otherScanName);
//                                if (otherScanNum < thisScanNum+2000 && otherScanNum > thisScanNum-2000 && otherScanNum != thisScanNum) {
//                                    for (Peptide pep : ptmOnlyCandiMap.get(otherScanName)) {
//                                        if (Math.abs(pep.getTheoMass()-precursorMassMap.get(thisScanName)) < 0.02) {
//                                            realPtmFreeList.add(pep.clone());
//                                        } else {
//                                            realPtmOnlyList.add(pep.clone());
//                                        }
//                                        realPeptideInfoMap.put(pep.getFreeSeq(), otherScanPeptideInfoMap.get(pep.getFreeSeq()).clone());
//                                    }
//                                    for (Peptide pep : ptmFreeCandiMap.get(otherScanName)) {
//                                        if (Math.abs(pep.getTheoMass()-precursorMassMap.get(thisScanName)) < 0.02) {
//                                            realPtmFreeList.add(pep.clone());
//                                        } else {
//                                            realPtmOnlyList.add(pep.clone());
//                                        }
//                                        realPeptideInfoMap.put(pep.getFreeSeq(), otherScanPeptideInfoMap.get(pep.getFreeSeq()).clone());
//                                    }
//                                }
//                            }
//                        }
//                    }
//
//                    int precursorCharge = precursorChargeMap.get(thisScanName);
//                    int precursorScanNo = 0;
//                    double precursorMass = precursorMassMap.get(thisScanName);
//                    submitTimePtm++;
////                    System.out.println(thisScanName+" ptmSearch:"+realPtmOnlyList.size()+","+realPtmFreeList.size()+","+realPeptideInfoMap.size());
//                    taskListPTM.add(threadPoolPtm.submit(new PtmSearch(thisScanNum, buildIndex, massTool, ms1Tolerance
//                            , leftInverseMs1Tolerance, rightInverseMs1Tolerance, ms1ToleranceUnit, ms2Tolerance
//                            , inferPTM.getMinPtmMass(), inferPTM.getMaxPtmMass(), Math.min(precursorCharge > 1 ? precursorCharge - 1 : 1, 3)
//                            , spectraParserArray[thisFileId], minClear, maxClear, lockPtm, thisScanName, precursorCharge
//                            , precursorMass, inferPTM, specProcessor, sqlPath, binomial, precursorScanNo, realPtmOnlyList, realPtmFreeList, realPeptideInfoMap)));
//                }
//            }
//        }
//        System.out.println("submit times in Ptm" + submitTimePtm);
//        // check progress every minute, record results,and delete finished tasks.
//        int lastProgressPtm = 0;
//        int resultCountPtm = 0;
//        int totalCountPtm = taskListPTM.size();
//        int countPtm = 0;
//
//        Map<String, List<ExRes>> exResForScanNameMap = new HashMap<>();
////        System.out.println("lsz totalCount, " + totalCount);
//        while (countPtm < totalCountPtm) {
//            // record search results and delete finished ones.
//            List<Future<PtmSearch.Entry>> toBeDeleteTaskList = new ArrayList<>(totalCountPtm - countPtm);
//            for (Future<PtmSearch.Entry> task : taskListPTM) {
//                if (task.isDone()) {
//                    if (task.get() != null) {
//                        PtmSearch.Entry entry = task.get();
//                        sqlPreparedStatement.setInt(1, entry.scanNum);
//                        sqlPreparedStatement.setString(2, entry.scanName);
//                        sqlPreparedStatement.setInt(3, entry.precursorCharge);
//                        sqlPreparedStatement.setDouble(4, entry.precursorMass);
//                        sqlPreparedStatement.setString(5, mgfTitleMap.get(entry.scanName));
//                        sqlPreparedStatement.setInt(6, isotopeCorrectionNumMap.get(entry.scanName));
//                        sqlPreparedStatement.setDouble(7, ms1PearsonCorrelationCoefficientMap.get(entry.scanName));
//                        sqlPreparedStatement.setString(8, entry.labelling);
//                        sqlPreparedStatement.setString(9, entry.peptide);
//                        sqlPreparedStatement.setDouble(10, entry.theoMass);
//                        sqlPreparedStatement.setInt(11, entry.isDecoy);
//                        sqlPreparedStatement.setInt(12, entry.globalRank);
//                        sqlPreparedStatement.setDouble(13, entry.normalizedCorrelationCoefficient);
//                        sqlPreparedStatement.setDouble(14, entry.score);
//                        sqlPreparedStatement.setDouble(15, entry.deltaLCn);
//                        sqlPreparedStatement.setDouble(16, entry.deltaCn);
//                        sqlPreparedStatement.setInt(17, entry.matchedPeakNum);
//                        sqlPreparedStatement.setDouble(18, entry.ionFrac);
//                        sqlPreparedStatement.setDouble(19, entry.matchedHighestIntensityFrac);
//                        sqlPreparedStatement.setDouble(20, entry.explainedAaFrac);
//                        sqlPreparedStatement.setString(21, entry.otherPtmPatterns);
//                        sqlPreparedStatement.setString(22, entry.aScore);
//                        sqlPreparedStatement.setString(23, entry.candidates);
//                        sqlPreparedStatement.setString(24, entry.peptideSet);
//                        sqlPreparedStatement.setInt(25, entry.whereIsTopCand);
//                        sqlPreparedStatement.setInt(26, entry.shouldPtm);
//                        sqlPreparedStatement.setInt(27, entry.hasPTM);
//                        sqlPreparedStatement.setInt(28, entry.ptmNum);
//                        sqlPreparedStatement.setInt(29, entry.isSettled);
//
//                        sqlPreparedStatement.executeUpdate();
//                        ++resultCountPtm;
//                    }
//
//                    toBeDeleteTaskList.add(task);
//                } else if (task.isCancelled()) {
//                    toBeDeleteTaskList.add(task);
//                }
//            }
//            countPtm += toBeDeleteTaskList.size();
//            taskListPTM.removeAll(toBeDeleteTaskList);
//            taskListPTM.trimToSize();
//
//            sqlConSpecCoder.commit();
//
//            int progress = countPtm * 20 / totalCountPtm;
//            if (progress != lastProgressPtm) {
//                logger.info("Searching {}%...", progress * 5);
////                System.out.println(toBeDeleteTaskList.size()+ ","+ taskListPTM.size() +"," + count2 + ", " + totalCount + "," + lastProgress);
//                lastProgressPtm = progress;
//            }
//
//            if (countPtm == totalCountPtm) {
//                break;
//            }
//            Thread.sleep(6000);
//        }
//        System.out.println("exResForScanNameMap, "+exResForScanNameMap.size());
//        // shutdown threads.
//        threadPoolPtm.shutdown();
//        if (!threadPoolPtm.awaitTermination(60, TimeUnit.SECONDS)) {
//            threadPoolPtm.shutdownNow();
//            if (!threadPoolPtm.awaitTermination(60, TimeUnit.SECONDS))
//                throw new Exception("Pool did not terminate");
//        }
//
//        sqlConSpecCoder.commit();
//        sqlConSpecCoder.setAutoCommit(true);
//        sqlConSpecCoder.close();
//        if (lockPtm.isLocked()) {
//            lockPtm.unlock();
//        }
//
//        if (resultCountPtm == 0) {
//            throw new Exception("There is no useful results.");
//        }

        String percolatorInputFileName = spectraPath + "." + labelling + ".input.temp";
        pfm(outputDir, allPeptideInfoMap, sqlPath, buildIndex.protSeqMap, massTool, hostName);

        writePercolator(percolatorInputFileName, allPeptideInfoMap, sqlPath);
        Map<String, PercolatorEntry> percolatorResultMap = null;

        if (parameterMap.get("add_decoy").contentEquals("0")) {
            logger.warn("add_decoy = 0. Don't estimate FDR.");
        } else {
            logger.info("Estimating FDR...");
            String percolatorOutputFileName = spectraPath + "." + labelling + ".output.temp";
            String percolatorDecoyOutputFileName = spectraPath + "." + labelling + ".Decoyoutput.temp";
            String percolatorProteinOutputFileName = spectraPath + "." + labelling + ".protein.tsv";
            percolatorResultMap = runPercolator(percolatorPath, percolatorInputFileName, percolatorOutputFileName, percolatorDecoyOutputFileName, percolatorProteinOutputFileName, parameterMap.get("db") + ".TD.fasta", parameterMap.get("enzyme_name_1"));
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
//        writeDebugResults(sqlPath, spectraPath);

        writeFinalResult(percolatorResultMap, outputDir + "." + hostName + ".pipi.csv", allPeptideInfoMap, sqlPath, outputDir, hostName);
//        new WritePepXml(spectraPath + "." + labelling + ".pipi.pep.xml", spectraPath, parameterMap, massTool.getMassTable(), percolatorResultMap, buildIndex.getPepInfoMap(), buildIndex.returnFixModMap(), sqlPath);
    }

    private boolean isTarget(Collection<String> proteinIds) {
        for (String protein : proteinIds) {
            if (!protein.startsWith("DECOY_")) {
                return true;
            }
        }
        return false;
    }
    class ExRes implements Cloneable{
        public ExRes clone() throws CloneNotSupportedException {
            super.clone();
            ExRes other = new ExRes();
            other.scanNum = scanNum;
            other.scanName = scanName;
            other.precursorCharge = precursorCharge;
            other.precursorMass = precursorMass;
            other.mgfTitle = mgfTitle;
            other.isotopeCorrectionNum = isotopeCorrectionNum;
            other.ms1PearsonCorrelationCoefficient = ms1PearsonCorrelationCoefficient;
            other.labelling = labelling;
            other.peptide = peptide;
            other.theoMass = theoMass;
            other.isDecoy = isDecoy;
            other.globalRank = globalRank;
            other.normalizedCorrelationCoefficient = normalizedCorrelationCoefficient;
            other.score = score;
            other.deltaLCn = deltaLCn;
            other.deltaCn = deltaCn;
            other.matchedPeakNum = matchedPeakNum;
            other.ionFrac = ionFrac;
            other.matchedHighestIntensityFrac = matchedHighestIntensityFrac;
            other.explainedAaFrac = explainedAaFrac;
            other.otherPtmPatterns = otherPtmPatterns;
            other.aScore = aScore;
            other.candidates = candidates;
            other.peptideSet = peptideSet;
            other.whereIsTopCand = whereIsTopCand;
            other.shouldPtm = shouldPtm;
            other.hasPTM = hasPTM;
            other.ptmNum = ptmNum;
            other.isSettled = isSettled;

            return other;
        }
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
        public double getScore() {return score;}
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

    private boolean isHomo(ExRes r1, ExRes r2, Map<String, PepInfo> peptide0Map) {
        if (!peptide0Map.containsKey(r1.peptide.replaceAll("\\(-?(\\d+)(\\.\\d+)\\)","")) || !peptide0Map.containsKey(r2.peptide.replaceAll("\\(-?(\\d+)(\\.\\d+)\\)",""))) {
            int  cc = 1;
        }
        PepInfo pep01 = peptide0Map.get(r1.peptide.replaceAll("\\(-?(\\d+)(\\.\\d+)\\)",""));
        PepInfo pep02 = peptide0Map.get(r2.peptide.replaceAll("\\(-?(\\d+)(\\.\\d+)\\)",""));

        String[] prots1 = pep01.proteins;
        String[] prots2 = pep02.proteins;

        HashSet<String> set = new HashSet<>(Arrays.asList(prots1));
        set.retainAll(Arrays.asList(prots2));
        if (set.isEmpty()) return false;
//        pep01.code.dot(pep02.code) > 0.5*pep01.code.union(pep02.code)
        return pep01.code.dot(pep02.code) > 0.3*Math.min(r1.peptide.replaceAll("\\(-?(\\d+)(\\.\\d+)\\)","").length(), r2.peptide.replaceAll("\\(-?(\\d+)(\\.\\d+)\\)","").length());
    }

    private boolean isHomoTest(ExRes r1, ExRes r2, Map<String, PepInfo> peptide0Map) {
        PepInfo pep01 = peptide0Map.get(r1.peptide.replaceAll("\\(-?(\\d+)(\\.\\d+)\\)",""));
        PepInfo pep02 = peptide0Map.get(r2.peptide.replaceAll("\\(-?(\\d+)(\\.\\d+)\\)",""));

        return pep01.code.dot(pep02.code) > 0.2*pep01.code.union(pep02.code);
    }
    private void writePercolator(String resultPath, Map<String, PeptideInfo> allPeptideInfoMap, String sqlPath) throws IOException, SQLException , CloneNotSupportedException{
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

                PeptideInfo pepInfo = allPeptideInfoMap.get(peptide.replaceAll("[^A-Z]+", ""));
                TreeSet<String> proteinIdSet = new TreeSet<>();
                for (String protein : pepInfo.protIdSet) {
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
                int isDecoy = pepInfo.isDecoy ? 1 : 0;
                int globalRank = sqlResultSet.getInt("globalRank");
                double normalizedCorrelationCoefficient = sqlResultSet.getDouble("normalizedCorrelationCoefficient");
                double score = sqlResultSet.getDouble("score");
                double deltaLCn = sqlResultSet.getDouble("deltaLCn");
                double deltaCn = sqlResultSet.getDouble("deltaCn");
                double ionFrac = sqlResultSet.getDouble("ionFrac");
                double matchedHighestIntensityFrac = sqlResultSet.getDouble("matchedHighestIntensityFrac");
                double explainedAaFrac = sqlResultSet.getDouble("explainedAaFrac");

//                if (!hasCo) {
                if (isDecoy == 1) {
                    writer.write(scanName + "\t-1\t" + scanNum + "\t" + expMass + "\t" + theoMass + "\t" + score + "\t" + deltaCn + "\t" + deltaLCn + "\t" + normalizedCorrelationCoefficient + "\t" + 1 + "\t" + Math.abs(massDiff * 1e6 / theoMass) + "\t" + ionFrac + "\t" + matchedHighestIntensityFrac + "\t" + sb.toString() + explainedAaFrac + "\t" + pepInfo.leftFlank + "." + peptide.replaceAll("\\(", "[").replaceAll("\\)", "]") + "." + pepInfo.rightFlank + "\t" + String.join("\t", proteinIdSet) + "\n"); // Percolator only recognize "[]".
                } else {
                    writer.write(scanName +  "\t1\t" + scanNum + "\t" + expMass + "\t" + theoMass + "\t" + score + "\t" + deltaCn + "\t" + deltaLCn + "\t" + normalizedCorrelationCoefficient + "\t" + 1 + "\t" + Math.abs(massDiff * 1e6 / theoMass) + "\t" + ionFrac + "\t" + matchedHighestIntensityFrac + "\t" + sb.toString() + explainedAaFrac + "\t" + pepInfo.leftFlank + "." + peptide.replaceAll("\\(", "[").replaceAll("\\)", "]") + "." + pepInfo.rightFlank + "\t" + String.join("\t", proteinIdSet) + "\n"); // Percolator only recognize "[]".
                }
//                }
            }
        }
        sqlResultSet.close();
        sqlStatement.close();

        writer.close();
        sqlConnection.close();
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
    class CandiScore{
        public String ptmContainingSeq;
        public PeptideInfo peptideInfo;
        public double pepScore;
        public double protScore = 0;
        CandiScore(PeptideInfo peptideInfo, double pepScore, String ptmContainingSeq) {
            this.peptideInfo = peptideInfo;
            this.pepScore = pepScore;
            this.ptmContainingSeq = ptmContainingSeq;
        }
    }
    private void pfm(String outputDir, Map<String, PeptideInfo> allPeptideInfoMap, String sqlPath, Map<String, String> protSeqMap, MassTool massTool, String hostName) throws IOException, SQLException , CloneNotSupportedException{
        //collect data
        Connection sqlConnection = DriverManager.getConnection(sqlPath);
        Statement sqlStatement = sqlConnection.createStatement();
        ResultSet sqlResultSet = sqlStatement.executeQuery("SELECT scanNum, scanName, precursorCharge, precursorMass, peptide, theoMass, isDecoy, globalRank, normalizedCorrelationCoefficient, score, deltaLCn, deltaCn, matchedPeakNum, ionFrac, matchedHighestIntensityFrac, explainedAaFrac, peptideSet FROM spectraTable");
        Map<String, Map<String, Double>> protPepScoreMap = new HashMap<>();

        List<ScanRes> scanResList = new ArrayList<>();
        while (sqlResultSet.next()) {
            String peptide = sqlResultSet.getString("peptide"); //this is ptmContaining Seq without n or c
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

                PeptideInfo pepInfo = allPeptideInfoMap.get(peptide.replaceAll("[^A-Z]+", ""));
                List<CandiScore> candiScoreList = new ArrayList<>();
                double firstScore = 0;
                for (int i = 0; i < numPep; i++) {
                    //dont put in the impossible ones
                    String ptmContainingSeq = candiSetStr[3*i+0];
                    PeptideInfo candiPeptideInfo = allPeptideInfoMap.get(ptmContainingSeq.replaceAll("[^A-Z]+", ""));
                    double thisScore = Double.valueOf(candiSetStr[3*i+1]);
                    if (i == 0) {
                        firstScore = thisScore;
                    } else { //if any non-top candidates have poor score, break.
                        if (thisScore < firstScore*0.5) {
                            break;
                        }
                    }
                    candiScoreList.add(new CandiScore(candiPeptideInfo, thisScore, ptmContainingSeq)); //peptideInfo and their score
                }
                Collections.sort(candiScoreList, Comparator.comparing(o -> o.pepScore, Comparator.reverseOrder()));
                scanResList.add(new ScanRes(scanName, scanNum, expMass, candiScoreList, charge));

                for (String protId : pepInfo.protIdSet){
                    if (protPepScoreMap.containsKey(protId)){
                        Map<String, Double> pepScoreMap = protPepScoreMap.get(protId);
                        if (pepScoreMap.containsKey(peptide)){
                            pepScoreMap.put(peptide, Math.max(pepScoreMap.get(peptide), score)); // use max score of repeated peptides for any protein
                        } else {
                            pepScoreMap.put(peptide, score);
                        }
                    } else {
                        Map<String, Double> pepScoreMap = new HashMap<>();
                        pepScoreMap.put(peptide, score);
                        protPepScoreMap.put(protId, pepScoreMap);
                    }
                }
//                String scanName = sqlResultSet.getString("scanName");
//                int isDecoy = sqlResultSet.getInt("isDecoy");
//                int globalRank = sqlResultSet.getInt("globalRank");
//                double normalizedCorrelationCoefficient = sqlResultSet.getDouble("normalizedCorrelationCoefficient");
//                double deltaLCn = sqlResultSet.getDouble("deltaLCn");
//                double deltaCn = sqlResultSet.getDouble("deltaCn");
//                double ionFrac = sqlResultSet.getDouble("ionFrac");
//                double matchedHighestIntensityFrac = sqlResultSet.getDouble("matchedHighestIntensityFrac");
//                double explainedAaFrac = sqlResultSet.getDouble("explainedAaFrac");
            }
        }
        sqlResultSet.close();
        sqlStatement.close();
        sqlConnection.close();
        //calculate prot score
        Map<String, Double> protScoreMap = new HashMap<>();
        for (String protId : protPepScoreMap.keySet()) {
            Map<String,Double> pepScoreMap = protPepScoreMap.get(protId);
            for (String pep : pepScoreMap.keySet()){
                if (pepScoreMap.get(pep) > 3.5) {     //this peptide score threshold is empirical
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
//            protScoreMap.put(protId, protScoreMap.get(protId) / massTool.dummyDigest(protSeqMap.get(protId), 0).size());
//            System.out.println(protId + "," + protScoreMap.get(protId));
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
        for (ScanRes scanRes : scanResList) {
            Collections.sort(scanRes.peptideInfoScoreList, Comparator.comparing(o -> o.protScore, Comparator.reverseOrder())); // rank candidates using peptide score
        }

        Collections.sort(scanResList, Comparator.comparing(o -> (o.peptideInfoScoreList.get(0).pepScore)*(o.peptideInfoScoreList.get(0).protScore+1), Comparator.reverseOrder())); // should still use peptideScore to do FDR
        if (!usePfmAndReduceDb){  //for ptm simulation, use dont  do pfm , use pepscore as it was
            Collections.sort(scanResList, Comparator.comparing(o -> o.peptideInfoScoreList.get(0).pepScore, Comparator.reverseOrder()));
        }
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
        newWriter.write("scanNum,qValue,TorD,peptide,pepScore,proteins,protscore,peptide,pepScore,proteins,protscore,peptide,pepScore,proteins,protscore,peptide,pepScore,proteins,protscore,peptide,pepScore,proteins,protscore,peptide,pepScore,proteins,protscore,peptide,pepScore,proteins,protscore,peptide,pepScore,proteins,protscore,peptide,pepScore,proteins,protscore,peptide,pepScore,proteins,protscore,peptide,pepScore,proteins,protscore,peptide,pepScore,proteins,protscore,peptide,pepScore,proteins,protscore\n");

//        Map<Integer, TempRes> scanNumFinalScoreMap = new HashMap<>();
        List<Pair<Double, String>> finalExcelList = new ArrayList<>(scanResList.size());
        for (ScanRes scanRes : scanResList) {
            List<CandiScore> candiScoreList = scanRes.peptideInfoScoreList;
            CandiScore topCandi = candiScoreList.get(0);
            StringBuilder str = new StringBuilder();
            str.append(scanRes.scanNum+",").append(scanRes.qValue+",").append(topCandi.peptideInfo.isDecoy ? 0 : 1);
            for (CandiScore candiScore : candiScoreList){
                str.append(","+candiScore.peptideInfo.freeSeq).append(","+candiScore.pepScore).append(","+String.join(";", candiScore.peptideInfo.protIdSet)).append(","+candiScore.protScore);
            }
            str.append("\n");
            newWriter.write(str.toString());

            double theoMass = massTool.calResidueMass(topCandi.ptmContainingSeq) + massTool.H2O;
            double massDiff = getMassDiff(scanRes.expMass, theoMass, MassTool.C13_DIFF);
            double ppm = Math.abs(massDiff * 1e6 / theoMass);
            //TempRes(double pepScore, double protScore, double qValue, boolean isDecoy, String ptmPepSeq)
//            scanNumFinalScoreMap.put(scanRes.scanNum, new TempRes(topCandi.pepScore, topCandi.protScore, scanRes.qValue, !topCandi.peptideInfo.isTarget, topCandi.peptideInfo.seq, 0)); //isdecoy
            double finalScore = topCandi.pepScore*(topCandi.protScore+1);
            if (!usePfmAndReduceDb) {
                finalScore = topCandi.pepScore; // dont use pfm
            }
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
        BufferedWriter writer = new BufferedWriter(new FileWriter(outputDir+"."+hostName+".PFM.csv"));
        writer.write("scanName,scanNum,qValue,TorD,finalScore,peptide,freeSeq,pepScore,proteins,protScore,ppm,theoMass,expMass,charge\n");
        for (Pair<Double, String> pair : finalExcelList) {
            writer.write(pair.getSecond());
        }
        writer.close();
//        logger.info("Number of PSMs with q < 0.01: {}", numScansQ001);
    }

    class TempRes{
        public double pepScore;
        public double protScore;
        public double finalScore;
        public double qValue;
        public boolean isDecoy;
        public String purePepSeq;
        public String ptmPepSeq;
        public String protIdSet;
        public double theoMass;
        public double expMass;
        public double massDiffPpm;
        TempRes(double pepScore, double protScore, double qValue, boolean isDecoy, String ptmPepSeq, double expMass){
            this.pepScore = pepScore;
            this.protScore = protScore;
            this.finalScore = pepScore*Math.sqrt(protScore+1);
            this.qValue = qValue;
            this.isDecoy = isDecoy;
            this.ptmPepSeq = ptmPepSeq;
            this.purePepSeq = ptmPepSeq; //todo
            this.theoMass = 0d;
            this.expMass = expMass;
            this.massDiffPpm = 0d;
        }
    }
    private void writeFinalResult(Map<String, PercolatorEntry> percolatorResultMap, String outputPath, Map<String, PeptideInfo> allPeptideInfoMap, String sqlPath, String outputDir, String hostName) throws IOException, SQLException {
        TreeMap<Double, List<String>> tempMap = new TreeMap<>();
        String resultPath = outputDir + "." + hostName +".candidates.csv";
        BufferedWriter dbgWriter = new BufferedWriter(new FileWriter(resultPath));
        dbgWriter.write("scanNum,pep1,s1,prot1,pep2,s2,prot2,pep3,s3,prot3,pep4,s4,prot4,pep5,s5,prot5,pep6,s6,prot6,pep7,s7,prot7,pep8,s8,prot8,pep9,s9,prot9,pep10,s10,prot10,pep11,s11,prot11,pep12,s12,prot12,pep13,s13,prot13,pep14,s14,prot14,pep15,s15,prot15,pep16,s16,prot16,pep17,s17,prot17,pep18,s18,prot18,pep19,s19,prot19,pep20,s20,prot20\n");

        BufferedWriter writer = new BufferedWriter(new FileWriter(outputPath));

        if (percolatorResultMap == null) {
            writer.write("scanName,scan_num,peptide,isDecoy, shouldPtm, hasPTM, ptmNum, isSettled,charge,theo_mass,exp_mass,abs_ppm,A_score,protein_ID,score,delta_C_n,deltaLCn,globalRank, " +
                    "cosScore,ionFrac,matchedHighestIntensityFrac,explainedAaFrac, other_PTM_patterns,MGF_title,labelling,isotope_correction,MS1_pearson_correlation_coefficient\n");
        } else {
            writer.write("scanName,scan_num,peptide, isDecoy, q_value, score,delta_C_n, theo_mass,exp_mass,abs_ppm,shouldPtm, hasPTM, ptmNum, isSettled,charge,A_score,protein_ID,deltaLCn,globalRank, " +
                    "cosScore,ionFrac,matchedHighestIntensityFrac,explainedAaFrac,percolator_score,posterior_error_prob,other_PTM_patterns,MGF_title,labelling,isotope_correction,MS1_pearson_correlation_coefficient\n");
        }

        Connection sqlConnection = DriverManager.getConnection(sqlPath);
        Statement sqlStatement = sqlConnection.createStatement();
        ResultSet sqlResultSet = sqlStatement.executeQuery("SELECT scanName,scanNum, shouldPtm, hasPTM, ptmNum, isSettled, precursorCharge, precursorMass, mgfTitle, isotopeCorrectionNum, ms1PearsonCorrelationCoefficient, " +
                "labelling, peptide, theoMass, isDecoy, score, otherPtmPatterns, aScore, deltaCn,deltaLCn, isDecoy, globalRank, normalizedCorrelationCoefficient, matchedPeakNum, ionFrac, matchedHighestIntensityFrac, explainedAaFrac, peptideSet FROM spectraTable");
        int numScansQ001 = 0;

        while (sqlResultSet.next()) {
            int isDecoy = sqlResultSet.getInt("isDecoy");
            if (!sqlResultSet.wasNull()) {
//                if (isDecoy == 0) {
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
                String peptideSet = sqlResultSet.getString("peptideSet");

                PeptideInfo pepInfo = allPeptideInfoMap.get(peptide.replaceAll("[^A-Z]+", ""));
                TreeSet<String> proteinIdSet = new TreeSet<>();
                for (String protein : pepInfo.protIdSet) {
                    proteinIdSet.add(protein.trim());
                }
                String[] test = new String[]{"a","s"};
                String aScore = sqlResultSet.getString("aScore");
                dbgWriter.write(scanNum  + "," + peptideSet + "\n");

                if (percolatorResultMap == null) {
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
                    if ((isDecoy == 0 && percolatorEntry.isDecoy) || (isDecoy == 1 && !percolatorEntry.isDecoy)){
//                        System.out.println("wrong isDecoy" + scanName);
                    }

                    String str = String.format(Locale.US, "%s, %d,%s,%d,%s,%f,%f,%f,%f,%f,%d,%d,%d,%d,%d,%s,%s,%f,%d,%f,%f,%f,%f,%f,%s,%s,\"%s\",%s,%d,%f\n"
                            , scanName, scanNum, peptide,percolatorEntry.isDecoy ? 1 : 0, percolatorEntry.qValue, score, sqlResultSet.getDouble("deltaCn"),theoMass, expMass, ppm
                            , sqlResultSet.getInt("shouldPtm"),sqlResultSet.getInt("hasPTM"),sqlResultSet.getInt("ptmNum"),sqlResultSet.getInt("isSettled")
                            , sqlResultSet.getInt("precursorCharge"),  aScore, String.join(";", proteinIdSet).replaceAll(",", "~"), deltaLCn,globalRank,cosScore, ionFrac, matchedHighestIntensityFrac,explainedAaFrac
                            , percolatorEntry.percolatorScore, percolatorEntry.PEP, sqlResultSet.getString("otherPtmPatterns"), sqlResultSet.getString("mgfTitle")
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
        dbgWriter.close();
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
    private static Map<String, PercolatorEntry> runPercolator(String percolatorPath, String percolatorInputFileName, String percolatorOutputFileName, String percolatorDecoyOutputFileName, String percolatorProteinOutputFileName, String tdFastaPath, String enzymeName) throws Exception {
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

            String[] commands = new String[]{percolatorPath, "--only-psms", "--verbose", "0", "--no-terminate", "--search-input", "concatenated", "--protein-decoy-pattern", "DECOY_", "--picked-protein", tdFastaPath, "--protein-enzyme", percolatorEnzyme, "--protein-report-fragments", "--protein-report-duplicates", "--results-proteins", percolatorProteinOutputFileName, "--results-psms", percolatorOutputFileName,"--decoy-results-psms", percolatorDecoyOutputFileName, percolatorInputFileName};

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
                    percolatorResultMap.put(parts[0], new PercolatorEntry(Double.valueOf(parts[1]), parts[2], parts[3], false));
                }
            }
            BufferedReader decoyReader = new BufferedReader(new FileReader(percolatorDecoyOutputFileName));
            String decoyLine;
            while ((decoyLine = decoyReader.readLine()) != null) {
                decoyLine = decoyLine.trim();
                if (!decoyLine.startsWith("PSMId")) {
                    String[] parts = decoyLine.split("\t");
                    percolatorResultMap.put(parts[0], new PercolatorEntry(Double.valueOf(parts[1]), parts[2], parts[3], true));
                }
            }
            reader.close();
            decoyReader.close();
        } else {
            logger.error("Cannot find Percolator input file (from {}) for estimating Percolator Q-Value.", percolatorInputFileName);
            return percolatorResultMap;
        }

        return percolatorResultMap;
    }

//    private void writeDebugResults(String sqlPath, String spectraPath) throws IOException, SQLException {
//        Connection sqlConnection = DriverManager.getConnection(sqlPath);
//        Statement sqlStatement = sqlConnection.createStatement();
//        ResultSet sqlResultSet = sqlStatement.executeQuery("SELECT scanNum, candidates, whereIsTopCand FROM spectraTable");
//
//        String resultPath = spectraPath + "candidates.csv";
//        try (BufferedWriter writer = new BufferedWriter(new FileWriter(resultPath))) {
//            writer.write("scanNo,whereIsTopCand,pep1,s1,b1,m1,p1,pep2,s2,b2,m2,p2,pep3,s3,b3,m3,p3,pep4,s4,b4,m4,p4,pep5,s5,b5,m5,p5,pep6,s6,b6,m6,p6,pep7,s7,b7,m7,p7,pep8,s8,b8,m8,p8,pep9,s9,b9,m9,p9,pep10,s10,b10,m10,p10,pep11,s11,b11,m11,p11,pep12,s12,b12,m12,p12,pep13,s13,b13,m13,p13,pep14,s14,b14,m14,p14,pep15,s15,b15,m15,p15,pep16,s16,b16,m16,p16,pep17,s17,b17,m17,p17,pep18,s18,b18,m18,p18,pep19,s19,b19,m19,p19,pep20,s20,b20,m20,p20\n");
//            while (sqlResultSet.next()) {
//                int scanNum = sqlResultSet.getInt("scanNum");
//                int whereIsTopCand = sqlResultSet.getInt("whereIsTopCand");
//                String candidatesList = sqlResultSet.getString("candidates");
//                if (!sqlResultSet.wasNull()) {
//                    writer.write(scanNum + "," + whereIsTopCand + "," + candidatesList + "\n");
//                }
//            }
//        } catch (IOException ex) {
//            ex.printStackTrace();
//            logger.error(ex.getMessage());
//            System.exit(1);
//        }
//        sqlResultSet.close();
//        sqlStatement.close();
//        sqlConnection.close();
//    }


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
