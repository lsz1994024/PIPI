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
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.concurrent.*;
import java.util.concurrent.locks.ReentrantLock;

public class PIPI {
    private static final Logger logger = LoggerFactory.getLogger(PIPI.class);
    public static final String versionStr = "1.4.7";
    static final boolean useXcorr = false;
    static final int minTagLenToReduceProtDb = 5;
    public static final boolean isPtmSimuTest = false; //normal //todo
    static final boolean usePfmAndReduceDb = true;  //normal //todo
    public    static final double MIN_PEAK_SUM_INFER_AA = 0.5;
    static final int  maxNumVarPtmConsidered = 5;
    public static final int[] debugScanNumArray = new int[]{};
    public static ArrayList<Integer> lszDebugScanNum = new ArrayList<>(Arrays.asList(14075));//35581 16918, 16847,16457,16483,
    public static int neighborNum = 100;
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

        SpecProcessor specProcessor = new SpecProcessor(massTool);

        logger.info("Get long tags to reduce proteins...");
        int threadNum_0 = Integer.valueOf(parameterMap.get("thread_num"));
        if (threadNum_0 == 0) {
            threadNum_0 = 3 + Runtime.getRuntime().availableProcessors();
        }
        if (java.lang.management.ManagementFactory.getRuntimeMXBean().getInputArguments().toString().indexOf("jdwp") >= 0){
            //change thread 1
            threadNum_0 = 1;
        }
        System.out.println("thread NUM "+ threadNum_0);

        Map<String, String> mgfTitleMap = new HashMap<>();
        Map<String, Integer> isotopeCorrectionNumMap = new HashMap<>();
        Map<String, Double> ms1PearsonCorrelationCoefficientMap = new HashMap<>();

        BufferedReader bReader = new BufferedReader(new FileReader(outputDir+"truth.csv"));
        Map<Integer, String> backboneTruth = new HashMap<>();
        String line;
        while ((line = bReader.readLine()) != null) {
            line = line.trim();
            String[] splitRes = line.split(",");
            if (splitRes[0].contentEquals("scanNum") || splitRes[0].contentEquals("scanNo")) {
                continue;
            }
            String backbone = splitRes[1].replace('I','L');
            backboneTruth.put(Integer.valueOf(splitRes[0]), backbone);
        }

        //////////==================================================
        //   Get Long Tags
        ExecutorService threadPoolGetLongTag = Executors.newFixedThreadPool(threadNum_0);
        ArrayList<Future<PreSearch.Entry>> taskListGetLongTag = new ArrayList<>(datasetReader.getUsefulSpectraNum() + 10);
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

            String[] scanNameStr = scanName.split("\\.");
            boolean shouldRun = false;
            for (int debugScanNum : lszDebugScanNum) {
                if (Math.abs(scanNum-debugScanNum) < neighborNum) {
                    shouldRun = true;
                }
            }
            if (java.lang.management.ManagementFactory.getRuntimeMXBean().getInputArguments().toString().indexOf("jdwp") >= 0){
                if (!shouldRun) {
                    continue;
                }
            }
            int fileId = fileNameIdMap.get( scanNameStr[0] );
            taskListGetLongTag.add(threadPoolGetLongTag.submit(new PreSearch(scanNum, buildIndex, massTool, ms2Tolerance, ms1Tolerance, leftInverseMs1Tolerance, rightInverseMs1Tolerance
                    , ms1ToleranceUnit, inferPTM.getMinPtmMass(), inferPTM.getMaxPtmMass(), Math.min(precursorCharge > 1 ? precursorCharge - 1 : 1, 3)
                    , spectraParserArray[fileId], minClear, maxClear, lockGetLongTag, scanName, precursorCharge, precursorMass, specProcessor , backboneTruth.get(scanNum), Integer.valueOf(parameterMap.get("min_peptide_length")), Integer.valueOf(parameterMap.get("max_peptide_length")))));
        }
        System.out.println("totalSubmit in SpecCoder, "+ submitNumSpecCoderX);
        sqlResSetGetLongTag.close();
        sqlStateGetLongTag.close();

        int lastProgressGetLongTag = 0;
        int totalCountGetLongTag = taskListGetLongTag.size();
        int countGetLongTag = 0;
        Map<Integer, Map<String, Double>> scanVarLenTruthGoodnessMap = new HashMap<>();
        while (countGetLongTag < totalCountGetLongTag) {
            // record search results and delete finished ones.
            List<Future<PreSearch.Entry>> toBeDeleteTaskList = new ArrayList<>(totalCountGetLongTag - countGetLongTag);
            for (Future<PreSearch.Entry> task : taskListGetLongTag) {
                if (task.isDone()) {
                    if (task.get() != null ) {
                        PreSearch.Entry entry = task.get();
                        scanVarLenTruthGoodnessMap.put(entry.scanNum, entry.varLenTruthGoodnessMap);
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
        //write res
        List<Pair<Double, String>> finalExcelList = new ArrayList<>(scanVarLenTruthGoodnessMap.size());
        for (int scanNum : scanVarLenTruthGoodnessMap.keySet()) {
            Map<String, Double> varLenTruthGoodnessMap = scanVarLenTruthGoodnessMap.get(scanNum);
            String finalStr = String.format(Locale.US, "%d,%s,%s,%s,%s,%s,%s,%s,%s,%s\n", scanNum, backboneTruth.get(scanNum),
                    varLenTruthGoodnessMap.get("3"),varLenTruthGoodnessMap.get("4"),varLenTruthGoodnessMap.get("5"),varLenTruthGoodnessMap.get("6"),
                    varLenTruthGoodnessMap.get("7"),varLenTruthGoodnessMap.get("8"),varLenTruthGoodnessMap.get("9"),varLenTruthGoodnessMap.get("V"));
            finalExcelList.add(new Pair(varLenTruthGoodnessMap.get("V"), finalStr));
        }
        // official output with pfm
        Collections.sort(finalExcelList, Comparator.comparing(o -> o.getFirst(), Comparator.reverseOrder()));
        BufferedWriter writer = new BufferedWriter(new FileWriter(outputDir+"res.csv"));

        writer.write("scanNum,truthPep,res3,res4,res5,res6,res7,res8,res9,resV,\n");
        for (Pair<Double, String> pair : finalExcelList) {
            writer.write(pair.getSecond());
        }
        writer.close();
        logger.info("Saving results...");
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

}
