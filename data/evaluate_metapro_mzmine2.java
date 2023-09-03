package net.sf.mzmine.modules.peaklistmethods.alignment.ransac;

import com.google.common.collect.BoundType;
import com.google.common.collect.Range;
import com.opencsv.CSVWriter;
import net.sf.mzmine.datamodel.*;
import net.sf.mzmine.datamodel.impl.SimpleDataPoint;
import net.sf.mzmine.datamodel.impl.SimpleFeature;
import net.sf.mzmine.datamodel.impl.SimplePeakList;
import net.sf.mzmine.datamodel.impl.SimplePeakListRow;
import net.sf.mzmine.modules.peaklistmethods.alignment.join.JoinAlignerParameters;
import net.sf.mzmine.modules.peaklistmethods.alignment.join.JoinAlignerTask;
import net.sf.mzmine.parameters.ParameterSet;
import net.sf.mzmine.parameters.parametertypes.selectors.PeakListsSelectionType;
import net.sf.mzmine.parameters.parametertypes.tolerances.MZTolerance;
import net.sf.mzmine.parameters.parametertypes.tolerances.RTTolerance;
import net.sf.mzmine.project.impl.MZmineProjectImpl;
import net.sf.mzmine.project.impl.RawDataFileImpl;
import weka.core.pmml.jaxbbindings.False;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;

public class evaluate_mzmine2 {

    public static void eval(String[] csvFilePaths, String outputFilePath, double mz_tolerance, double ppm_mz_tolerance,
                            double rt_tolerance_before, double rt_tolerance_after, double margin, String method) {
        RawDataFile[] rawDataFiles = new RawDataFile[csvFilePaths.length];
        MZmineProject project = new MZmineProjectImpl();
        PeakList[] peakLists = new PeakList[csvFilePaths.length];
        try {
            for (int i = 0; i < csvFilePaths.length; i ++) {
                String[] pathSplit = csvFilePaths[i].split("\\\\");
                String fileName = pathSplit[pathSplit.length - 1].split("\\.")[0];
                RawDataFile rawDataFile = new RawDataFileImpl(fileName);
                rawDataFiles[i] = rawDataFile;
                PeakList peakList = new SimplePeakList(fileName, rawDataFile);
                BufferedReader reader = new BufferedReader(new FileReader(csvFilePaths[i]));
                int rowId = 0;
                String line = null;
                while ((line = reader.readLine()) != null) {
                    String[] lineSplit = line.strip().split(",");
                    PeakListRow row = new SimplePeakListRow(rowId);
                    Feature peak = new SimpleFeature(rawDataFile, Double.parseDouble(lineSplit[0]),
                            Double.parseDouble(lineSplit[1]), Double.parseDouble(lineSplit[2]),
                            Double.parseDouble(lineSplit[2]), null, new DataPoint[]{new SimpleDataPoint(0,0)},
                            Feature.FeatureStatus.DETECTED, 0, 0,
                            null, null, null,
                            Range.range(0.0, BoundType.CLOSED, Double.parseDouble(lineSplit[2]), BoundType.CLOSED));
                    row.addPeak(peak.getDataFile(), peak);
                    peakList.addRow(row);
                    rowId ++;
                }
                peakLists[i] = peakList;
            }
        } catch (Exception e) {
            e.printStackTrace();
        }

        MZTolerance mzTolerance  = new MZTolerance(mz_tolerance, ppm_mz_tolerance);
        RTTolerance rtToleranceBefore = new RTTolerance(true, rt_tolerance_before);
        RTTolerance rtToleranceAfter = new RTTolerance(true, rt_tolerance_after);
        switch (method) {
            case "ransac":
                ParameterSet ransacParams = new RansacAlignerParameters();
                ransacParams.getParameter(RansacAlignerParameters.MZTolerance).setValue(mzTolerance);
                ransacParams.getParameter(RansacAlignerParameters.RTToleranceBefore).setValue(rtToleranceBefore);
                ransacParams.getParameter(RansacAlignerParameters.RTToleranceAfter).setValue(rtToleranceAfter);

                ransacParams.getParameter(RansacAlignerParameters.peakListName).setValue("RANSAC");
                ransacParams.getParameter(RansacAlignerParameters.SameChargeRequired).setValue(false);

                ransacParams.getParameter(RansacAlignerParameters.Iterations).setValue(100000);
                ransacParams.getParameter(RansacAlignerParameters.NMinPoints).setValue(0.5);
                ransacParams.getParameter(RansacAlignerParameters.Margin).setValue(margin);
                ransacParams.getParameter(RansacAlignerParameters.Linear).setValue(false);

                RansacAlignerTask ransacTask = new RansacAlignerTask(project, peakLists, ransacParams);
                ransacTask.run();
                break;
            case "join":
                ParameterSet joinParams = new JoinAlignerParameters();
                joinParams.getParameter(JoinAlignerParameters.MZTolerance).setValue(mzTolerance);
                joinParams.getParameter(JoinAlignerParameters.MZWeight).setValue(1d);
                joinParams.getParameter(JoinAlignerParameters.RTTolerance).setValue(rtToleranceAfter);
                joinParams.getParameter(JoinAlignerParameters.RTWeight).setValue(1d);

                joinParams.getParameter(JoinAlignerParameters.peakListName).setValue("Join");
                joinParams.getParameter(JoinAlignerParameters.peakLists).setValue(PeakListsSelectionType.SPECIFIC_PEAKLISTS, peakLists);
                joinParams.getParameter(JoinAlignerParameters.SameChargeRequired).setValue(false);
                joinParams.getParameter(JoinAlignerParameters.SameIDRequired).setValue(false);
                JoinAlignerTask joinTask = new JoinAlignerTask(project, joinParams);
                joinTask.run();
                break;
        }



        PeakListRow[] resultList = project.getPeakLists()[0].getRows();
        try {
            CSVWriter csvWriter = new CSVWriter(new FileWriter(outputFilePath));
            String[] firstLine = new String[4 + 3 * csvFilePaths.length];
            firstLine[0] = "mz";
            firstLine[1] = "rt";
            firstLine[2] = "area";
            firstLine[3] = "#";
            for (int i = 0; i < rawDataFiles.length; i ++) {
                firstLine[4 + 3 * i] = rawDataFiles[i].getName() + "_mz";
                firstLine[5 + 3 * i] = rawDataFiles[i].getName() + "_rt";
                firstLine[6 + 3 * i] = rawDataFiles[i].getName() + "_area";
            }
            csvWriter.writeNext(firstLine, false);
            for (PeakListRow peakListRow: resultList) {
                if (peakListRow.getPeaks() == null || peakListRow.getPeaks().length == 0) {
                    continue;
                }
                String[] line = new String[4 + 3 * csvFilePaths.length];
                double mzSum = 0d, rtSum = 0d, areaSum = 0d;
                for (int i = 0; i < rawDataFiles.length; i ++) {
                    Feature peak = peakListRow.getPeak(rawDataFiles[i]);
                    if (peak != null) {
                        mzSum += peak.getMZ();
                        rtSum += peak.getRT();
                        areaSum += peak.getArea();
                        line[4 + 3 * i] = Double.toString(peak.getMZ());
                        line[5 + 3 * i] = Double.toString(peak.getRT());
                        line[6 + 3 * i] = Double.toString(peak.getArea());
                    } else {
                        line[4 + 3 * i] = "0.0";
                        line[5 + 3 * i] = "0.0";
                        line[6 + 3 * i] = "0.0";
                    }
                }
                line[0] = Double.toString(mzSum / peakListRow.getPeaks().length);
                line[1] = Double.toString(rtSum / peakListRow.getPeaks().length);
                line[2] = Double.toString(areaSum / peakListRow.getPeaks().length);
                line[3] = "0";
                csvWriter.writeNext(line, false);
            }
            csvWriter.flush();
        } catch (Exception e) {
            e.printStackTrace();
        }
        System.out.println("CSV saved");
    }

    public static void main(String[] args) {
//        String[] wiffFilePaths = new String[]{
//                "D:\\workspace\\GAligner\\metapro\\TripleTOF_6600\\SampleA_1.csv",
//                "D:\\workspace\\GAligner\\metapro\\TripleTOF_6600\\SampleA_2.csv",
//                "D:\\workspace\\GAligner\\metapro\\TripleTOF_6600\\SampleA_3.csv",
//                "D:\\workspace\\GAligner\\metapro\\TripleTOF_6600\\SampleA_4.csv",
//                "D:\\workspace\\GAligner\\metapro\\TripleTOF_6600\\SampleB_1.csv",
//                "D:\\workspace\\GAligner\\metapro\\TripleTOF_6600\\SampleB_2.csv",
//                "D:\\workspace\\GAligner\\metapro\\TripleTOF_6600\\SampleB_3.csv",
//                "D:\\workspace\\GAligner\\metapro\\TripleTOF_6600\\SampleB_4.csv"
//        };
//        String wiffOutputPath = "D:\\workspace\\GAligner\\metapro\\TripleTOF_6600_results\\TripleTOF_6600_aligned_mzmine2_ransac1.csv";
//        long startTime = System.currentTimeMillis();
//        eval(wiffFilePaths, wiffOutputPath, 0.01, 0.0, 0.3, 0.5, 0.05, "ransac");
//        System.out.println(System.currentTimeMillis() - startTime);
//        startTime = System.currentTimeMillis();
//        wiffOutputPath = "D:\\workspace\\GAligner\\metapro\\TripleTOF_6600_results\\TripleTOF_6600_aligned_mzmine2_join1.csv";
//        eval(wiffFilePaths, wiffOutputPath, 0.01, 0.0, 0.3, 0.5, 0.05, "join");
//        System.out.println(System.currentTimeMillis() - startTime);
//        startTime = System.currentTimeMillis();
//
//        String[] rawFilePaths = new String[]{
//                "D:\\workspace\\GAligner\\metapro\\QE_HF\\SA1.csv",
//                "D:\\workspace\\GAligner\\metapro\\QE_HF\\SA2.csv",
//                "D:\\workspace\\GAligner\\metapro\\QE_HF\\SA3.csv",
//                "D:\\workspace\\GAligner\\metapro\\QE_HF\\SA4.csv",
//                "D:\\workspace\\GAligner\\metapro\\QE_HF\\SA5.csv",
//                "D:\\workspace\\GAligner\\metapro\\QE_HF\\SB1.csv",
//                "D:\\workspace\\GAligner\\metapro\\QE_HF\\SB2.csv",
//                "D:\\workspace\\GAligner\\metapro\\QE_HF\\SB3.csv",
//                "D:\\workspace\\GAligner\\metapro\\QE_HF\\SB4.csv",
//                "D:\\workspace\\GAligner\\metapro\\QE_HF\\SB5.csv"
//        };
//        String rawOutputPath = "D:\\workspace\\GAligner\\metapro\\QE_HF_results\\QE_HF_aligned_mzmine2_ransac1.csv";
//        eval(rawFilePaths, rawOutputPath, 0.005, 0.0, 0.15, 0.3, 0.025, "ransac");
//        System.out.println(System.currentTimeMillis() - startTime);
//        startTime = System.currentTimeMillis();
//        rawOutputPath = "D:\\workspace\\GAligner\\metapro\\QE_HF_results\\QE_HF_aligned_mzmine2_join1.csv";
//        eval(rawFilePaths, rawOutputPath, 0.005, 0.0, 0.15, 0.3, 0.025, "join");
//        System.out.println(System.currentTimeMillis() - startTime);

        String[] mtblsFilePaths = new String[]{
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\12W-1.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\12W-2.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\12W-3.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\12W-4.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\12W-5.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\12W-6.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\12W-7.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\12W-8.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\24W-1.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\24W-2.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\24W-3.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\24W-4.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\24W-5.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\24W-6.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\24W-7.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\24W-8.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\32W-1.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\32W-2.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\32W-3.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\32W-4.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\32W-5.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\32W-6.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\32W-7.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\32W-8.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\4W-1.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\4W-2.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\4W-3.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\4W-4.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\4W-5.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\4W-6.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\4W-7.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\4W-8.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\52W-1.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\52W-2.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\52W-3.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\52W-4.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\52W-5.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\52W-6.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\52W-7.csv",
                "D:\\workspace\\GAligner\\data\\MTBLS562\\metapro\\52W-8.csv",
        };
        String mtblsOutputPath = "D:\\workspace\\GAligner\\data\\MTBLS562_results_metapro\\MTBLS562_aligned_mzmine2_ransac.csv";
        Long startTime = System.currentTimeMillis();
        eval(mtblsFilePaths, mtblsOutputPath, 0.015, 0.0, 0.1, 0.3, 0.05, "ransac");
        System.out.println(System.currentTimeMillis() - startTime);
        startTime = System.currentTimeMillis();
        mtblsOutputPath = "D:\\workspace\\GAligner\\data\\MTBLS562_results_metapro\\MTBLS562_aligned_mzmine2_join.csv";
        eval(mtblsFilePaths, mtblsOutputPath, 0.015, 0.0, 0.1, 0.3, 0.05, "join");
        System.out.println(System.currentTimeMillis() - startTime);
    }
}
