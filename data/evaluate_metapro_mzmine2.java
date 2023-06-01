package net.sf.mzmine.modules.peaklistmethods.alignment.ransac;

import com.google.common.collect.BoundType;
import com.google.common.collect.Range;
import com.opencsv.CSVWriter;
import net.sf.mzmine.datamodel.*;
import net.sf.mzmine.datamodel.impl.SimpleDataPoint;
import net.sf.mzmine.datamodel.impl.SimpleFeature;
import net.sf.mzmine.datamodel.impl.SimplePeakList;
import net.sf.mzmine.datamodel.impl.SimplePeakListRow;
import net.sf.mzmine.parameters.ParameterSet;
import net.sf.mzmine.parameters.parametertypes.tolerances.MZTolerance;
import net.sf.mzmine.parameters.parametertypes.tolerances.RTTolerance;
import net.sf.mzmine.project.impl.MZmineProjectImpl;
import net.sf.mzmine.project.impl.RawDataFileImpl;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;

public class evaluate_mzmine2 {

    public static void eval(String[] csvFilePaths, String outputFilePath, double mz_tolerance, double ppm_mz_tolerance,
                            double rt_tolerance, double margin) {
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

        ParameterSet parameters = new RansacAlignerParameters();
        MZTolerance mzTolerance  = new MZTolerance(mz_tolerance, ppm_mz_tolerance);
        RTTolerance rtToleranceBefore = new RTTolerance(true, rt_tolerance);
        RTTolerance rtToleranceAfter = new RTTolerance(true, rt_tolerance);

        parameters.getParameter(RansacAlignerParameters.MZTolerance).setValue(mzTolerance);
        parameters.getParameter(RansacAlignerParameters.RTToleranceBefore).setValue(rtToleranceBefore);
        parameters.getParameter(RansacAlignerParameters.RTToleranceAfter).setValue(rtToleranceAfter);

        parameters.getParameter(RansacAlignerParameters.peakListName).setValue("RANSAC");
        parameters.getParameter(RansacAlignerParameters.SameChargeRequired).setValue(false);

        parameters.getParameter(RansacAlignerParameters.Iterations).setValue(100000);
        parameters.getParameter(RansacAlignerParameters.NMinPoints).setValue(0.5);
        parameters.getParameter(RansacAlignerParameters.Margin).setValue(margin);
        parameters.getParameter(RansacAlignerParameters.Linear).setValue(false);

        RansacAlignerTask task = new RansacAlignerTask(project, peakLists, parameters);
        task.run();
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
        String[] wiffFilePaths = new String[]{
                "D:\\workspace\\GAligner\\data\\TripleTOF_6600\\metapro\\SampleA_1.csv",
                "D:\\workspace\\GAligner\\data\\TripleTOF_6600\\metapro\\SampleA_2.csv",
                "D:\\workspace\\GAligner\\data\\TripleTOF_6600\\metapro\\SampleA_3.csv",
                "D:\\workspace\\GAligner\\data\\TripleTOF_6600\\metapro\\SampleA_4.csv",
                "D:\\workspace\\GAligner\\data\\TripleTOF_6600\\metapro\\SampleB_1.csv",
                "D:\\workspace\\GAligner\\data\\TripleTOF_6600\\metapro\\SampleB_2.csv",
                "D:\\workspace\\GAligner\\data\\TripleTOF_6600\\metapro\\SampleB_3.csv",
                "D:\\workspace\\GAligner\\data\\TripleTOF_6600\\metapro\\SampleB_4.csv"
        };
        String wiffOutputPath = "D:\\workspace\\GAligner\\data\\TripleTOF_6600_results_metapro\\TripleTOF_6600_aligned_mzmine2.csv";
        eval(wiffFilePaths, wiffOutputPath, 0.01, 0.0, 0.5, 0.1);

        String[] rawFilePaths = new String[]{
                "D:\\workspace\\GAligner\\data\\QE_HF\\metapro\\SA1.csv",
                "D:\\workspace\\GAligner\\data\\QE_HF\\metapro\\SA2.csv",
                "D:\\workspace\\GAligner\\data\\QE_HF\\metapro\\SA3.csv",
                "D:\\workspace\\GAligner\\data\\QE_HF\\metapro\\SA4.csv",
                "D:\\workspace\\GAligner\\data\\QE_HF\\metapro\\SA5.csv",
                "D:\\workspace\\GAligner\\data\\QE_HF\\metapro\\SB1.csv",
                "D:\\workspace\\GAligner\\data\\QE_HF\\metapro\\SB2.csv",
                "D:\\workspace\\GAligner\\data\\QE_HF\\metapro\\SB3.csv",
                "D:\\workspace\\GAligner\\data\\QE_HF\\metapro\\SB4.csv",
                "D:\\workspace\\GAligner\\data\\QE_HF\\metapro\\SB5.csv"
        };
        String rawOutputPath = "D:\\workspace\\GAligner\\data\\QE_HF_results_metapro\\QE_HF_aligned_mzmine2.csv";
        eval(rawFilePaths, rawOutputPath, 0.005, 0.0, 0.3, 0.05);

    }
}
