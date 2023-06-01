import os
import glob
import time

from pyteomics import mzml


class Spectrum:
    def __init__(self, rt, mzs, intensities, base_peak_intensity, base_peak_mz, total_ion_current):
        self.rt = rt
        self.mzs = mzs
        self.intensities = intensities
        self.base_peak_intensity = base_peak_intensity
        self.base_peak_mz = base_peak_mz
        self.total_ion_current = total_ion_current


class RawFileReader:
    def __init__(self, raw_params):
        self.min_intensity = max(0, raw_params.min_intensity)

    def load_ms1_spectra(self, result_path):
        glob_paths = glob.glob(result_path.split('.')[0] + '*.[mM][zZ][mM][lL]')
        assert len(glob_paths) == 1, 'File Error! Result and RAW file name and number should be the same. '
        raw_path = glob_paths[0]

        ms1_spectra = []
        ms_file = mzml.read(raw_path)
        for spectrum in ms_file:
            if spectrum.get('ms level') != 1:
                continue
            mzs = spectrum.get('m/z array')
            intensities = spectrum.get('intensity array')
            rt = spectrum.get('scanList').get('scan')[0].get('scan start time')
            idx = intensities > self.min_intensity
            mzs = mzs[idx]
            intensities = intensities[idx]
            base_peak_intensity = spectrum.get('base peak intensity')
            base_peak_mz = spectrum.get('base peak m/z')
            total_ion_current = spectrum.get('total ion current')
            ms1_spectra += [Spectrum(rt, mzs, intensities, base_peak_intensity, base_peak_mz, total_ion_current)]
        return ms1_spectra
