# TimingCalibration
Timing calibration repository for AARDVARCv3 at Syracuse

The code "TimingCalibration.py" generates two time calibration files (even and odd channels) for a single channnel of the AARDVARCv3
The data used to generate time calibration is created by sending a fast rise Gaussian pulse (1ns rise time currently) to the AARDVARC with a frequency of ~10kHz. Other frequencies can be used, however, it is important not to match the frequency of this pulse to the internal clock of the AARDVARC. The AARDVARC is set to internally trigger on this pulse such that the the pulse crosses a threshold at a random location within a 128 sample window.
Once generated, the time calibration files can be used in "TOA_withCalibration.py" to get a calibrated TOA using either a pulse to pulse measurement or a single pulse measurement.
