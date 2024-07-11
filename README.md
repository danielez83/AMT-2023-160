# AMT-2023-160
Support code for reproducing the figures of the article:

*A versatile water vapor generation module for vapor isotope calibration and liquid isotope measurements*

by Hans Christian Steen-Larsen and Daniele Zannoni

Accepted on Atmospheric Measurement Techniques on 25/05/2024 (Preprint available at https://doi.org/10.5194/amt-2023-160)

Note: To run the script the raw Picarro data and some extra support files are required. The extra support files are available in this reporisotory. The raw Picarro data can be retrieved here at http:....

# Figure 3
Flash evaporator 92 h stability test of the BER standard.

**Script to run:** plot_ADEV_w_wo_memory.py

**Support files required:**

- ADEV_BER17k_withmemory_CALIBRATED_R1.csv   --> Allan Deviation data without removing memory effect
- Timings_Picarro.xlsx                       --> Excel spreadsheet with time and dates of experiments, used as a lookup table to retrieve the raw data correctly

# Figure 4
Normalized step change between standards SP and BER

**Script to run:** load_and_plot_long_step_WITH_INJ.py

**Support files required:**

- SP_BER_inj_time.csv                         --> injection timings
- HKDS2156_IsoWater_20221116_165037.csv       --> Results of liquid injections with Picarro vaporizer

# Figure 5
Humidity levels selected to test the short-term performances of the vapor generation module.

**Script to run:** run_allandev_diff_hum_levels.py

**Support files required:**

- Timings_Picarro.xlsx                       --> Excel spreadsheet with time and dates of experiments, used as a lookup table to retrieve the raw data correctly

# Figure 6
Overlapping Allan deviations of $\Delta^{17}$O measured at different humidity levels and for different integration times.

**Script to run:** run_allandev_diff_hum_levels.py

**Support files required:**

- Timings_Picarro.xlsx                       --> Excel spreadsheet with time and dates of experiments, used as a lookup table to retrieve the raw data correctly

# Figure 7
Humidity–isotope characterization of the HKDS2092 analyzer between 500 and 3500 ppmv using three different isotopic standards.

**Script to run:** run_WIFVOS_test.py

**Support files required:**

- Timings_Picarro.xlsx                       --> Excel spreadsheet with time and dates of experiments, used as a lookup table to retrieve the raw data correctly

# Figure 8
The 10 h analysis of the same water vapor source using two water isotope analyzers (HKDS2092, HKDS2156).

**Script to run:** load_and_plot_PICARROs_comparison_data.py

**Support files required:**

- Timings_Picarro.xlsx                       --> Excel spreadsheet with time and dates of experiments, used as a lookup table to retrieve the raw data correctly

# Figure 9
Variability of repeated injections of BER and SP standards using the multiple ovens configuration (raw data). 

**Script to run:** plot_cal_pulses.py

**Support files required:**
- Cal_Pulses_Selector.csv                    --> Data obtained with the VICI selector configuration (only one oven working)
- Cal_Pulses_MultiOven_new20230609.csv       --> Data obtained with the multioven configuration

# Figure 10
Calibrated and normalized step change for $\Delta^{17}$O. 

**Script to run:** load_and_plot_long_step_WITH_INJ.py

**Support files required:**

- SP_BER_inj_time.csv                         --> injection timings
- HKDS2156_IsoWater_20221116_165037.csv       --> Results of liquid injections with Picarro vaporizer
