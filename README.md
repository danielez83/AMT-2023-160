# AMT-2023-160

![image](https://github.com/danielez83/AMT-2023-160/assets/7271070/5e46e5fa-aaa9-4121-a7ac-66c601523bda)
---


Support code to reproduce Figure 3 - Figure 10 from the article:

*A versatile water vapor generation module for vapor isotope calibration and liquid isotope measurements*

by Hans Christian Steen-Larsen and Daniele Zannoni

Accepted in Atmospheric Measurement Techniques on 25/05/2024 (Preprint available at https://doi.org/10.5194/amt-2023-160)

Figure numbers refer to the final (peer-reviewed and accepted) version of the manuscript.

Note: Support data is required to run the scripts. The latter can be retrieved on Zenodo at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12730778.svg)](https://doi.org/10.5281/zenodo.12730778)

---
# Suggested directory structure
---
Python scripts and the data should be organized in the following way to run the code seamlessly:

```
.
└── Main/
    ├── AMT-2023-160/
    │   ├── article_2.mplstyle
    │   ├── article.mplstyle
    │   ├── defines.py
    │   ├── import_picarro_raw_data.py
    │   ├── load_and_plot_long_step_WITH_INJ.py
    │   ├── load_and_plot_PICARROs_comparison_data.py
    │   ├── plot_ADEV_w_wo_memory.py
    │   ├── plot_cal_pulses.py
    │   ├── run_allandev_diff_hum_levels.py
    │   ├── run_WIFVOS_test.py
    │   ├── ADEV_BER17k_withmemory_CALIBRATED_R1.csv *
    │   ├── Cal_Pulses_MultiOven_new20230609.csv *
    │   ├── Cal_Pulses_Selector.csv *
    │   ├── HKDS2156_IsoWater_20221116_165037.csv *
    │   ├── SP_BER_inj_time.csv *
    │   └── Timings_Picarro.xlsx *
    └── DATA/
        ├── HKDS2092/ *
        │   └── 2022/
        │       └── .../
        └── HKDS2156/ *
            └── 2022/
                └── .../
```

Starred items indicate files and directories that can be downloaded at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.12730778.svg)](https://doi.org/10.5281/zenodo.12730778).

---
# Python packages reguired
---
The following packages are required:
- numpy
- pandas
- matplotlib
- scipy
- allantools

---
# Description of Figures
---

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
