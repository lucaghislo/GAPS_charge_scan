# GAPS charge scan analysis tool

## Python CLI (Command Line Interface)

This implementation allows to compute **charge scan** or **threshold scan** data acquired via the manual test tab of the GAPS_ModuleTester.py script or via the beckend DAQ.

Charge scan data file **must** be compliant to the following configuration:

```
# Event data block #1
# Threshold DAC Events  Triggered   CH#
205 0   1000    0   0
[...]
205 0   1000    0   31
```

### Download

Download the [latest release](https://github.com/lucaghislo/GAPS_charge_scan/releases/) suitable for you operating system (Windows or Linux).

### Usage

Once launched, the script requires to input the **complete filepath** of the necessary files required to compute the charge scan (or threshold scan). The filepath can be supplied either sorrounded by single ('') or double ("") quotation marks, but can also be written without. The program check for user input correctness and, in case it is not, the program asks for user input again.

#### Charge scan (or threshold scan) filepath

The first file required is the one produced by either the manual test tab of the GAPS_ModuleTester.py Python script or the main GAPS beckend DAQ. In case of threshold scan, the currently only supported file configuration is the one produced by the manual test of the  GAPS_ModuleTester.py script.

```
*** GAPS CHARGE SCAN TOOL v1.0 ***

    Charge or threshold scan filepath: C:\path\charge_threshold_scan_dummy.txt
```

#### Channel range

Next, the user is required to input the first and last channel for which the script should compute the charge scan. The minimum accepted value is **0** and the maximum is **31**.

```
                        First channel: 0
                         Last channel: 31                  
```

#### Excluded channels

Next, a list of excluded channels can be supplied. This can be done by listing the channels to be excluded separated by commas (with or without spacing). In case no channel needs to be removed, leave it blank.

```
 Excluded channels (comma separated): 7, 8, 9
```

In case all specified channels need to be analysed, leave the space blank and press enter.

```
 Excluded channels (comma separated):
```

#### Output folder

The user is required to input the path of the folder where the results will be stored.

```
               Output folder filepath: C:\path\output_folder  
```

#### Channel deactivation based on minimum threshold

The user can input the threshold, in keV, below which channels are deactivated.

```
      Deactivate channels below [keV]: 30 
```

#### Channel deactivation based on maximum ENC

The user can input the ENC, in keV, above which channels are deactivated.

```
 Deactivate channels with ENC > [keV]: 5
```

#### Parasitic injection compensation

Next, the user can decide to compensate the parasitic injection (**y**) or not (**n**). In case the user decides not to not compensate it, computation starts and the results are stored in the specified output path.

```
Compensate parasitic injection? (y/n): n

CHARGE SCAN

Working on it, be patient...
```

In case parasitic injection has to be compensated, the user is required to specify the **pedestal** and **transfer function** filepaths for the module under test. These two files can be found in the data/ folder of the automated test performed on the module.

```
Compensate parasitic injection? (y/n): y
```

The two required files are **Pedestals.dat** and **TransferFunction.dat** obtained after computing the automated test with the specific MATLAB or Python script.

```
         Pedestal from automated test: C:\path\Pedestals.dat
Transfer function from automated test: C:\path\TransferFunction.dat
```

Next, the peaking time at which the charge scan has been performed has to be specified. The only acceptable values range from **0** to **7**. After that, charge scan is computed and the results are stored in the specified output path.

```
                Peaking time (0 to 7): 5

CHARGE SCAN

Working on it, be patient...
```

### Results

Results are stored in the specified output path. In case the user decided not to compensate the parasitic injection, the output folder is organised as follows:

```
output_folder/
├── ENC_THR/
│   ├── ch#-#_ENC.pdf
│   ├── ch#-#_ENC.png
│   ├── ch#-#_THR_ENC.dat
│   ├── ch#-#_THR_hist.pdf
│   ├── ch#-#_THR_hist.png
│   ├── ch#-#_THR_plot.pdf
│   └── ch#-#_THR_plot.png
├── single_channels/
│   ├── data/
│   │   ├── ch_#_THR_###.dat
│   │   ├── [...]
│   │   └── ch_#_THR_###.dat
│   └── plots/
│       ├── charge_scan_ch#_THR_###.pdf
│       ├── charge_scan_ch#_THR_###.png
│       ├── [...]
│       ├── charge_scan_ch#_THR_###.pdf
│       └── charge_scan_ch#_THR_###.png
├── ch#-#_activation_mask.txt
├── charge_scan_ch#-#.pdf
└── charge_scan_ch#-#.png
```

In case parasitic injection has been compensated, the script adds the compensated version of every file listed above, identifiable by the "_inj" notation added to every filename.

```
output_folder/
├── ENC_THR/
│   ├── ch#-#_ENC.pdf
│   ├── ch#-#_ENC.png
│   ├── ch#-#_ENC_inj.pdf
│   ├── ch#-#_ENC_inj.png
│   ├── ch#-#_THR_ENC.dat
│   ├── ch#-#_THR_ENC_inj.dat
│   ├── ch#-#_THR_hist.pdf
│   ├── ch#-#_THR_hist.png
│   ├── ch#-#_THR_hist_inj.pdf
│   ├── ch#-#_THR_hist_inj.png
│   ├── ch#-#_THR_plot.pdf
│   ├── ch#-#_THR_plot.png
│   ├── ch#-#_THR_plot_inj.pdf
│   └── ch#-#_THR_plot_inj.png
├── single_channels/
│   ├── data/
│   │   ├── ch_#_THR_###.dat
│   │   ├── ch_#_THR_###_inj.dat
│   │   ├── [...]
│   │   ├── ch_#_THR_###.dat
│   │   └── ch_#_THR_###_inj.dat
│   └── plots/
│       ├── charge_scan_ch#_THR_###.pdf
│       ├── charge_scan_ch#_THR_###_inj.pdf
│       ├── charge_scan_ch#_THR_###.png
│       ├── charge_scan_ch#_THR_###_inj.png
│       ├── [...]
│       ├── charge_scan_ch#_THR_###.pdf
│       ├─ charge_scan_ch#_THR_###_inj.pdf
│       ├── charge_scan_ch#_THR_###.png
│       └── charge_scan_ch#_THR_###_inj.png
├── ch#-#_activation_mask.txt
├── charge_scan_ch#-#.pdf
├── charge_scan_ch#-#_inj.pdf
├── charge_scan_ch#-#.png
├── charge_scan_ch#-#_inj.png
└── summary_inj_ch#-#.dat
```

#### Output folder content

Results are organised as follows:

- `charge_scan_ch#-#.pdf`: charge scan overview for all channels in the specified range.

![FTHR_THR_205_pt5_ch_0-31_1](https://user-images.githubusercontent.com/36998696/216663752-847b4119-5e77-4d80-94c3-25ea902bf05c.png)

- `charge_scan_ch#-#_inj.pdf`: charge scan overview with parasitic injection removed for all channels in the specified range.

![FTHR_THR_205_pt5_ch_0-31_inj_1](https://user-images.githubusercontent.com/36998696/216663898-0e7cc452-ee2c-46e3-a66c-d68308d760f6.png)

- `single_channels/data/ch_#_THR_###.dat` and `single_channels/data/ch_#_THR_###_inj.dat`: raw charge scan data with and without parasitic injection compensation for channel # at DAC_thr_code ###. This is done for all channels in the specified range. First column is the energy [keV] and the second column is the trigger probability [%].

```
   0.000000     0.000000
   0.841000     0.000000
   [...]
   251.459000   100.000000
   252.300000   100.000000

```

- `single_channels/plots/charge_scan_ch#_THR_###.pdf` and `single_channels/plots/charge_scan_ch#_THR_###_inj.pdf`: raw charge scan data with and without parasitic injection compensation for channel # at DAC_thr_code ###. This is done for all channels in the specified range.

![charge_scan_ch0_THR_205_1](https://user-images.githubusercontent.com/36998696/216658156-3af2fe82-0a2f-48aa-a78c-fee9abf448a8.png)

- `ENC_THR/ch#-#_ENC.pdf` and `ENC_THR/ch#-#_ENC_inj.pdf`: estimated ENC plot for all channels in the specified range, with and without parasitic injection compensation.

![ch0-31_ENC_1](https://user-images.githubusercontent.com/36998696/216660475-0f0f9cf3-0a71-4284-b845-717a58717a7f.png)

- `ENC_THR/ch#-#_ENC.dat` and `ENC_THR/ch#-#_ENC_inj.dat`: estimated threshold and ENC values for all channels in specified range, with and without parasitic injection compensation. The first column is the channel, second column is the threshold [keV] and the third column is the ENC [keV].

```
ch  thr           enc
0   69.736702     5.965096
[...]
31  112.185906    4.861941

```

- `ENC_THR/ch#-#_THR_hist.pdf` and `ENC_THR/ch#-#_THR_hist_inj.pdf`: estimated threshold histogram for all channels in specified range, with and without parasitic injection compensation.

![ch0-31_THR_hist_1](https://user-images.githubusercontent.com/36998696/216660727-ce348ced-3ab5-42e1-820b-1757da88f878.png)

- `ENC_THR/ch#-#_THR_plot.pdf` and `ENC_THR/ch#-#_THR_plot_inj.pdf`: estimated threshold plot for all channels in specified range, with and without parasitic injection compensation.

![ch0-31_THR_plot_1](https://user-images.githubusercontent.com/36998696/216660939-9184ca9f-463d-4aff-8fd8-d3c7a6e68dab.png)

- `summary_inj_ch#-#.dat` and `summary_inj_ch#-#_inj.dat`: estimated transfer function gain *tf_gain* [ADU/keV] and pedestal *tf_pedestal* [ADU] for all channels in specified range, with and without parasitic injection compensation. For every channel it's also provided the pedestal from automated test *auto_pedestal* [ADU]. The last column, *par_inj*, refers to the estimated parasitic injection [keV].

```
ch  tf_gain  tf_pedestal auto_pedestal   par_inj
0   1.01966  150.30685   139.243         11.06389
[...]
31  1.08302  156.381     141.505         14.876
```

- `ch#-#_activation_mask.txt` and `ch#-#_activation_mask_inj.txt`: channel activation mask. Channels below the threshold value provided by the user are automatically disabled. Unresponsive channels are also set to 0. Mask is expressed in hexadecimal format with **channels listed from 31 to 0**.

The example below referers to binary mask `10011111111111111111111111111111` where channels 29 and 30 have been deactivated based on the above mentioned criteria.

```
0x9fffffff
```
