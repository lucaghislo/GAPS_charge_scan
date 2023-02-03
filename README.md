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

Once launched, the script requires to input the **complete filepath** of the necessary files required to compute the charge scan. The filepath can be supplied either sorrounded by single ('') or double ("") quotation marks, but can also be written without.

#### Charge scan (or threshold scan) filepath

The first file required is the one produced by either the manual test tab of the GAPS_ModuleTester.py Python script or the main GAPS beckend DAQ.

```
*** GAPS CHARGE SCAN TOOL v1.0 ***

    Charge or threshold scan filepath:
```
