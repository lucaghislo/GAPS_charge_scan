# GAPS charge scan analysis tool

## Python CLI

This implementation allows to compute **charge scan** or **threshold scan** data acquired via the manual test tab of the GAPS_ModuleTester.py script or via the beckend DAQ.

Charge scan data file **must** be compliant to the following configuration:

```
# Event data block #1
# Threshold DAC Events  Triggered   CH#
205 0   1000    0   0
[...]
205 0   1000    0   31
```

### Instructions

Run [main.py](python_script_CLI/main.py) located in [python_script_CLI](python_script_CLI) folder.
