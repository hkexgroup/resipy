Example: simple flat 2D ERI survey

1. Open th GUI (pyR2):
	- Option 1: use standalone version
	- Option 2: use command line: 	cd /pyr2/src
									python ui.py

2. Select "Inverse" (checked by default)
3. Select working directory: /examples/2DCrossBorehole
4. Check "Borehole Survey" checkbox in the "Data" tab
5. Select data type from drop down menu: Protocol (for this example)
6. Import data: protocolXbh.dat
7. Go to "Electrodes (XYZ/Topo)" tab and hit "Import from CSV files (no headers)"
8. Import borehole electrode configuration file: elecXbh.csv
9. Invert data and follow the rest of the instructions from pyR2 manual.