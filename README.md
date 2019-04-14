# Newton_Rhapson_Load_Flow_Analysis
This repository contains the Matlab code for load flow analysis when the line data and bus data are available in the specified format.

## Format for input data
The input data is provided in the text files *bus_dat.txt* and *line_dat.txt* in the given format.

### Bus Data

| Bus No. | Bus Type | Voltage(p.u.) | Angle(deg) | Pg(MW) | Qg(MVAR) | Pload(MW) | Qload(MVAR) |
|---------|----------|---------------|------------|--------|----------|-----------|-------------|
| 1 | 101 | 1.00 | 0 | 0 | 0 | 0.55 | 0.13 |
| 2 | 101 | 1.00 | 0 | 0 | 0 | 0 | 0 |
| 3 | 101 | 1.00 | 0 | 0 | 0 | 0.30 | 0.18 |
| 4 | 101 | 1.00 | 0 | 0 | 0 | 0.50 | 0.05 |
| 5 | 102 | 1.03 | 0 | 0.75 | 0 | 0.30 | 0.10 |
| 6 | 103 | 1.02 | 0 | 0 | 0 | 0 | 0 |
