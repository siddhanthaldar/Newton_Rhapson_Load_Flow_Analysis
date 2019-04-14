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

**101 : P-Q Bus;    102 : P-|V| Bus;    103 : |V|-$\theta$ Bus**

Number of buses :  
P-Q Bus = 4  
P-|V| Bus = 1  
|V|-$\theta$ Bus = 1  

### Line Data

| From Bus | To Bus | Resistance r(p.u.) | Reactance x(p.u.) | Line Charging B(p.u.) | ONR |
|---------|----------|---------------|------------|--------|----------|
| 6 | 2  | 0.080 | 0.370 | 0.280 | 1.000 |
| 6 | 4 | 0.123 | 0.518 | 0.400 | 1.000 |
| 5 | 1 | 0.723 | 1.050 | 0.200 | 1.000 |
| 5 | 3 | 0.282 | 0.640 | 0.300 | 1.000 |
| 2 | 4 | 0.097 | 0.407 | 0.240 | 1.000 |
| 2 | 1 | 0.0 | 0.133 | 0.000 | 1.000 |
| 4 | 3 | 0.0 | 0.300 | 0.000 | 1.000 |


