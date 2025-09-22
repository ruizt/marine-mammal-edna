# predicting marine mammal density from eDNA

Code repository for: *Microbial and small zooplankton communities predict density of baleen whales in the southern California Current Ecosystem.*

Citation:

> E.V. Satterthwaite, T.D. Ruiz, N.V. Patin, M.N. Alksne, L. Thomas, J. Dinasquet, R.H. Lampe, K.G. Chan, N.A. Patrick, A.E. Allen, S. Baumann-Pickering, B.X. Semmens. Microbial and small zooplankton communities predict density of baleen whales in the southern California Current Ecosystem.

Repository contributors: T.D. Ruiz, N. Patrick, K.G. Chan

## Notes

The file `~/scripts/analysis.R` is a high-level script showing the order of execution of individual steps in the analysis described in the paper: data processing, aggregation, and transformation; stability selection; nested validation; model fitting; figure and table generation. Please note:

- The scripts in `~/scripts/processing/` will not execute without the raw data files, which are not included with this repository. Raw data files are not necessary to reproduce the analysis, but may be requested from the authors along with an explanation of the purpose of the request.
- The results of the scripts in `~/scripts/processing/` needed for further analysis **are** available in this repository and are stored in `~/data/processed/`, with the exception of the data partitions, which can be obtained from [10.5281/zenodo.15678927](https://doi.org/10.5281/zenodo.15678927).
- To reproduce the analysis in the paper starting from processed data, unzip the files `_combined-partitions.zip` and `_partitions.zip` in the `~/data/processed/` directory. Then execute the scripts as shown in `~/scripts/analysis.R` (not including processing scripts). Outputs of analyses are stored in `~/rslt/`. Note that execution times are long (approximately 3 hours per script).

The map figures utilize shapefiles from GSHHG Release v2.3.7, URL [https://www.soest.hawaii.edu/pwessel/gshhg/](https://www.soest.hawaii.edu/pwessel/gshhg/) accessed December 2024.
- Wessel, P., and W. H. F. Smith (1996), A global, self-consistent, hierarchical, high-resolution shoreline database, J. Geophys. Res., 101(B4), 8741–8743, doi:10.1029/96JB00104.

## Software

RStudio Version 2024.12.1+563 (2024.12.1+563)

R version 4.4.3 (2025-02-28)

Platform: x86_64-apple-darwin20

Running under: macOS Sequoia 15.3.2

Packages (not including dependencies):

- sf 1.0-18        
- mapdata 2.3.1
- maps 3.4.2       
- ggmap 4.0.0     
- patchwork 1.2.0  
- openxlsx 4.2.5.2 
- readxl 1.4.3      
- lubridate 1.9.3
- tidyverse 2.0.0 
- collapse 2.0.15 
- pls 2.8-3
- spls 2.2-3 
- modelr 0.1.11  
- fs 1.6.4
- magrittr 2.0.3 
- vegan 2.6-6.1
