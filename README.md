# predicting marine mammal density from eDNA

Code repository for: *Microbial and small plankton environmental DNA predicts density of blue, fin, and humpback whales in the southern California bight.*

Citation:

> E.V. Satterthwaite, T.D. Ruiz, K.G. Chan, N. Patrick, M.N. Alksne, N.V. Patin, J. Dinasquet, R.H. Lampe, A.O. Shelton, L. Thomas, B. Semmens. Microbial and small plankton environmental DNA predicts density of blue, fin, and humpback whales in the southern California bight.

Repository contributors: T.D. Ruiz, N. Patrick, K.G. Chan

## Notes

The file `analysis.R` is a high-level script showing the order of execution of individual analyses. 

- The scripts in `~/scripts/processing/` will not execute without the raw data files, which are not included with this repository; the results, however, are stored in `~/data/processed/`.
- The remaining scripts reproduce the analysis in the paper from this processed data; outputs are stored in `~/rslt/`. Note that execution times are long (approximately 3 hours per script).
- The stability selection and nested validation portions of the analysis require data partitions not included in the repository due to file size. These can be found here [ADD LINK]

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
