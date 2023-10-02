# HIV drug resistance in the Rakai Community Cohort Study, 2012 - 2019
This repository contains the data and code needed to replicate the analyses in the manuscript "Population dynamics of HIV drug resistance among pre-treatment and treatment-experienced PLHIV during treatment scale-up in Uganda: a population-based longitudinal study". The raw data is located in the ./data directory and the ./scripts/run_all.sh shell script runs all of the analyses. All code was written by Michael A. Martin, PhD who can be reached at mmart108@jhmi.edu with any questions about the code. Inquiries about the research itself can be directed Dr. M. Kate Grabowski, PhD at mgrabow2@jhmi.edu, who is corresponding author on the manuscript. 

The analysis is all run in the R programming language with the following packages: 

tidyverse v.2.0.0
cowplot v.1.1.1
geepack v.1.3.9
ggpattern v.1.0.1
ggplo2 v.3.4.3
ggrepel v.0.9.3
haven v.2.5.3
mtest v.0.9.40
Readxl v.1.4.3
sandwich v.3.0.2.

A couple of notes: 
1. For privacy reasons, we do not share individual level community data. As such, the data needed to replicate Figure 3A and B and Table S28 have not been shared, although the code is shared here. Requests to use these data can be directed to Dr. M. Kate Grabowski at mgrabow2@jhmi.edu
2. Because many models are stratified by survey round and use "missing vl" as a predictor and missing vl only exists in round 16, many models return a rank-deficiency warning. Because we are doing prediction for missing vl in R15 this does not affect the reliability of the results. 

