# To do #
# install multi-plot package that isn't cowplot, forget which it is

#### NOTES ####
# 1. Because many models are stratified by round and use 
#	"missing vl" as a prediction and missing vl only exists 
#	in R15, many models return a rank-deficiency warning
#	but, because we are doing prediction for missing vl in R15
#	this does not affect the reliability of the results 

#### --------------------- ####
#### 0. CREATE DIRECTORIES ####
#### --------------------- ####
mkdir -p figures
mkdir -p figures/final
mkdir -p models
mkdir -p tables

#### -------------------------------- ####
#### 1. FORMAT DR DATA AND CATEGORIZE ####
#### -------------------------------- ####
# categorizing results into susceptible, low, intemrediate, and high resistance
Rscript scripts/categorize_dr_dat.R
# format mutation data
Rscript scripts/format_dr_muts.R

#### ------------------------------------------------------------------------- ####
#### 2. PREVALENCE OF ALL HIV, ALL VIREMIC HIV, AND PRE-TREATMENT HIV BY ROUND ####
#### ------------------------------------------------------------------------- ####
# prevalence of hiv by round and viremic prevalence stratified by community type
Rscript scripts/hiv_prev_by_round.R
# prevalence of hiv stratified by community number
# community-level data not being made public
#Rscript scripts/hiv_prev_by_round_comm_num.R

#### ---------------------------------------- ####
#### 3. CALCULATE INVERSE PROBABILITY WEIGHTS ####
#### ---------------------------------------- ####
Rscript scripts/pretreat_sequencing_weights.R
Rscript scripts/all_sequencing_weights.R
Rscript scripts/treat_sequencing_weights.R

#### ---------------------------------------- ####
#### 4. TEST CORRELATION STRUCTURE ####
#### ---------------------------------------- ####
Rscript scripts/test_correlation_structure.R

#### ----------------------------------------- ####
#### 5. PREVALENCE OF PRE-TREATMENT RESISTANCE #### 
#### ----------------------------------------- ####
## 3A. Analaysis amongst everyone ##
# now model with pretreatment PLWHIV as population
Rscript scripts/pretreat_dr_amongPLWHIV_by_round.R
# and model with all participants as population
Rscript scripts/pretreat_dr_amongPar_by_round.R
## 3B. Incident case analysis ##
# first calculate weights
Rscript scripts/pretreat_inc_sequencing_weights.R
# now run analysis
Rscript scripts/pretreat_dr_amongInc_by_round.R

#### ---------------------------------- ####
#### 6. PREVALENCE OF RESISTANT VIREMIA ####
#### ---------------------------------- ####
# then run analysis with viremic PLWHIV as population
Rscript scripts/all_dr_amongPLWHIV_by_round.R
# per community data not being made public
#Rscript scripts/all_dr_amongPLWHIV_by_comm.R
# now run analysis with participants as population
Rscript scripts/all_dr_amongPar_by_round.R

#### ------------------------------------------------- ####
#### 7. PREVALENCE OF TREATMENT-EXPERIENCED RESISTANCE ####
#### ------------------------------------------------- ####
# now run analysis with viremic treatment-experienced PLWHIV as population
Rscript scripts/treat_dr_amongPLWHIV_by_round.R
# run analysis with participants as population
Rscript scripts/treat_dr_amongPar_by_round.R

#### -------------------------------------- ####
#### 8. PREVALENCE OF MULTICLASS RESISTANCE ####
#### -------------------------------------- ####
Rscript scripts/all_multi_dr_amongPLWHIV.R

#### --------------------------------------------------- ####
#### 9. PREVALENCE OF RESISTANCE STRATIFIED BY COMMUNITY ####
#### --------------------------------------------------- ####
# community-level data not being made public
#Rscript scripts/all_dr_amongPar_by_comm.R

#### ------------------------------------------ ####
#### 10. PREVALENCE OF INDIVIDUAL SUBSTITUTIONS ####
#### ------------------------------------------ ####
mkdir -p models/mutations
# prevalence of substitutions among pre-treatment individuals 
Rscript scripts/pretreat_dr_muts_amongPLWHIV_by_round.R
# k103 stratified by sex
Rscript scripts/pretreat_k103n_amongPLWHIV_by_round_sex.R
# prevalence of substitutions among treatment-experienced individuals
Rscript scripts/treat_dr_muts_amongPLWHIV_by_round.R
# just T97A
Rscript scripts/pretreat_t97a_amongPLWHIV_by_round.R
# test T97A correlation structure
Rscript scripts/pretreat_t97a_amongPLWHIV_by_round_gee.R

#### ---------------- ####
#### 11. MAKE FIGURES ####
#### ---------------- ####
Rscript scripts/figure_1.R
Rscript scripts/figure_2.R
#Rscript scripts/figure_3.R
Rscript scripts/figure_4.R
Rscript scripts/figure_s1.R
Rscript scripts/figure_s2.R
Rscript scripts/figure_s3.R
Rscript scripts/figure_s4.R
Rscript scripts/figure_s5.R
Rscript scripts/figure_s6.R
Rscript scripts/figure_s7.R
Rscript scripts/figure_s8.R

#### ---------------- ####
#### 12. MAKE TABLES ####
#### ---------------- ####
Rscript scripts/table_1.R
Rscript scripts/table_2.R
Rscript scripts/table_3.R
# table S1 data not being made public
#Rscript scripts/table_s1.R
Rscript scripts/table_s2.R
Rscript scripts/table_s3.R
Rscript scripts/table_s5.R
Rscript scripts/table_s6.R
Rscript scripts/table_s7.R
Rscript scripts/table_s8.R
Rscript scripts/table_s9.R
Rscript scripts/table_s10.R
Rscript scripts/table_s11.R
Rscript scripts/table_s12.R
Rscript scripts/table_s13.R
Rscript scripts/table_s14.R
Rscript scripts/table_s15.R
Rscript scripts/table_s16.R
Rscript scripts/table_s17.R
Rscript scripts/table_s18.R
Rscript scripts/table_s19.R
Rscript scripts/table_s20.R
Rscript scripts/table_s21.R
Rscript scripts/table_s22.R
Rscript scripts/table_s23.R
Rscript scripts/table_s24.R
Rscript scripts/table_s25.R
Rscript scripts/table_s26.R
Rscript scripts/table_s27.R
# individual level communityd data not being shared publicly 
#Rscript scripts/table_s28.R
Rscript scripts/table_s29.R
Rscript scripts/table_s30.R
Rscript scripts/table_s31.R
Rscript scripts/table_s32.R
Rscript scripts/table_s33.R

# copy over figures and tables with individual level community data
cp ../tables/table_s28.tsv ./tables/table_s28.tsv
cp ../figures/final/figure_3.pdf ./figures/final/figure_3.pdf

# copy over table that requires full RCCS dataset
cp ../tables/table_s1.tsv ./tables/table_s1.tsv

