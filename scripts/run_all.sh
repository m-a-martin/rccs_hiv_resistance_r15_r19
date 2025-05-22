#### -------------------------------- ####
#### 1. FORMAT DR DATA AND CATEGORIZE ####
#### -------------------------------- ####
# assign 
# format data how we want it, including linking rakai ids and incorporating metadata
Rscript scripts/format_dr_dat.R
# categorizing results into susceptible, low, intemrediate, and high resistance
Rscript scripts/categorize_dr_dat.R
# format mutation data
Rscript scripts/format_dr_muts.R
# first calculate weights
Rscript scripts/calc_sequencing_probs.R


#### ------------- ####
#### 2. RUN MODELS ####
#### ------------- ####
# 2.1A. Prev of HIV, viremic HIV, viremic pre-treatment HIV, viremic treatment-experienced HIV among participants by round
Rscript scripts/hiv_prev_by_round.R

# 2.1B. Prev of viremic resistance among participants by round
Rscript scripts/dr_amongPar_by_round.R

# 2.1C. Prev of pre-treatment viremic resistance among participants by round
Rscript scripts/pretreat_dr_amongPar_by_round.R

# 2.1D. Prev of treatment-experienced viremic resistance among participants by round
Rscript scripts/treat_dr_amongPar_by_round.R

# 2.1E. Prev. of pre-treatment among resistant by round
Rscript scripts/pretreat_amongDr_by_round.R

# 2.2A. Prev. of viremic multi-class resistance among participants by round
Rscript scripts/multi_dr_amongPar_by_round.R

# 2.3A. Prev. of resistance among pre-treatment viremic participants by round
Rscript scripts/dr_amongPretreat_by_round.R

# 2.3B. Prev. of resistance mutations among pre-treatment viremic participants by round
Rscript scripts/dr_muts_amongPretreat_by_round.R

# 2.4A. Prev. of resistance among treatment-experienced viremic participants by round
Rscript scripts/dr_amongTreat_by_round.R

# 2.4B. Prev. of resistance mutations among treatment-experienced viremic participants by round
Rscript scripts/dr_muts_amongTreat_by_round.R

# 2.5 Organize
mkdir -p models/rds
mv models/*.rds models/rds/.

#### --------------- ####
#### 3. PLOT FIGURES ####
#### --------------- ####
# 3.1. Prev of HIV, viremic HIV, viremic pre-treatment HIV, viremic treatment-experienced HIV among participants by round
Rscript scripts/plot_prev_among_par.R
# 3.2.  Prev. of viremic multi-class resistance among participants by round
Rscript scripts/plot_multiclass_prev_among_par.R
# 3.3. Prev. of resistance among pre-treatment viremic participants by round
Rscript scripts/plot_prev_among_pretreat.R
# 3.4. Prev. of resistance among treatment-experienced viremic participants by round
Rscript scripts/plot_prev_among_treat.R
# 3.5 Within host frequency of resistance mutations
Rscript scripts/plot_mut_freq.R
# 3.6 repeat visits
Rscript scripts/plot_repeat_visits.R
# 3.7 Sensitivity to imputation of viremia status
Rscript scripts/plot_pretreat_viremia_sensitivity.R


#### ------------------ ####
#### 4. GENERATE TABLES ####
#### ------------------ ####
# 4.1. Table 1
Rscript scripts/table_1.R

# 4.2. IPW tables
# IPW tables
for r in 15 16 17 18 19; do
Rscript scripts/tabulate_regression_coeffs.R \
	--models \
		models/seq_probs/all_R${r}_nnrti_probs.tsv \
		models/seq_probs/all_R${r}_nrti_probs.tsv \
		models/seq_probs/all_R${r}_pi_probs.tsv \
		models/seq_probs/all_R${r}_nnrti_nrti_pi_probs.tsv \
		models/seq_probs/all_R${r}_insti_nnrti_nrti_pi_probs.tsv \
	--out seq_prob_coeffs_R${r}
done

# 4.3. Demographics of study pop over time
Rscript scripts/tabulate_time_stratified_count.R \
	--dat data/rakai_drug_resistance_categorized.tsv \
	--tableVars ageyrs age_cat comm_type sex \
	--out par_round_stratified_count

Rscript scripts/tabulate_time_stratified_count.R \
	--dat data/rakai_drug_resistance_categorized.tsv \
	--tableVars ageyrs age_cat comm_type sex \
	--out plhiv_round_stratified_count \
	--filter "finalhiv == 'P'"

# 4.4. Treatment-experience over time
Rscript scripts/tabulate_stratified_outcome.R \
	--dat data/rakai_drug_resistance_categorized.tsv \
	--tableVars round \
	--outcome 'pre_treatment == FALSE' \
	--outcomeLabel treatment-experienced \
	--tableCats \
		'finalhiv == "P"' \
	--tableCatLabels PLHIV \
	--out n_plhiv_treatment

# 4.5. Prevalence of HIV, HIV viremia, pre-treatment HIV viremia, and treatment-experienced HIV viremia by round
cat \
	<(awk '{print $0"\tclass"}' <(head -n 1 models/HIV.tsv)) \
	<(awk -F'\t' '{print $0"\tPLHIV"}' <(tail -n +2 models/HIV.tsv)) \
	<(awk -F'\t' '{print $0"\tvPLHIV"}' <(tail -n +2 models/HIV_viremic.tsv)) \
	<(awk -F'\t' '{print $0"\tptvPLHIV"}' <(tail -n +2 models/HIV_pretreat.tsv)) \
	<(awk -F'\t' '{print $0"\ttvPLHIV"}' <(tail -n +2 models/HIV_treat.tsv)) \
	> models/tmp_rr.tsv
cat \
	<(awk '{print $0"\tclass"}' <(head -n 1 models/HIV_prev_pred.tsv)) \
	<(awk -F'\t' '{print $0"\tPLHIV"}' <(tail -n +2 models/HIV_prev_pred.tsv)) \
	<(awk -F'\t' '{print $0"\tvPLHIV"}' <(tail -n +2 models/HIV_viremic_prev_pred.tsv)) \
	<(awk -F'\t' '{print $0"\tptvPLHIV"}' <(tail -n +2 models/HIV_pretreat_prev_pred.tsv)) \
	<(awk -F'\t' '{print $0"\ttvPLHIV"}' <(tail -n +2 models/HIV_treat_prev_pred.tsv)) \
	> models/tmp_prev_pred.tsv
Rscript scripts/tabulate_round_prev_rr.R \
	--prev models/tmp_prev_pred.tsv \
	--rr models/tmp_rr.tsv \
	--classOrder PLHIV vPLHIV ptvPLHIV tvPLHIV \
	--nOut one \
	--out hiv_round_prev_rr
rm -rf models/tmp_rr.tsv
rm -rf models/tmp_prev_pred.tsv

# 4.6. Samples on whcih sequencing was attempted
Rscript scripts/tabulate_stratified_outcome.R \
	--dat data/rakai_drug_resistance_categorized.tsv \
	--tableVars age_cat comm_type sex round vl_cat \
	--outcome '!is.na(sequenced)' \
	--outcomeLabel sequenced \
	--tableCats \
		'finalhiv == "P" & viremic'\
		'finalhiv == "P" & viremic & pre_treatment'\
		'finalhiv == "P" & viremic & !pre_treatment' \
	--tableCatLabels viremic pt_viremic t_viremic \
	--out n_drm_attempted_stratified

# 4.7. Sequencing technology used to generate data
Rscript scripts/tabulate_time_stratified_count.R \
	--dat data/rakai_drug_resistance_categorized.tsv \
	--tableVars sequenced \
	--out n_drm_attempted_approach \
	--filter 'finalhiv == "P" & viremic == TRUE & !is.na(sequenced)'

# 4.8. Sequencing success 
# 4.8.1. >= 1 drug
Rscript scripts/tabulate_stratified_outcome.R \
	--dat data/rakai_drug_resistance_categorized.tsv \
	--tableVars age_cat comm_type sex round vl_cat sequenced\
	--outcome 'valid_dr_dat == TRUE' \
	--outcomeLabel sequenced \
	--filter '!is.na(sequenced)' \
	--tableCats \
		'finalhiv == "P" & viremic'\
		'finalhiv == "P" & viremic & pre_treatment'\
		'finalhiv == "P" & viremic & !pre_treatment' \
	--tableCatLabels viremic pt_viremic t_viremic \
	--out n_drm_any_stratified

# 4.8.2. all INSTIs
Rscript scripts/tabulate_stratified_outcome.R \
	--dat data/rakai_drug_resistance_categorized.tsv \
	--tableVars age_cat comm_type sex round vl_cat sequenced \
	--outcome '!is.na(insti)' \
	--outcomeLabel sequenced \
	--filter '!is.na(sequenced)' \
	--tableCats \
		'finalhiv == "P" & viremic'\
		'finalhiv == "P" & viremic & pre_treatment'\
		'finalhiv == "P" & viremic & !pre_treatment' \
	--tableCatLabels viremic pt_viremic t_viremic \
	--out n_drm_insti_stratified

# 4.8.3. all NNRTIs
Rscript scripts/tabulate_stratified_outcome.R \
	--dat data/rakai_drug_resistance_categorized.tsv \
	--tableVars age_cat comm_type sex round vl_cat sequenced\
	--outcome '!is.na(nnrti)' \
	--outcomeLabel sequenced \
	--filter '!is.na(sequenced)' \
	--tableCats \
		'finalhiv == "P" & viremic'\
		'finalhiv == "P" & viremic & pre_treatment'\
		'finalhiv == "P" & viremic & !pre_treatment' \
	--tableCatLabels viremic pt_viremic t_viremic \
	--out n_drm_nnrti_stratified

# 4.8.4. all NRTIs
Rscript scripts/tabulate_stratified_outcome.R \
	--dat data/rakai_drug_resistance_categorized.tsv \
	--tableVars age_cat comm_type sex round vl_cat sequenced \
	--outcome '!is.na(nrti)' \
	--outcomeLabel sequenced \
	--filter '!is.na(sequenced)' \
	--tableCats \
		'finalhiv == "P" & viremic'\
		'finalhiv == "P" & viremic & pre_treatment'\
		'finalhiv == "P" & viremic & !pre_treatment' \
	--tableCatLabels viremic pt_viremic t_viremic \
	--out n_drm_nrti_stratified

# 4.8.5.  all PIs
Rscript scripts/tabulate_stratified_outcome.R \
	--dat data/rakai_drug_resistance_categorized.tsv \
	--tableVars age_cat comm_type sex round vl_cat sequenced\
	--outcome '!is.na(pi)' \
	--outcomeLabel sequenced \
	--filter '!is.na(sequenced)' \
	--tableCats \
		'finalhiv == "P" & viremic'\
		'finalhiv == "P" & viremic & pre_treatment'\
		'finalhiv == "P" & viremic & !pre_treatment' \
	--tableCatLabels viremic pt_viremic t_viremic \
	--out n_drm_pi_stratified

# 4.8.6. All NNRTIs, NRTIs, and PIs
Rscript scripts/tabulate_stratified_outcome.R \
	--dat data/rakai_drug_resistance_categorized.tsv \
	--tableVars age_cat comm_type sex round vl_cat sequenced \
	--outcome '!is.na(nnrti) & !is.na(nrti) & !is.na(pi)' \
	--outcomeLabel sequenced \
	--filter '!is.na(sequenced)' \
	--tableCats \
		'finalhiv == "P" & viremic'\
		'finalhiv == "P" & viremic & pre_treatment'\
		'finalhiv == "P" & viremic & !pre_treatment' \
	--tableCatLabels viremic pt_viremic t_viremic \
	--out n_seq_nnrti_nrti_pi_stratified

# 4.8.7. All INSTIs, NNRTIs, NRTIs, and PIs
Rscript scripts/tabulate_stratified_outcome.R \
	--dat data/rakai_drug_resistance_categorized.tsv \
	--tableVars age_cat comm_type sex round vl_cat sequenced \
	--outcome '!is.na(insti) & !is.na(nnrti) & !is.na(nrti) & !is.na(pi)' \
	--outcomeLabel sequenced \
	--filter '!is.na(sequenced)' \
	--tableCats \
		'finalhiv == "P" & viremic'\
		'finalhiv == "P" & viremic & pre_treatment'\
		'finalhiv == "P" & viremic & !pre_treatment' \
	--tableCatLabels viremic pt_viremic t_viremic \
	--out n_seq_insti_nnrti_nrti_pi_stratified

# 4.9. Number of samples with resistance
Rscript scripts/tabulate_time_stratified_dr_count.R \
	--dat data/rakai_drug_resistance_categorized.tsv  \
	--out n_genotyped_resistant

# 4.10. Prevalence of viremic resistance among participants by round
# 4.10.1 non-stratified
Rscript scripts/tabulate_round_prev_rr.R \
	--prev models/dr_amongPar_prev_pred.tsv \
	--rr models/dr_amongPar_rr.tsv \
	--classOrder nnrti nrti pi \
	--out dr_amongPar_prev_rr

# 4.10.2 bivariate coefficients
Rscript scripts/tabulate_bivar_coeffs.R \
	--coeffs models/dr_amongPar_bivar_coeffs.tsv \
	--classOrder nnrti nrti pi \
	--out dr_amongPar_bivar_coeffs

# 4.10.3 stratified 
Rscript scripts/tabulate_stratified_round_prev_rr.R \
	--prev models/dr_amongPar_prev_pred_stratified.tsv \
	--rr models/dr_amongPar_rr_stratified.tsv \
	--classOrder nnrti nrti pi \
	--out dr_amongPar_prev_rr_stratified

# 4.11. Prevalence of pre-treatment viremic resistance among participants by round
Rscript scripts/tabulate_round_prev_rr.R \
	--prev models/pretreat_amongPar_prev_pred.tsv \
	--rr models/pretreat_amongPar_rr.tsv \
	--classOrder nnrti nrti pi \
	--out pretreat_amongPar_prev_rr

# 4.12. Prevalence of treatment-experienced viremic resistance among participants by round
Rscript scripts/tabulate_round_prev_rr.R \
	--prev models/treat_amongPar_prev_pred.tsv \
	--rr models/treat_amongPar_rr.tsv \
	--classOrder nnrti nrti pi \
	--out treat_amongPar_prev_rr

# 4.13. Prevalence of multi-class viremic resistance among participants by round
Rscript scripts/tabulate_round_prev_rr.R \
	--prev models/multi_amongPar_prev_pred.tsv \
	--rr models/multi_amongPar_rr.tsv \
	--classOrder nnrti nrti pi nnrti_nrti nnrti_pi nrti_pi nnrti_nrti_pi \
	--nOut one \
	--out multi_amongPar_prev_rr

# 4.14. Prevalence of resistance among pre-treatment viremic participants by round
# 4.14.1 non-stratified
Rscript scripts/tabulate_round_prev_rr.R \
	--prev models/dr_amongPretreat_prev_pred.tsv \
	--rr models/dr_amongPretreat_rr.tsv \
	--classOrder nnrti nrti pi \
	--out dr_amongPretreat_prev_rr
# 4.14.2 bivariate coefficients
Rscript scripts/tabulate_bivar_coeffs.R \
	--coeffs models/dr_amongPretreat_bivar_coeffs.tsv \
	--classOrder nnrti nrti pi \
	--out dr_amongPretreat_bivar_coeffs 
# 4.14.3 stratified
Rscript scripts/tabulate_stratified_round_prev_rr.R \
	--prev models/dr_amongPretreat_prev_pred_stratified.tsv \
	--rr models/dr_amongPretreat_rr_stratified.tsv \
	--classOrder nnrti nrti pi \
	--out dr_amongPretreat_prev_rr_stratified
# 4.14.3 Prevalence of resistance mutations
Rscript scripts/tabulate_round_prev_rr.R \
	--prev models/dr_muts_amongPretreat_prev_pred.tsv \
	--rr models/dr_muts_amongPretreat_rr_pred.tsv \
	--out dr_muts_amongPretreat_prev_rr \
	--nOut one \
	--endRound 18

# 4.15 Prevalence of resistance among treatment-experienced viremic participants by round
# 4.15.1 non-stratified
Rscript scripts/tabulate_round_prev_rr.R \
	--prev models/dr_amongTreat_prev_pred.tsv \
	--rr models/dr_amongTreat_rr.tsv \
	--classOrder nnrti nrti pi \
	--out dr_amongTreat_prev_rr
# 4.15.2 bivariate
Rscript scripts/tabulate_bivar_coeffs.R \
	--coeffs models/dr_amongTreat_bivar_coeffs.tsv \
	--classOrder nnrti nrti pi \
	--out dr_amongTreat_bivar_coeffs 
# 4.15.2 stratified
Rscript scripts/tabulate_stratified_round_prev_rr.R \
	--prev models/dr_amongTreat_prev_pred_stratified.tsv \
	--rr models/dr_amongTreat_rr_stratified.tsv \
	--classOrder nnrti nrti pi \
	--out dr_amongTreat_prev_rr_stratified
# 4.15.3 Prevalence of resistance mutations
Rscript scripts/tabulate_round_prev_rr.R \
	--prev models/dr_muts_amongTreat_prev_pred.tsv \
	--rr models/dr_muts_amongTreat_rr_pred.tsv \
	--nOut one \
	--out dr_muts_amongTreat_prev_rr \
	--endRound 18

# 4.6. Count of observed resistance mutations
Rscript scripts/tabulate_n_muts.R

