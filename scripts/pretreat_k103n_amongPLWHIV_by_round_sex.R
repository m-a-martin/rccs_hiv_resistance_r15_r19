library(tidyverse)
source('scripts/utils.R')
library(lmtest)
library(sandwich)

# todo read in as arguments
dr_cat_file = 'data/rakai_drug_resistance_categorized.tsv'
dat_file = 'data/rakai_drug_resistance_mut_formatted.tsv'
w_file = 'models/pretreat_weights.tsv'

hiv_dr_cat = read_tsv(dr_cat_file, show_col_types=FALSE)

# read in data
hiv_muts = read_tsv(dat_file, show_col_types=FALSE)

# just HIV+ pretreatment individuals
pretreat_muts = hiv_muts %>% 
	mutate(round = as.character(round)) %>%
	filter(viremic & pre_treatment & finalhiv == "P") 

# metadata for all pretreatment samples 
# that have resistance prediction for at least one drug
pretreat_dat = hiv_dr_cat %>% 
	mutate(round = as.character(round)) %>%
	filter(
		viremic & pre_treatment & finalhiv == "P" & 
		!is.na(nnrti) & !is.na(nrti) & !is.na(pi) & !is.na(insti))

# sampling weights
pretreat_weights = read_tsv(w_file, show_col_types=FALSE)
# merge em
pretreat_dat = pretreat_dat %>% 
	left_join(
		pretreat_weights %>% 
			mutate(round = as.character(round)) %>%
			select(study_id, round, w_all), by=c("study_id", "round")) %>%
	rename(w = w_all)



stratified_pred = list()
stratified_rr = list()
# fit a model for each mutation
muts = c('rtK103N')
for (i_mut in muts){
	# which samples ahve this mutation
	mut_dat = pretreat_muts %>% 
		filter(mut == i_mut)
	# need to merge in epi data
	# if they're not in the mut dat then they lack this mutation
	mut_dat = mut_dat %>% full_join(pretreat_dat, by=c('round', 'study_id')) %>%
		mutate(
			m = if_else(is.na(mut), FALSE, TRUE))
	base_output = 
		run_process_model(mut_dat, m ~ round*sex, 
			paste(c('models/mutations/pretreat_', i_mut, '_sex_amongPLWHIV'), collapse=''))
	stratified_pred[[i_mut]] = as_tibble(predict_glm(
		base_output[['m']], 
		mut_dat %>% select(round, sex) %>% unique() %>% filter(!is.na(sex)),
		~round*sex, v=base_output[['vcov']], type="response")) %>%
		mutate(mut = i_mut, s = 'sex')
	stratified_rr[[i_mut]] = tibble(
		mut = i_mut,
		s = 'sex',
		var = rownames(base_output[['o']]),
		RR = base_output[['o']][,1],
		LCI = base_output[['o']][,2],
		UCI = base_output[['o']][,3],
		P = base_output[['o']][,4]
		)
}

stratified_pred = bind_rows(stratified_pred)
write_tsv(stratified_pred, 'models/pretreat_dr_muts_amongPLWHIV_prev_pred_stratified.tsv')

stratified_rr = bind_rows(stratified_rr)
write_tsv(stratified_rr, 'models/pretreat_dr_muts_amongPLWHIV_rr_pred_stratified.tsv')
