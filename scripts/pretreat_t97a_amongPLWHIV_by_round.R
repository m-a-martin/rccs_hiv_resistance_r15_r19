library(tidyverse)
source('scripts/utils.R')
library(lmtest)
library(sandwich)


# todo read in as arguments
dr_cat_file = 'data/rakai_drug_resistance_categorized.tsv'
dat_file = 'data/rakai_drug_resistance_mut_formatted.tsv'
w_file = 'models/pretreat_weights.tsv'
i_mut = 'inT97A'
strata = c('sex', 'age_cat', 'comm_type', 'subtype')
drop_stratum = c('subtype:B', 'subtype:C', 'subtype:other (non-recombinant)', 'subtype:missing')
# read in data
hiv_muts = read_tsv(dat_file)

# just HIV+ pretreatment individuals
pretreat_muts = hiv_muts %>% 
	mutate(round = as.character(round)) %>%
	filter(pre_treatment & finalhiv == "P" & viremic) 


# metadata for all pretreatment samples that were sequenced
hiv_dr_cat = read_tsv(dr_cat_file)
pretreat_dat = hiv_dr_cat %>% 
	mutate(round = as.character(round)) %>%
	filter(
		pre_treatment & finalhiv == "P" & viremic & 
		valid_dr_dat)

# sampling weights
pretreat_weights = read_tsv(w_file)
# merge em
pretreat_dat = pretreat_dat %>% 
	left_join(
		pretreat_weights %>% 
			mutate(round = as.character(round)) %>%
			select(study_id, round, w_all), by=c("study_id", "round")) %>%
	rename(w = w_all)


mut_dat = pretreat_muts %>% 
		filter(mut == i_mut)
# need to merge in epi data
# if they're not in the mut dat then they lack this mutation
mut_dat = mut_dat %>% full_join(pretreat_dat, by=c('round', 'study_id')) %>%
	mutate(
		m = if_else(is.na(mut), FALSE, TRUE))



stratified_pred = list()
stratified_rr = list()
for (s in strata){
	s_dat = mut_dat %>%
			rename(s_var = !!s) %>%
			filter(!paste(s, s_var, sep=':') %in% drop_stratum) 
	s_output = run_process_model(s_dat, m ~ round*s_var, 
			paste(c('models/pretreat_', i_mut, '_amongPLWHIV_', s), collapse=''))
	stratified_pred[[s]] = as_tibble(predict_glm(
			s_output$m, 
			s_dat %>% select(round, s_var) %>% unique(),
			~round*s_var,
			v=s_output$vcov, type="response")) %>%
			mutate(s = s, s_var = s_var,
				upr=if_else(upr > 1, 1, upr)) 
	stratified_rr[[s]] = as_tibble(s_output[['o']]) %>%
			mutate(var = rownames(s_output[['o']]),
				s = s) 
}

stratified_pred = bind_rows(stratified_pred)
write_tsv(stratified_pred, 'models/pretreat_int97a_amongPLWHIV_prev_pred.tsv')

write_tsv(bind_rows(stratified_rr), 'models/pretreat_int97a_amongPLWHIV_rr.tsv')


