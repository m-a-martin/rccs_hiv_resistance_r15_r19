library(tidyverse)
source('scripts/utils.R')
library(lmtest)
library(sandwich)

# todo read in as arguments
dr_cat_file = 'data/rakai_drug_resistance_categorized.tsv'
dat_file = 'data/rakai_drug_resistance_mut_formatted.tsv'
w_file = 'models/treat_weights.tsv'

hiv_dr_cat = read_tsv(dr_cat_file)

# read in data
hiv_muts = read_tsv(dat_file)

# just HIV+ viremic treated individuals
treat_muts = hiv_muts %>% 
	mutate(round = as.character(round)) %>%
	filter(
		!pre_treatment & viremic & finalhiv == "P" & 
		round > 16 & round < 19)

# metadata for all pretreatment samples that were sequenced
treat_dat = hiv_dr_cat %>%
	mutate(round = as.character(round)) %>%
	filter(!pre_treatment & viremic & finalhiv == 'P' & 
		round > 16 & round < 19 & age_cat != '(50, 100]' & dr_dat)

# sampling weights
treat_weights = read_tsv(w_file)

# merge em
treat_dat = treat_dat %>% left_join(
	treat_weights %>% 
		mutate(round = as.character(round)) %>%
		select(study_id, round, w_gen), by=c("study_id", "round")) %>%
	rename(w = w_gen)

all_rr = list()
all_pred = list()
# fit a model for each mutation
for (i_mut in unique(treat_muts$mut)){
	mut_dat = treat_muts %>% 
		filter(mut == i_mut)
	# need to merge in epi data
	# if they're not in the mut dat then they lack this mutation
	mut_dat = mut_dat %>% full_join(treat_dat, by=c('round', 'study_id')) %>%
		mutate(
			m = if_else(is.na(mut), FALSE, TRUE))
	base_output = 
		run_process_model(mut_dat, m ~ round, 
			paste(c('models/mutations/treat_', i_mut, '_amongPLWHIV'), collapse=''))
	all_pred[[i_mut]] = as_tibble(predict_glm(
		base_output[['m']], 
		mut_dat %>% select(round) %>% unique(),
		~round, v=base_output[['vcov']], type="response")) %>%
		mutate(mut = i_mut)
	all_rr[[i_mut]] = tibble(
		mut = i_mut,
		var = rownames(base_output[['o']]),
		rate_ratio = base_output[['o']][,1],
		lwr = base_output[['o']][,2],
		upr = base_output[['o']][,3])
}

all_pred = bind_rows(all_pred)
write_tsv(all_pred, 'models/treat_dr_muts_amongPLWHIV_prev_pred.tsv')


all_rr = bind_rows(all_rr)
write_tsv(all_rr, 'models/treat_dr_muts_amongPLWHIV_rr_pred.tsv')



