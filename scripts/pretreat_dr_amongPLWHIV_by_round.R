library(tidyverse)
library(lmtest)
library(sandwich)
source('scripts/utils.R')


# todo read in as arguments
dat_file = 'data/rakai_drug_resistance_categorized.tsv'
w_file = 'models/pretreat_weights.tsv'
strata = c('sex', 'age_cat', 'comm_type')
drop_stratum = c()
strata_class = c('nnrti')

# read in data
hiv_dr_cat = read_tsv(dat_file)
# sampling weights
pretreat_weights = read_tsv(w_file)
# merge em
hiv_dr_cat = hiv_dr_cat %>% left_join(pretreat_weights, by=c("study_id", "round"))
# get dates for each round
round_dates = hiv_dr_cat %>% select(round, round_mid_date) %>% unique() %>% 
	mutate(round = as.character(round))
# just viremic HIV+ pretreatment individuals
pretreat_dat = hiv_dr_cat %>% 
	filter(
		viremic & pre_treatment & finalhiv == "P")


# fit a model per round
# using a gee model for flexibility later on if needed,
# but independent correlation for now so could as well be glm
# using gee also runs sandwich variance estimators by default
base_pred = list()
stratified_pred = list()
base_rr = list()
stratified_rr = list()
for (class in c("nnrti", "nrti", "pi")){
	print(class)
	class_dat = pretreat_dat %>% 
		rename(
			r = !!class, 
			w = !!paste("w_", class, sep='')) %>%
		filter(!is.na(r)) %>%
		mutate(
			r = viremic & (r == 'intermediate/high'), 
			round = as.character(round)) %>%
		group_by(study_id) %>%
		mutate(idx = cur_group_id()) %>%
		ungroup()
	# base model
	base_output = 
		run_process_model(class_dat %>% arrange(idx, round), r ~ round, 
			paste(c('models/pretreat_', class, '_amongPLWHIV'), collapse=''))
	# finally, get average prevalence in each round
	base_pred[[class]] = as_tibble(predict_glm(
			base_output$m, 
			tibble(round = unique(class_dat$round)),
			~round,
			v=base_output$vcov, type="response")) %>%
		mutate(class = class)
	base_rr[[class]] = as_tibble(base_output[['o']]) %>%
		mutate(var = rownames(base_output[['o']]),
			class = class)
	if (class %in% strata_class) {
		for (s in strata){
			print(s)
			s_class_dat = class_dat %>%
				rename(s_var = !!s) %>%
				filter(!paste(s, s_var, sep=':') %in% drop_stratum)
			s_output = run_process_model(s_class_dat, r ~ round*s_var, 
				paste(c('models/pretreat_', class, '_amongPLWHIV_', s), collapse=''))
			stratified_pred[[paste(class, s, sep='_')]] = 
				as_tibble(predict_glm(
					s_output$m, 
					s_class_dat %>% select(round, s_var) %>% unique(),
					~round*s_var,
					v=s_output$vcov, type="response")) %>% arrange(round, s_var)  %>%
				# collapse anythign > 1 to 1
				mutate(
					s = s,
					class = class,
					upr = if_else(upr > 1, 1, upr))
			stratified_rr[[paste(class, s, sep='_')]] = as_tibble(s_output[['o']]) %>%
				mutate(var = rownames(s_output[['o']]),
					class = class, 
					s = s) 
		}
	}
}

# combine output and save
base_pred = bind_rows(base_pred)
write_tsv(base_pred, 'models/pretreat_amongPLWHIV_prev_pred.tsv')

stratified_pred = bind_rows(stratified_pred)
write_tsv(stratified_pred, 'models/pretreat_amongPLWHIV_prev_pred_stratified.tsv')

base_rr = bind_rows(base_rr)
write_tsv(base_rr, 'models/pretreat_amongPLWHIV_rr.tsv')

stratified_rr = bind_rows(stratified_rr)
write_tsv(stratified_rr, 'models/pretreat_amongPLWHIV_rr_stratified.tsv')









