library(tidyverse)
source('scripts/utils.R')
library(lmtest)
library(sandwich)


# todo read in as arguments
dat_file = 'data/rakai_drug_resistance_categorized.tsv'
w_file = 'models/all_weights.tsv'
strata = c('comm_type')
drop_stratum = c()
# read in data
hiv_dr_cat = read_tsv(dat_file)
# sampling weights
all_weights = read_tsv(w_file)
# merge em
hiv_dr_cat = hiv_dr_cat %>% left_join(all_weights, by=c("study_id", "round"))
# get dates for each round
round_dates = hiv_dr_cat %>% select(round, round_mid_date) %>% unique() %>% 
	mutate(round = as.character(round))

# all participants!
# drop weird age categories
# just the rounds we want
dat = hiv_dr_cat %>% 
	filter(
		age_cat != '(50, 100]' & 
		round > 16 & round < 19)

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
	class_dat = dat %>% 
		rename(
			r = !!class, 
			w = !!paste("w_", class, sep='')) %>%
		mutate(
			round = as.character(round)) %>%
		group_by(study_id) %>%
		mutate(idx = cur_group_id()) %>%
		ungroup() %>%
		mutate(outcome_cat = case_when(
				is.na(finalhiv) ~ 'none',
				finalhiv == 'N' ~ 'none',
				finalhiv == 'P' & viremic & is.na(r) ~ 'no_data',
				finalhiv == 'P' & viremic & !is.na(r) ~ 'data'),
			w = case_when(
				outcome_cat == 'data' ~ w,
				outcome_cat == 'no_data' ~ 0,
				outcome_cat == 'none' ~ 1),
			r = (finalhiv == 'P') & viremic & (!is.na(r) & (r == 'intermediate/high')))
	# base model
	base_output = 
		run_process_model(class_dat %>% arrange(idx, round), r ~ round, 
			paste(c('models/all_', class, '_amongPar'), collapse=''))
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
	for (s in strata){
			print(s)
			s_class_dat = class_dat %>%
				rename(s_var = !!s) %>%
				filter(!paste(s, s_var, sep=':') %in% drop_stratum)
			s_output = run_process_model(s_class_dat, r ~ round*s_var, 
				paste(c('models/all_', class, '_amongPar_', s), collapse=''))
			stratified_pred[[paste(class, s, sep='_')]] = 
				as_tibble(predict_glm(
					s_output$m, 
					s_class_dat %>% select(round, s_var) %>% unique(),
					~round*s_var,
					v=s_output$vcov, type="response")) %>% arrange(round, s_var) %>%
					# collapse anythign > 1 to 1
					mutate(
						class = class,
						s = s,
						upr = if_else(upr > 1, 1, upr)) 
			stratified_rr[[paste(class, s, sep='_')]] = as_tibble(s_output[['o']]) %>%
				mutate(var = rownames(s_output[['o']]),
					class = class,
					s = s)
	}
}

# combine output and save
base_pred = bind_rows(base_pred)
write_tsv(base_pred, 'models/all_amongPar_prev_pred.tsv')

base_rr = bind_rows(base_rr)
write_tsv(base_rr, 'models/all_amongPar_rr.tsv')

stratified_pred = bind_rows(stratified_pred)
write_tsv(stratified_pred, 'models/all_amongPar_prev_pred_stratified.tsv')


stratified_rr = bind_rows(stratified_rr)
write_tsv(stratified_rr, 'models/all_amongPar_rr_stratified.tsv')







