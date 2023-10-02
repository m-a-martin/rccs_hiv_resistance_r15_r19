library(tidyverse)
source('scripts/utils.R')
library(lmtest)
library(sandwich)


# todo read in as arguments
dat_file = 'data/rakai_drug_resistance_categorized.tsv'
w_file = 'models/treat_weights.tsv'

# read in data
hiv_dr_cat = read_tsv(dat_file)
# sampling weights
treat_weights = read_tsv(w_file)
# merge em
hiv_dr_cat = hiv_dr_cat %>% left_join(treat_weights, by=c("study_id", "round"))

# just viremic HIV+ treated individuals
treat_dat = hiv_dr_cat %>% 
	filter(
		!pre_treatment & 
		finalhiv == "P" & 
		viremic &
		round > 16 & 
		round < 19)

# fit a model per round
# using a gee model for flexibility later on if needed,
# but independent correlation for now so could as well be glm
# using gee also runs sandwich variance estimators by default
base_pred = list()
stratified_pred = list()
base_rr = list()
for (class in c("nnrti", "nrti", "pi")){
	print(class)
	class_dat = treat_dat %>% 
		rename(
			r = !!class, 
			w = !!paste("w_", class, sep='')) %>%
		filter(!is.na(r)) %>%
		mutate(
			r = (finalhiv == "P") & viremic & (r == 'intermediate/high'), 
			round = as.character(round)) %>%
		group_by(study_id) %>%
		mutate(idx = cur_group_id()) %>%
		ungroup()
	# base model
	base_output = 
		run_process_model(class_dat %>% arrange(idx, round), r ~ round, 
			paste(c('models/treat_', class, '_amongPLWHIV'), collapse=''))
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

}

# bind rows and save
# combine output and save
base_pred = bind_rows(base_pred)
write_tsv(base_pred, 'models/treat_amongPLWHIV_prev_pred.tsv')

base_rr = bind_rows(base_rr)
write_tsv(base_rr, 'models/treat_amongPLWHIV_rr.tsv')
