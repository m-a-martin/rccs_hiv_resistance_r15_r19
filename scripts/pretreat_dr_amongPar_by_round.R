library(tidyverse)
source('scripts/utils.R')
library(lmtest)
library(sandwich)

# todo read in as arguments
dat_file = 'data/rakai_drug_resistance_categorized.tsv'
w_file = 'models/pretreat_weights.tsv'

# read in data
hiv_dr_cat = read_tsv(dat_file)
# sampling weights
pretreat_weights = read_tsv(w_file)
# merge em
hiv_dr_cat = hiv_dr_cat %>% left_join(pretreat_weights, by=c("study_id", "round"))

# get round dates
round_dates = hiv_dr_cat %>% select(round, round_mid_date) %>% unique() %>% mutate(round = as.character(round))
# just filter age categories and rounds we don't want
# no further subsetting because the population here is all participants 
# no further
dat = hiv_dr_cat

# fit a model per round
# using a gee model for flexibility later on if needed,
# but independent correlation for now so could as well be glm
# using gee also runs sandwich variance estimators by default
base_pred = list()
base_rr = list()
for (class in c("nnrti", "nrti", "pi")){
	print(class)
	# individuals with no probability of the outcome
	# regardless of resistance profiling get a weight of 1
	# individuals with a chance of the outcome dependent on resistance profiling
	# get their weight as defined by the sequencing weight if they have data
	# and a weight of 0 otherwise
	# in other words, among pre-treatment HIV+ individuals we upweight
	# the measurements from individuals with resistance profiling to make up 
	# for the fact that we do not have resistance profiling on all individuals 
	class_dat = dat %>% 
		rename(
			r = !!class, 
			w = !!paste("w_", class, sep='')) %>%
		mutate(
			round = as.character(round)) %>%
		mutate(
			outcome_cat = case_when(
				is.na(finalhiv) ~ 'none',
				finalhiv == 'N' ~ 'none',
				finalhiv == 'P' & !pre_treatment ~ 'none',
				finalhiv == 'P' & pre_treatment & !viremic ~ 'none',
				finalhiv == 'P' & pre_treatment & viremic & is.na(r) ~ 'no_data',
				finalhiv == 'P' & pre_treatment & viremic & !is.na(r) ~ 'data'),
			w = case_when(
				outcome_cat == 'data' ~ w,
				outcome_cat == 'no_data' ~ 0,
				outcome_cat == 'none' ~ 1),
			r = 1 * (!is.na(r) & viremic & pre_treatment & (r == 'intermediate/high')))
	# base model
	base_output = 
		run_process_model(class_dat, r ~ round, 
			paste(c('models/pretreat_', class, '_amongPar'), collapse=''))

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


# combine output and save
base_pred = bind_rows(base_pred)
write_tsv(base_pred, 'models/pretreat_amongPar_prev_pred.tsv')

base_rr = bind_rows(base_rr)
write_tsv(base_rr, 'models/pretreat_amongPar_rr.tsv')








