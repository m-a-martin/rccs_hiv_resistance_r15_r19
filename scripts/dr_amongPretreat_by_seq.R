suppressMessages(require(tidyverse))
suppressMessages(require(geepack))
suppressMessages(library(emmeans))
suppressMessages(source('scripts/utils.R'))

# todo read in as arguments
dat_file = 'data/rakai_drug_resistance_categorized.tsv'
w_file = 'models/all_sequencing_probs.tsv'
strata = c('sequenced')
drop_stratum = c()
min_stratum_obs = 1

# read in data
hiv_dr_cat = read_tsv(dat_file, show_col_types=FALSE) %>% filter(round == 16)
# sampling weights
all_weights = read_tsv(w_file, show_col_types=FALSE)
# merge em
hiv_dr_cat = hiv_dr_cat %>% left_join(all_weights, by=c("study_id", "round"))
# just viremic HIV+ pretreatment individuals
pretreat_dat = hiv_dr_cat %>% 
	filter(
		viremic & pre_treatment & finalhiv == "P") %>%
	mutate(round = as.character(round))

# only one round so fit with sequencing approach as predictor
base_pred = list()
base_rr = list()
for (class in c("nnrti", "nrti", "pi")){
	print(class)
	class_dat = pretreat_dat %>% 
		mutate(outcome_cat = case_when( 
				finalhiv == 'P' & viremic & is.na(!!as.name(class)) ~ 'no_data',
				finalhiv == 'P' & viremic & !is.na(!!as.name(class)) ~ 'data'))
	stopifnot(!any(is.na(class_dat$outcome_cat)))
	class_dat = class_dat %>%
		mutate(w = case_when(
				outcome_cat == 'data' ~ 1 / !!as.name(paste('p_', class, sep='')),
				outcome_cat == 'no_data' ~ 0),
			r = outcome_cat == 'data' & 
					(!!as.name(class) == 'intermediate' | !!as.name(class) == 'high')) %>%
		group_by(round, outcome_cat != 'none') %>%
		mutate(w = n() * w/sum(w)) %>%
		ungroup()
	counts = class_dat %>% group_by(round) %>% summarise(w = round(sum(w)), n =n())
	stopifnot(all(counts$w == counts$n))
	# no need to pass these observations to the model
	class_dat = class_dat %>% filter(w > 0)
	# tabulate number of outcomes in each category
	n_outcomes = class_dat %>% 
		group_by(sequenced) %>% 
		summarise(n_outcome_obs = sum(r))
	# base model
	base_output = 
		run_process_model(class_dat, r ~ sequenced, 
			paste(c('models/dr_', class, '_amongPretreat'), collapse=''))
	# finally, get average prevalence in each round
	base_pred[[class]] = 
		as_tibble(emmeans(
			base_output$m,  
			~ sequenced,
			type="response")) %>% 
		rename(c(
			'fit' = 'rate',
			'se.fit' = 'SE',
			'lwr'= 'lower.CL', 
			'upr'='upper.CL')) %>%
		mutate(
			corr = base_output$min_qic_corr,
			class = class) %>%
		left_join(n_outcomes, by='sequenced')
	base_rr[[class]] = as_tibble(base_output$o, rownames='var') %>% 
			mutate(
				corr = base_output$min_qic_corr,
				class = class)
}

# combine output and save
base_pred = bind_rows(base_pred)
write_tsv(base_pred, 'models/dr_amongPretreat_seq_prev_pred.tsv')

base_rr = bind_rows(base_rr)
write_tsv(base_rr, 'models/dr_amongPretreat_seq_rr.tsv')







