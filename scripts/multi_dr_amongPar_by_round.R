suppressMessages(require(tidyverse))
suppressMessages(require(geepack))
suppressMessages(library(emmeans))
suppressMessages(source('scripts/utils.R'))


# todo read in as arguments
dat_file = 'data/rakai_drug_resistance_categorized.tsv'
w_file = 'models/all_sequencing_probs.tsv'

# read in data
hiv_dr_cat = read_tsv(dat_file, show_col_types=FALSE) %>%
	filter(round > 16 & round < 19)

# sampling weights
all_weights = read_tsv(w_file, show_col_types=FALSE)

multi_resistance_profiles = c(
	'TRUE_TRUE_TRUE' = 'nnrti_nrti_pi',
	'TRUE_TRUE_FALSE' = 'nnrti_nrti',
	'TRUE_FALSE_TRUE' = 'nnrti_pi',
	'TRUE_FALSE_FALSE' = 'nnrti',
	'FALSE_TRUE_TRUE' = 'nrti_pi',
	'FALSE_TRUE_FALSE' = 'nrti',
	'FALSE_FALSE_TRUE' = 'pi',
	'FALSE_FALSE_FALSE' = 'none')

dat = hiv_dr_cat %>%
	mutate(
		all_gt = !is.na(nnrti) & !is.na(nrti) & !is.na(pi),
		nnrti = if_else(nnrti == 'intermediate' | nnrti == 'high', TRUE, FALSE),
		nrti = if_else(nrti == 'intermediate' | nrti == 'high', TRUE, FALSE),
		pi = if_else(pi == 'intermediate' | pi == 'high', TRUE, FALSE)) %>%
	unite("resistance_profile", nnrti:pi, remove=FALSE) %>%
	mutate(multi = multi_resistance_profiles[resistance_profile])
	
# merge in sampling probabilities
# calculate weights, same weigths for all outcomes
# here we are modeling resistant viremia among all participants
# people not living with HIV: w = 1
# people living with suppressed HIV: w = 1
# people living with viremic HIV without sequence data: w = 0
# people living with viremic HIV with sequence data: w = (1/p_i) / sum_i (1-p_i)
dat = dat %>% 
	left_join(all_weights, by = c('study_id', 'round')) %>%
	mutate(
		round = as.character(round),
		outcome_cat = case_when( 
			is.na(finalhiv) ~ 'none', # todo impute missing finalhiv measures to N 
			finalhiv == 'N' ~ 'none',
			finalhiv == 'P' & !viremic ~ 'none', 
			finalhiv == 'P' & viremic & is.na(multi) ~ 'no_data',
			finalhiv == 'P' & viremic & !is.na(multi) ~ 'data'))

stopifnot(!any(is.na(dat$outcome_cat)))
dat = dat %>%
		mutate(w = case_when(
				outcome_cat == 'data' ~ 1 / p_nnrti_nrti_pi,
				outcome_cat == 'no_data' ~ 0,
				outcome_cat == 'none' ~ 1)) %>%
		group_by(round, outcome_cat != 'none') %>%
		mutate(w = n() * w/sum(w)) %>%
		ungroup()
counts = dat %>% group_by(round) %>% summarise(w = round(sum(w)), n =n())
stopifnot(all(counts$w == counts$n))

base_pred = list()
base_rr = list()
for (class in multi_resistance_profiles[multi_resistance_profiles != "none"]){
	# define outcome variable for each class
	class_dat = dat %>%
		mutate(
			r = outcome_cat == 'data' & 
					(multi == !!class))
	# tabulate outcome n 
	n = class_dat %>% 
			group_by(round) %>% 
			summarise(n = n(), y=sum(r))
	# errors due to 0 obs in R17
	base_output = 
			run_process_model(class_dat, r ~ round, 
				paste(c('models/multi_', class, '_amongPar'), collapse=''))
	# finally, get average prevalence in each round
	base_pred[[class]] = as_tibble(emmeans(
			base_output$m,  
			~ round,
			type="response")) %>% 
		rename(c(
			'fit' = 'rate',
			'se.fit' = 'SE',
			'lwr'= 'lower.CL', 
			'upr'='upper.CL')) %>%
		mutate(
			corr = base_output$min_qic_corr,
			class = class)  %>%
		left_join(n, by='round')

	base_rr[[class]] = as_tibble(base_output$o, rownames='var') %>% 
			mutate(corr = base_output$min_qic_corr, class=class)
		#	 %>%
		#left_join(as_tibble(colSums(model.matrix(n_outcome_obs ~ round, 
		#				data=n_outcomes) * 
		#					n_outcomes$n_outcome_obs), 
		#				rownames='var'), 
		#	by='var')
}

# combine output and save
base_pred = bind_rows(base_pred)
write_tsv(base_pred, 'models/multi_amongPar_prev_pred.tsv')

base_rr = bind_rows(base_rr)
write_tsv(base_rr, 'models/multi_amongPar_rr.tsv')



