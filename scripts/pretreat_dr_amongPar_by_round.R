suppressMessages(require(tidyverse))
suppressMessages(require(geepack))
suppressMessages(library(emmeans))
suppressMessages(source('scripts/utils.R'))

# todo read in as arguments
dat_file = 'data/rakai_drug_resistance_categorized.tsv'
w_file = 'models/all_sequencing_probs.tsv'

# read in data
hiv_dr_cat = read_tsv(dat_file, show_col_types=FALSE)
# sampling weights
all_weights = read_tsv(w_file, show_col_types=FALSE)
# merge em
hiv_dr_cat = hiv_dr_cat %>% left_join(all_weights, by=c("study_id", "round"))

# get round dates
round_dates = hiv_dr_cat %>% select(round, round_mid_date) %>% unique() %>% mutate(round = as.character(round))

dat = hiv_dr_cat %>%
	mutate(round = as.character(round))

# fit a model with round as predictor
base_pred = list()
base_rr = list()
for (class in c("nnrti", "nrti", "pi")){
	print(class)
	# calculate weights
	# here we are modeling viremic pre-treatment resistance among all participants
	# people not living with HIV: w = 1
	# people living with suppressed HIV: w = 1
	# people living with viremic HIV who are treatment experienced: w = 1
	# people living with viremic HIV who are pre-treatment w/ seq data: w = (1/p_i) / sum_i (1-p_i)
	# people living with viremic HIV who are pre-treatment w/o seq data: w = 0 
	class_dat = dat %>% 
		mutate(outcome_cat = case_when( 
					is.na(finalhiv) ~ 'none', # todo impute missing finalhiv measures to N 
					finalhiv == 'N' ~ 'none',
					finalhiv == 'P' & !viremic ~ 'none', 
					finalhiv == 'P' & pre_treatment == FALSE ~ 'none',
					finalhiv == 'P' & viremic & pre_treatment == TRUE & is.na(!!as.name(class)) ~ 'no_data',
					finalhiv == 'P' & viremic & pre_treatment == TRUE & !is.na(!!as.name(class)) ~ 'data'))
	stopifnot(!any(is.na(class_dat$outcome_cat)))
	class_dat = class_dat %>%
		mutate(w = case_when(
				outcome_cat == 'data' ~ 1 / !!as.name(paste('p_', class, sep='')),
				outcome_cat == 'no_data' ~ 0,
				outcome_cat == 'none' ~ 1),
			r = outcome_cat == 'data' & 
				(!!as.name(class) == 'intermediate' | !!as.name(class) == 'high')) %>%
		group_by(round, outcome_cat != 'none') %>%
		mutate(w = n() * w/sum(w)) %>% # normalize weights to population size
		ungroup()
	counts = class_dat %>% group_by(round) %>% summarise(w = round(sum(w)), n =n())
	stopifnot(all(counts$w == counts$n))
	# no need to pass these observations to the model
	class_dat = class_dat %>% filter(w > 0)
	# tabulate outcome n 
	n = class_dat %>% 
			group_by(round) %>% 
			summarise(n = n(), y=sum(r))
	# base model
	base_output = 
		run_process_model(class_dat, r ~ round, 
			paste(c('models/pretreat_', class, '_amongPar'), collapse=''))
	# finally, get average prevalence in each round
	base_pred[[class]] = 
		as_tibble(emmeans(
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
				class = class) %>%
		left_join(n, by='round')
	base_rr[[class]] = as_tibble(base_output$o, rownames='var') %>% 
			mutate(corr = base_output$min_qic_corr) %>%
		mutate(class = class)
}

# combine output and save
base_pred = bind_rows(base_pred)
write_tsv(base_pred, 'models/pretreat_amongPar_prev_pred.tsv')

base_rr = bind_rows(base_rr)
write_tsv(base_rr, 'models/pretreat_amongPar_rr.tsv')








