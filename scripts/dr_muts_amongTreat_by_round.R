suppressMessages(require(tidyverse))
suppressMessages(require(geepack))
suppressMessages(library(emmeans))
suppressMessages(source('scripts/utils.R'))

# todo read in as arguments
dr_cat_file = 'data/rakai_drug_resistance_categorized.tsv'
dat_file = 'data/rakai_drug_resistance_mut_formatted.tsv'
w_file = 'models/all_sequencing_probs.tsv'

# read in mutation data
hiv_muts = read_tsv(dat_file, show_col_types=FALSE)
# just HIV+ treatment individuals
treat_muts = hiv_muts %>% 
	mutate(round = as.character(round)) %>%
	filter(viremic & !pre_treatment & finalhiv == "P" & round > 16 & round < 19) %>%
	group_by(mut) %>%
	filter(n() >= 5) %>%
	ungroup()

# read in data
hiv_dr_cat = read_tsv(dr_cat_file, show_col_types=FALSE)
# sampling weights
all_weights = read_tsv(w_file, show_col_types=FALSE)
# merge em
hiv_dr_cat = hiv_dr_cat %>% left_join(all_weights, by=c("study_id", "round"))
# just viremic HIV+ treatment individuals
treat_dat = hiv_dr_cat %>% 
	filter(
		viremic & !pre_treatment & finalhiv == "P" & round > 16 & round < 19) %>%
	mutate(
		round = as.character(round),
		outcome_cat = case_when( 
			finalhiv == 'P' & viremic & 
				(is.na(insti) | is.na(nnrti) | is.na(nrti) | is.na(pi)) ~ 'no_data',
			finalhiv == 'P' & viremic & 
				!is.na(insti) & !is.na(nnrti) & !is.na(nrti) & !is.na(pi) ~ 'data'))

stopifnot(!any(is.na(treat_dat$outcome_cat)))

# calculate weights
treat_dat = treat_dat %>%
		mutate(w = case_when(
				outcome_cat == 'data' ~ 1 / p_insti_nnrti_nrti_pi,
				outcome_cat == 'no_data' ~ 0,
				outcome_cat == 'none' ~ 1)) %>%
		group_by(round, outcome_cat != 'none') %>%
		mutate(w = n() * w/sum(w)) %>%
		ungroup()

# no need to pass these observations to the model
treat_dat = treat_dat %>% filter(w > 0)

if (!dir.exists('models/mutations')){dir.create('models/mutations')}

all_rr = list()
all_pred = list()
stratified_pred = list()
# fit a model for each mutation
for (i_mut in (treat_muts %>% group_by(mut) %>% summarise(n=n()) %>% arrange(-n))$mut){
	print(i_mut)
	# which samples have this mutation
	mut_dat = treat_muts %>% 
		filter(mut == i_mut)
	# need to merge in epi data
	# if they're not in the mut dat then they lack this mutation
	mut_dat = mut_dat %>% right_join(treat_dat, by=c('round', 'study_id')) %>%
		mutate(
			m = if_else(is.na(mut), FALSE, TRUE))
	stopifnot(all(!is.na(mut_dat$w)))
	# tabulate number of outcomes in each category
	n = mut_dat %>% 
			group_by(round) %>% 
			summarise(n = n(), y=sum(m))
	base_output = 
		run_process_model(mut_dat, m ~ round, 
			paste(c('models/mutations/', i_mut, '_amongTreat'), collapse=''), corstr='independence')
	all_pred[[i_mut]] = as_tibble(emmeans(
			base_output$m,  
			~ round,
			type="response")) %>% 
		rename(c(
			'fit' = 'rate',
			'se.fit' = 'SE',
			'lwr'= 'lower.CL', 
			'upr'='upper.CL')) %>%
		mutate(
			class = i_mut,
			corr = 'independence') %>%
		left_join(n, by='round')
	all_rr[[i_mut]] = as_tibble(base_output$o, rownames='var') %>% 
			mutate(
				class = i_mut,
				corr = 'independence')
}

all_pred = bind_rows(all_pred) 
write_tsv(all_pred, 'models/dr_muts_amongTreat_prev_pred.tsv')

all_rr = bind_rows(all_rr)
write_tsv(all_rr, 'models/dr_muts_amongTreat_rr_pred.tsv')
