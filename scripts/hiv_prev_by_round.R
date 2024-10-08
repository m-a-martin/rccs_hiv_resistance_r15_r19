suppressMessages(library(tidyverse))
suppressMessages(library(geepack))
suppressMessages(library(emmeans))
suppressMessages(source('scripts/utils.R'))

# todo read in as arguments
dat_file = 'data/rakai_drug_resistance_categorized.tsv'

# read in data
hiv_dr_cat = read_tsv(dat_file, show_col_types=FALSE) 

# process dat 
prev_dat = hiv_dr_cat %>% mutate(
	round = as.character(round),
	h = (!is.na(finalhiv) & finalhiv == 'P'),
	v = 
		case_when(
			round == 15 ~ NA,
			round > 15 ~ h & viremic),
	p = h & viremic & pre_treatment,
	t = 
		case_when(
			round == 15 ~ NA,
			round > 15 ~ h & viremic & !pre_treatment),
	w = 1) %>%
	group_by(study_id) %>%
	mutate(idx = cur_group_id()) %>% 
	ungroup() %>%
	arrange(idx, round)

round_dates = hiv_dr_cat %>% select(round, round_mid_date) %>% unique() %>% 
	mutate(round = as.character(round))

out_names = c(
	'h' = 'models/HIV',
	'v' = 'models/HIV_viremic',
	'p' = 'models/HIV_pretreat',
	't' = 'models/HIV_treat')

if (!dir.exists('models')) {dir.create('models')}


for (y in c('h', 'v', 'p', 't')){
	# first overall HIV prevalence
	prev_output = run_process_model(
		prev_dat[!is.na(prev_dat[[y]]),], 
		as.formula(paste(y, '~ round', sep='')), 
		out_names[y])
	# tabulate outcome n 
	n = prev_dat[!is.na(prev_dat[[y]]),] %>% 
			group_by(round) %>% 
			summarise(n = n(), y=sum(!!as.name(y)))
	# predict
	prev_pred = 
		as_tibble(emmeans(
			prev_output$m,  
			~ round,
			type="response")) %>%
		rename(c(
			'fit' = 'rate',
			'se.fit' = 'SE',
			'lwr'= 'lower.CL', 
			'upr'='upper.CL')) %>%
		mutate(corr = prev_output$min_qic_corr) %>%
		left_join(n, by='round')
	write_tsv(prev_pred, 
		paste(out_names[y], '_prev_pred.tsv', sep=''))
	write_tsv(as_tibble(prev_output$o, rownames='var') %>% 
			mutate(corr = prev_output$min_qic_corr), 
		paste(out_names[y], '.tsv', sep=''))
	saveRDS(prev_pred, paste(c(out_names[y], '_prev_pred.rds'), collapse=''))
}


