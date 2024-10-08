suppressMessages(require(tidyverse))
suppressMessages(require(geepack))
suppressMessages(library(emmeans))
suppressMessages(source('scripts/utils.R'))

# todo read in as arguments
dat_file = 'data/rakai_drug_resistance_categorized.tsv'
w_file = 'models/all_sequencing_probs.tsv'
strata = c('round')
drop_stratum = c()
signif = 0.05


# read in data
hiv_dr_cat = read_tsv(dat_file, show_col_types=FALSE)
# sampling weights
all_weights = read_tsv(w_file, show_col_types=FALSE)
# merge em
hiv_dr_cat = hiv_dr_cat %>% left_join(all_weights, by=c("study_id", "round"))
# get dates for each round
round_dates = hiv_dr_cat %>% select(round, round_mid_date) %>% unique() %>% 
	mutate(round = as.character(round))

# just viremic participants
# just the rounds we want
dat = hiv_dr_cat %>% 
	filter(
		round > 16 & round < 19) %>%
	mutate(round = as.character(round))


stratified_rr = list()
stratified_pred = list()
for (class in c("nnrti", "nrti", "pi")){
	print(class)
	# calculate weights
	# here we are modeling resistant viremia among all participants
	# people not living with HIV: w = 1
	# people living with suppressed HIV: w = 1
	# people living with viremic HIV without sequence data: w = 0
	# people living with viremic HIV with sequence data: w = (1/p_i) / sum_i (1-p_i)
	class_dat = dat %>% 
		mutate(outcome_cat = case_when( 
				is.na(finalhiv) ~ 'none', # todo impute missing finalhiv measures to N 
				finalhiv == 'N' ~ 'none',
				finalhiv == 'P' & !viremic ~ 'none', 
				finalhiv == 'P' & viremic & is.na(!!as.name(class)) ~ 'no_data',
				finalhiv == 'P' & viremic & !is.na(!!as.name(class)) ~ 'data'))
	stopifnot(!any(is.na(class_dat$outcome_cat)))
	class_dat = class_dat %>%
		mutate(w = case_when(
				outcome_cat == 'data' ~ 1 / !!as.name(paste('p_', class, sep='')),
				outcome_cat == 'no_data' ~ 0,
				outcome_cat == 'none' ~ 1),
			r = outcome_cat == 'data' & 
					(!!as.name(class) == 'intermediate' | !!as.name(class) == 'high')) %>%
		group_by(round, outcome_cat != 'none') %>%
		mutate(w = n() * w/sum(w)) %>%
		ungroup() 
	counts = class_dat %>% group_by(round) %>% summarise(w = round(sum(w)), n =n())
	stopifnot(all(counts$w == counts$n))
	# no need to pass these observations to the model
	class_dat = class_dat %>% filter(w > 0)
	# only interested in those with resistance
	class_dat = class_dat %>% 
		filter(r == 1)
	for (s in strata){
		# tabulate outcome n 
		s_n = bind_rows(
			class_dat %>% 
				group_by(across(all_of(c(s)))) %>% 
				summarise(n = n(), y=sum(pre_treatment), .groups='drop') %>%
				mutate(s = !!s))
		for (s_var in unique(class_dat[[s]])){
			s_var_class_dat = class_dat[class_dat[[s]] == s_var,]
			# independent correlation structure because stratified by round
			s_var_output = run_process_model(s_var_class_dat, pre_treatment ~ 1, 
					paste(c('models/tx_among', class, '_', s, "_", s_var), collapse=''),
					corstr='independence')
			stratified_rr[[paste(c(class, s, s_var), collapse='_')]] =
				as_tibble(s_var_output$o, rownames='var') %>% 
					mutate(corr = 'independence') %>%
					mutate(
						class = class,
						s = s, !!s := s_var)
		}
	}
}


stratified_rr = bind_rows(stratified_rr)
write_tsv(stratified_rr, 'models/tx_amongDr_rr_stratified.tsv')

