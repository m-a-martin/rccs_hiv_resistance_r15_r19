suppressMessages(require(tidyverse))
suppressMessages(require(geepack))
suppressMessages(library(emmeans))
suppressMessages(source('scripts/utils.R'))

# todo read in as arguments
dat_file = 'data/rakai_drug_resistance_categorized.tsv'
w_file = 'models/all_sequencing_probs.tsv'

strata = c('comm_type', 'age_cat', 'sex')
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

# all participants!
# just the rounds we want
dat = hiv_dr_cat %>% 
	filter(
		round > 16 & round < 19) %>%
	mutate(round = as.character(round))

# fit a model with round as predictor
base_pred = list()
stratified_pred = list()
base_rr = list()
stratified_rr = list()
bivar_coeffs = list()

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
	# tabulate outcome n 
	n = class_dat %>% 
			group_by(round) %>% 
			summarise(n = n(), y=sum(r))
	# base model
	base_output = 
		run_process_model(class_dat, r ~ round, 
			paste(c('models/dr_', class, '_amongPar'), collapse=''))
	# finally, get average prevalence in each round
	base_pred[[class]] = 
		as_tibble(emmeans(
			base_output$m,  
			~ round,
			type="response")) %>% 
				mutate(corr = base_output$min_qic_corr) %>%
			rename(c(
				'fit' = 'rate',
				'se.fit' = 'SE',
				'lwr'= 'lower.CL', 
				'upr'='upper.CL')) %>%
			mutate(class = class)  %>%
		left_join(n, by='round')
	base_rr[[class]] = as_tibble(base_output$o, rownames='var') %>% 
			mutate(corr = base_output$min_qic_corr) %>%
		mutate(class = class)
	for (s in strata){
		print(s)
		s_bi_output = run_process_model(class_dat, as.formula(paste('r~round + ', s, sep='')), 
						paste(c('models/dr_', class, '_amongPar_', s), collapse=''),
						corstr=base_output$min_qic_corr)
		bivar_coeffs[[paste(c(class, s), collapse='_')]] =
				as_tibble(s_bi_output$o, rownames='var') %>% 
					mutate(corr = base_output$min_qic_corr) %>%
					mutate(
						class = class,
						s = s)
		# if any significant do stratified analysis
		if (any((bivar_coeffs[[paste(c(class, s), collapse='_')]] %>% 
				filter(var != "(Intercept)" & !grepl('round', var)))$P < signif)){
			print('significant')
			# tabulate outcome n 
			s_n = bind_rows(
				class_dat %>% 
					group_by(across(all_of(c('round', s)))) %>% 
					summarise(n = n(), y=sum(r), .groups='drop') %>%
					mutate(s = !!s),
				class_dat %>% 
					group_by(across(all_of(c('round')))) %>% 
					summarise(n = n(), y=sum(r), .groups='drop') %>%
					mutate(s = !!s),
				class_dat %>% 
					group_by(across(all_of(c(s)))) %>% 
					summarise(n = n(), y=sum(r), .groups='drop') %>%
					mutate(s = !!s))
			for (s_var in unique(class_dat[[s]])){
				s_var_class_dat = class_dat[class_dat[[s]] == s_var,]
				s_var_output = run_process_model(s_var_class_dat, r ~ round, 
						paste(c('models/dr_', class, '_amongPar_', s, "_", s_var), collapse=''),
						corstr=base_output$min_qic_corr)
				stratified_pred[[paste(c(class, s, s_var), collapse='_')]] = 
					as_tibble(emmeans(
							s_var_output$m,  
							~ round,
							type="response")) %>% 
						mutate(corr = base_output$min_qic_corr) %>%
						rename(c(
							'fit' = 'rate',
							'se.fit' = 'SE',
							'lwr'= 'lower.CL', 
							'upr'='upper.CL')) %>%
						mutate(s = s, !!s := s_var, class = class) %>%
						left_join(s_n, by=c('round', 's', s))
				stratified_rr[[paste(c(class, s, s_var), collapse='_')]] =
					as_tibble(s_var_output$o, rownames='var') %>% 
						mutate(corr = base_output$min_qic_corr) %>%
						mutate(
							class = class,
							s = s, !!s := s_var)
			}
		}
	}
}

# combine output and save
base_pred = bind_rows(base_pred)
write_tsv(base_pred, 'models/dr_amongPar_prev_pred.tsv')

base_rr = bind_rows(base_rr)
write_tsv(base_rr, 'models/dr_amongPar_rr.tsv')

stratified_pred = bind_rows(stratified_pred)
write_tsv(stratified_pred, 'models/dr_amongPar_prev_pred_stratified.tsv')

stratified_rr = bind_rows(stratified_rr)
write_tsv(stratified_rr, 'models/dr_amongPar_rr_stratified.tsv')

write_tsv(bind_rows(bivar_coeffs), 'models/dr_amongPar_bivar_coeffs.tsv')
