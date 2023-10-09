library(tidyverse)
library(lmtest)
library(sandwich)
source('scripts/utils.R')


# todo read in as arguments
dat_file = 'data/rakai_drug_resistance_categorized.tsv'
w_file = 'models/all_weights.tsv'

# read in data
hiv_dr_cat = read_tsv(dat_file)
w = read_tsv(w_file) %>% select(study_id, round, w_all)

# filter
d = hiv_dr_cat %>% 
	filter(round == 18 & finalhiv == 'P' & viremic) %>%
	left_join(w, by = c('study_id', 'round')) %>%
	mutate(
		all = finalhiv == 'P' & viremic & (!is.na(nnrti) & !is.na(nrti) & !is.na(pi) & !is.na(insti)),
		outcome_cat = case_when(
			is.na(finalhiv) ~ 'none',
			finalhiv == 'N' ~ 'none',
			finalhiv == 'P' & !viremic ~ 'none',
			finalhiv == 'P' & !all ~ 'no_data',
			finalhiv == 'P' & all ~ 'data'),
		w = case_when(
			outcome_cat == 'data' ~ w_all,
			outcome_cat == 'no_data' ~ 0,
			outcome_cat == 'none' ~ 1)) %>%
	filter(w > 0)

# define outcome variable
d = d %>% mutate(
	multi = case_when(
		finalhiv == 'P' & viremic ~ 
			(nnrti == 'intermediate/high') & 
				(pi == 'intermediate/high' | nrti == 'intermediate/high'),
		finalhiv == 'N' | !viremic ~ 
			FALSE),
	single = case_when(
		finalhiv == 'P' & viremic ~ 
			(nnrti == 'intermediate/high') & 
				(pi != 'intermediate/high' & nrti != 'intermediate/high'),
		finalhiv == 'N' | !viremic ~ 
			FALSE)) 

multi_m = glm(multi~1, data=d, family=poisson(link="log"), 
			weights=w)

multi_r = coeftest(multi_m, vcov=sandwich)

multi_prev = tibble(
	type = 'multi', 
	fit = exp(multi_r[,1]),
	se.fit = exp(multi_r[,2]),
	lwr = exp(multi_r[,1] + qnorm(0.05/2) * multi_r[,2]),
	upr = exp(multi_r[,1] - qnorm(0.05/2) * multi_r[,2]))


single_m = glm(single~1, data=d, family=poisson(link="log"), 
			weights=w)
single_r = coeftest(single_m, vcov=sandwich)


single_prev = tibble(
	type = 'single', 
	fit = exp(single_r[,1]),
	se.fit = exp(single_r[,2]),
	lwr = exp(single_r[,1] + qnorm(0.05/2) * single_r[,2]),
	upr = exp(single_r[,1] - qnorm(0.05/2) * single_r[,2]))

prev_pred = bind_rows(single_prev, multi_prev)

write_tsv(prev_pred, 'models/all_multi_amongPLWHIV_prev_pred.tsv')



