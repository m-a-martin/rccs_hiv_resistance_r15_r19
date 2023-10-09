library(tidyverse)
source('scripts/utils.R')
library(lmtest)
library(sandwich)
library(geepack)


# todo read in as arguments
dr_cat_file = 'data/rakai_drug_resistance_categorized.tsv'
dat_file = 'data/rakai_drug_resistance_mut_formatted.tsv'
w_file = 'models/pretreat_weights.tsv'
i_mut = 'inT97A'
# read in data
hiv_muts = read_tsv(dat_file)

# just HIV+ pretreatment individuals
pretreat_muts = hiv_muts %>% 
	mutate(round = as.character(round)) %>%
	filter(pre_treatment & finalhiv == "P") 


# metadata for all pretreatment samples that were sequenced
hiv_dr_cat = read_tsv(dr_cat_file)
pretreat_dat = hiv_dr_cat %>% 
	mutate(round = as.character(round)) %>%
	filter(
		pre_treatment & finalhiv == "P" & viremic & 
		valid_dr_dat)

# sampling weights
pretreat_weights = read_tsv(w_file)
# merge em
pretreat_dat = pretreat_dat %>% 
	left_join(
		pretreat_weights %>% 
			mutate(round = as.character(round)) %>%
			select(study_id, round, w_all), by=c("study_id", "round")) %>%
	rename(w = w_all)


mut_dat = pretreat_muts %>% 
		filter(mut == i_mut)
# need to merge in epi data
# if they're not in the mut dat then they lack this mutation
mut_dat = mut_dat %>% full_join(pretreat_dat, by=c('round', 'study_id')) %>%
	mutate(
		m = if_else(is.na(mut), FALSE, TRUE)) %>%
	select(study_id, round, m, w) %>%
	arrange(study_id, round) %>% 
	group_by(study_id) %>%
	mutate(idx = cur_group_id()) %>% 
	ungroup()

corr_pred = list()
corr_rr = list()
for (corr in c("independence", "exchangeable", "ar1", "unstructured")){

	corr_output = geeglm(m ~ round, data=mut_dat,
	 	id = idx, family=poisson(link="log"), corstr=corr, weights=w)
	corr_rr[[corr]] = tibble(
		c = corr,
		RR = exp(corr_output$coefficients[2:length(corr_output$coefficients)]),
		LCI = exp(corr_output$coefficients[2:length(corr_output$coefficients)] + 
				qnorm(0.05/2)*
					summary(corr_output)$coefficients[,2][2:length(corr_output$coefficients)]), 
		UCI = 
			exp(corr_output$coefficients[2:length(corr_output$coefficients)] - 
				qnorm(0.05/2)*
					summary(corr_output)$coefficients[,2][2:length(corr_output$coefficients)]), 
		P = summary(corr_output)$coefficients[2:length(corr_output$coefficients),4],
		var = names(corr_output$coefficients[2:length(corr_output$coefficients)]))
	corr_pred[[corr]] =  as_tibble(predict_glm(
		corr_output, 
		mut_dat %>% select(round) %>% unique(),
		~round, v=corr_output$geese$vbeta, type="response")) %>% mutate(c = corr)
}

corr_pred = bind_rows(corr_pred)
write_tsv(corr_pred, 'models/pretreat_int97a_amongPLWHIV_gee_prev_pred.tsv')

write_tsv(bind_rows(corr_rr), 'models/pretreat_int97a_amongPLWHIV_gee_rr.tsv')


