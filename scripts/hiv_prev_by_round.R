library(tidyverse)
source('scripts/utils.R')
library(sandwich)
library(lmtest)
library(sandwich)

# todo read in as arguments
dat_file = 'data/rakai_drug_resistance_categorized.tsv'

# read in data
hiv_dr_cat = read_tsv(dat_file) 

# process dat 
prev_dat = hiv_dr_cat %>% mutate(
	round = as.character(round),
	h = (!is.na(finalhiv) & finalhiv == 'P'),
	v = h & viremic,
	p = h & viremic & pre_treatment,
	t = h & viremic & !pre_treatment,
	w = 1) %>%
	group_by(study_id) %>%
	mutate(idx = cur_group_id()) %>% 
	ungroup() %>%
	arrange(idx, round)

round_dates = hiv_dr_cat %>% select(round, round_mid_date) %>% unique() %>% 
	mutate(round = as.character(round))

# first overall HIV prevalence
h_prev_output = run_process_model(prev_dat, h ~ round, 
		paste(c('models/HIV'), collapse=''))

h_prev_pred = as_tibble(predict_glm(
		h_prev_output$m, 
		tibble(round = unique(prev_dat$round)),
		~round,
		v=h_prev_output$vcov, type="response"))

write_tsv(h_prev_pred, 'models/HIV_prev_pred.tsv')
write_tsv(as_tibble(h_prev_output$o_format), 'models/HIV_format.tsv')

# viremic hiv prevalence
v_dat = prev_dat %>% filter(round >= 16)
v_prev_output = 
	run_process_model(v_dat, v ~ round, 
		paste(c('models/HIV_viremic'), collapse=''))

v_prev_pred = as_tibble(predict_glm(
		v_prev_output$m, 
		tibble(round = unique(v_dat$round)),
		~round,
		v=v_prev_output$vcov, type="response"))

write_tsv(v_prev_pred, 'models/HIV_viremic_prev_pred.tsv')
write_tsv(as_tibble(v_prev_output$o_format), 'models/HIV_viremic_format.tsv')

# viremic hiv prevalence stratified by community type
v_prev_output_stratified = 
	run_process_model(v_dat, v ~ round*comm_type, 
		paste(c('models/HIV_viremic'), collapse=''))

v_prev_pred_stratified = as_tibble(predict_glm(
		v_prev_output_stratified$m, 
		v_dat %>% select(round, comm_type) %>% unique(),
		~round*comm_type,
		v=v_prev_output_stratified$vcov, type="response"))

write_tsv(v_prev_pred_stratified, 'models/HIV_viremic_prev_pred_stratified.tsv')


# pretreatment viremic hiv prevalence
p_prev_output = 
	run_process_model(prev_dat, p ~ round, 
		paste(c('models/HIV_pretreat'), collapse=''))

p_prev_pred = as_tibble(predict_glm(
		p_prev_output$m, 
		tibble(round = unique(prev_dat$round)),
		~round,
		v=p_prev_output$vcov, type="response"))

write_tsv(p_prev_pred, 'models/HIV_pretreat_viremic_prev_pred.tsv')
write_tsv(as_tibble(p_prev_output$o_format), 'models/HIV_pretreat_format.tsv')


# treatment experienced viremic hiv prevalence
t_dat = prev_dat %>% filter(round >= 16)
t_prev_output = 
	run_process_model(t_dat, t ~ round, 
		paste(c('models/HIV_treat'), collapse=''))

t_prev_pred = as_tibble(predict_glm(
		t_prev_output$m, 
		tibble(round = unique(t_dat$round)),
		~round,
		v=t_prev_output$vcov, type="response"))

write_tsv(t_prev_pred, 'models/HIV_treat_viremic_prev_pred.tsv')
write_tsv(as_tibble(t_prev_output$o_format), 'models/HIV_treat_format.tsv')


