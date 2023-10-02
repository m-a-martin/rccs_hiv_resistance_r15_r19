library(tidyverse)
source('scripts/utils.R')
library(lmtest)
library(sandwich)


# todo read in as arguments
dat_file = 'data/rakai_drug_resistance_categorized.tsv'

# read in data
hiv_dr_cat = read_tsv(dat_file) %>%
	filter(age_cat != '(50, 100]' &  round > 14) %>%
	mutate(comm_num = as.character(comm_num))

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

# list of all unique community numbers and rounds
comms = unique(prev_dat$comm_num)
rounds = unique(prev_dat$round)

# first overall HIV prevalence
h_prev_pred = list()
for (i_round in unique(prev_dat$round)){
	h_prev_output = 
		run_process_model(prev_dat %>% filter(round == i_round), 
			 h ~ comm_num, 
			paste(c('models/HIV_comm_num_R', i_round), collapse=''))
	h_prev_pred[[i_round]] = as_tibble(predict_glm(
			h_prev_output$m, 
			prev_dat %>% filter(round == i_round) %>% select(comm_num) %>% unique(),
			~comm_num,
			v=h_prev_output$vcov, type="response")) %>% mutate(round = i_round)
}

h_prev_pred = bind_rows(h_prev_pred)
write_tsv(h_prev_pred, 'models/HIV_comm_num_prev_pred.tsv')

# viremic hiv prevalence
v_prev_pred = list()
for (i_round in unique(prev_dat$round)){
	v_prev_output = 
		run_process_model(prev_dat %>% filter(round == i_round), 
			 v ~ comm_num, 
			paste(c('models/HIV_comm_num_viremic_R', i_round), collapse=''))
	v_prev_pred[[i_round]] = as_tibble(predict_glm(
			v_prev_output$m, 
			prev_dat %>% filter(round == i_round) %>% select(comm_num) %>% unique(),
			~comm_num,
			v=v_prev_output$vcov, type="response")) %>% mutate(round = i_round)
}

write_tsv(bind_rows(v_prev_pred), 'models/HIV_comm_num_viremic_prev_pred.tsv')


