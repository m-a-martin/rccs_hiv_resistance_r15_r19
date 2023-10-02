library(tidyverse)
source('scripts/utils.R')
library(lmtest)
library(sandwich)

# todo read in as arguments
dat_file = 'data/rakai_drug_resistance_categorized.tsv'
w_file = 'models/all_weights.tsv'

# read in data
hiv_dr_cat = read_tsv(dat_file)
# sampling weights
all_weights = read_tsv(w_file)
# merge em
hiv_dr_cat = hiv_dr_cat %>% left_join(all_weights, by=c("study_id", "round"))

i_round = 18
# all participants!
# just the rounds we want
dat = hiv_dr_cat %>% 
	filter(
		finalhiv == "P" & 
		viremic &
		age_cat != '(50, 100]' & 
		round > 16 & round < 19) %>%
	mutate(comm_num = as.character(comm_num))

comms = unique(dat$comm_num)
rounds = unique(as.character(dat$round))
# fit a model per round
base_pred = list()
for (class in c("nnrti", "nrti", "pi")){
	print(class)
	class_dat = dat %>% 
		rename(
			r = !!class, 
			w = !!paste("w_", class, sep='')) %>%
		filter(!is.na(r)) %>%
		mutate(
			r = viremic & (r == 'intermediate/high'), 
			round = as.character(round)) %>%
		group_by(study_id) %>%
		mutate(idx = cur_group_id()) %>%
		ungroup()
	base_output = 
			run_process_model(class_dat, r ~ comm_num, 
				paste(c('models/all_', class, '_amongPar_comm_num_R', i_round), collapse=''))
	# finally, get average prevalence in each round
	base_pred[[class]] = 
		as_tibble(predict_glm(
			base_output$m, 
			class_dat %>% select(comm_num) %>% unique(),
			~comm_num,
			v=base_output$vcov, type="response")) %>%
		mutate(class = class, round = i_round)
}


# combine output and save
base_pred = bind_rows(base_pred)
write_tsv(base_pred, 'models/all_amongPLWHIV_comm_num_prev_pred.tsv')




