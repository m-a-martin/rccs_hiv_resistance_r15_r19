library(tidyverse)
source('scripts/utils.R')
library(lmtest)
library(sandwich)

# todo read in as arguments
dat_file = 'data/rakai_drug_resistance_categorized.tsv'
w_file = 'models/pretreat_inc_weights.tsv'
strata = c('sex', 'age_cat', 'comm_type', 'subtype')
drop_stratum = c('subtype:B')
# read in data
hiv_dr_cat = read_tsv(dat_file)
pretreat_dat = hiv_dr_cat %>% 
	filter(
		viremic &
		incident & 
		pre_treatment & 
		finalhiv == "P" &
		round == first_pos) %>%
	mutate(nnrti_dat = !is.na(nnrti),
		missing_vl = vl_cat == 'missing',
		log10vl = replace_na(log10vl, 0),
		round = as.character(round))

print(paste(c('there are ', nrow(pretreat_dat), ' sero-incident PLWHIV'), collapse=''))
	
# now merge with weights
# one weight per INDIVIDUAL
w = read_tsv(w_file)
pretreat_dat = pretreat_dat %>% left_join(w, by="study_id")

# fit a model per round
# using a gee model for flexibility later on if needed,
# but independent correlation for now so could as well be glm
# using gee also runs sandwich variance estimators by default
base_pred = list()
base_rr = list()
for (class in c("nnrti")){
	class_dat = pretreat_dat %>% 
		rename(
			r = !!class, 
			w = !!paste("w_", class, sep='')) %>%
		mutate(
			r = viremic & (r == 'intermediate/high')) %>%
		# remove those with no data for this drug class
		filter(!is.na(r)) %>%
		# first sample per person
		group_by(study_id) %>%
		filter(round == min(round)) %>%
		ungroup()
	print(paste(c('there are ', nrow(class_dat), ' sero-incident PLWHIV with ', class, ' resistance data'), collapse=''))
	# base model
	base_output = 
		run_process_model(class_dat, r ~ round, 
			paste(c('models/pretreat_', class, '_amongInc'), collapse=''))
	rounds = class_dat %>% select(round) %>% arrange(round) %>% unique()
	pred = as_tibble(predict_glm(base_output$m, rounds,  
		~round, v=base_output[['vcov']], type="response")) %>%
		mutate(class = class)
	# finally, get average prevalence in each round
	base_pred[[class]] = pred
	base_rr[[class]] = as_tibble(base_output[['o']]) %>%
		mutate(var = rownames(base_output[['o']]),
			class = class)
}

base_pred = bind_rows(base_pred)
write_tsv(base_pred, 'models/pretreat_amongInc_pred.tsv')

base_rr = bind_rows(base_rr)
write_tsv(base_rr, 'models/pretreat_amongInc_rr.tsv')









