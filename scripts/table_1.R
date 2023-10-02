library(tidyverse)

source('scripts/utils.R')



format_col = function(x1, x2){
		return(paste(paste(paste(x1, " (", sep=''),
			round((x1 / x2)*100, 2), sep=''),
			'%)', sep=''))
}


dat_file = 'data/rakai_drug_resistance_categorized.tsv'
table_vars = c("age_cat", "comm_type", "sex")

# read in data and filter out bad rounds and age categories
hiv_dr_cat = read_tsv(dat_file) %>% 
	mutate(round = as.character(round))

hiv_dr = read_tsv('data/rakai_drug_resistance_formatted.tsv')

t1_c1 = 
	bind_rows(
		hiv_dr_cat %>% summarise(n=n()) %>% 
			mutate(round = 'overall') %>% 
			select(round, n),
		tibble(round = 'Survey round', n = NA), 
		hiv_dr_cat %>% 
			group_by(round) %>% 
			summarise(n=n())) %>%
		rename(all_participant_visits=n)

t1_c2 = 
	bind_rows(
		hiv_dr_cat %>% filter(finalhiv == 'P') %>%
			summarise(n=n()) %>%
			mutate(round = "overall") %>%
			select(round, n),
		tibble(round = "Survey round", n=NA),
		hiv_dr_cat %>% group_by(round) %>% 
			filter(finalhiv == 'P') %>%
			summarise(n=n())) %>%
	rename(n_hiv_participant_visits=n)

t1_c3 = 
	bind_rows(
		hiv_dr_cat %>% filter(finalhiv == 'P' & viremic) %>%
			summarise(n=n()) %>%
			mutate(round = "overall") %>%
			select(round, n),
		tibble(round = "Survey round", n=NA),
		hiv_dr_cat %>% group_by(round) %>% 
			filter(finalhiv == 'P' & viremic) %>%
			summarise(n=n())) %>%
	rename(n_viremic_participant_visits=n) 

t1_c4 = 
	bind_rows(
		hiv_dr_cat %>% 
			filter(finalhiv == 'P' & viremic & pre_treatment)  %>%
			summarise(n=n()) %>%
			mutate(round = "overall") %>%
			select(round, n),
		tibble(round = "Survey round", n=NA),
		hiv_dr_cat %>% 
			group_by(round) %>% 
			filter(finalhiv == 'P' & viremic & pre_treatment) %>%
			summarise(n=n())) %>% 
	rename(n_pretreat_participant_visits=n)


t1_c5 = 
	bind_rows(
		hiv_dr_cat %>% 
			filter(finalhiv == 'P' & viremic & !pre_treatment)  %>%
			summarise(n=n()) %>%
			mutate(round = "overall") %>%
			select(round, n),
		tibble(round = "Survey round", n=NA),
		hiv_dr_cat %>% 
			group_by(round) %>% 
			filter(finalhiv == 'P' & viremic & pre_treatment) %>%
			summarise(n=n())) %>% 
	rename(n_treat_participant_visits=n)



#t1 = t1_c1 %>% 
#	left_join(t1_c2, by='round') %>% 
#	left_join(t1_c3, by='round') %>%
#	left_join(t1_c4, by="round") %>%
#	mutate_at(vars(hiv_participant_visits, viremic_participant_visits, pretreat_participant_visits), 
#		~if_else(!is.na(.), format_col(., all_participant_visits), '')) %>%
#	mutate(
#		all_participant_visits = 
#			if_else(!is.na(all_participant_visits), 
#				as.character(all_participant_visits), ''),
#		viremic_participant_visits = if_else(round == 15, '..', as.character(viremic_participant_visits)))


t1 = t1_c1 %>% 
	left_join(t1_c2, by='round') %>% 
	left_join(t1_c3, by='round') %>%
	left_join(t1_c4, by="round") %>%
	mutate(
		hiv_participant_visits = if_else(!is.na(n_hiv_participant_visits), format_col(n_hiv_participant_visits, all_participant_visits), ''),
		viremic_participant_visits = if_else(!is.na(n_viremic_participant_visits), format_col(n_viremic_participant_visits, n_hiv_participant_visits), ''),
		pretreat_participant_visits = if_else(!is.na(n_pretreat_participant_visits), format_col(n_pretreat_participant_visits, n_viremic_participant_visits), ''),
		all_participant_visits = 
			if_else(!is.na(all_participant_visits), 
				as.character(all_participant_visits), ''),
		viremic_participant_visits = if_else(round == 15, '..', as.character(viremic_participant_visits))) %>%
	select(-c(n_hiv_participant_visits, n_viremic_participant_visits, n_pretreat_participant_visits))



write_tsv(t1, "tables/table_1.tsv")



