library(tidyverse)

source('scripts/utils.R')



format_col = function(x1, x2){
		return(paste(paste(paste(x1, " (", sep=''),
			round((x1 / x2)*100, 2), sep=''),
			'%)', sep=''))
}


chisq_p = function(x, a, b){
	return(chisq.test(rbind(na.omit(x[[a]]), na.omit(x[[b]])))$p.value)
}





dat_file = 'data/rakai_drug_resistance_categorized.tsv'
table_vars = c("age_cat", "comm_type", "sex", "vl_cat")

# read in data and filter out bad rounds and age categories
hiv_dr_cat = read_tsv(dat_file) %>% 
	filter(
		age_cat != '(50, 100]' & 
		round > 14)


# viremic hiv+ participants
v_par = list()
v_par[['round']] = hiv_dr_cat %>%
		filter(finalhiv == 'P' & viremic) %>%
		rename(v = round) %>%
		mutate(v = as.character(v)) %>%
		group_by(v) %>%
		summarise(
			n_viremic_participants = n(), 
			n_viremic_sequenced = sum(!is.na(pi)), 
			.group='drop') %>%
		mutate(var = 'Survey round') %>%
		select(var, v, n_viremic_participants, n_viremic_sequenced)

for (i_var in table_vars){
	v_par[[i_var]] = hiv_dr_cat %>%
		filter(finalhiv == 'P' & viremic) %>%
		rename(v = !!i_var) %>%
		group_by(v) %>%
		summarise(
			n_viremic_participants=n(), 
			n_viremic_sequenced = sum(!is.na(pi)), 
			.groups='drop') %>%
		mutate('var' = str_to_sentence(var_rename[i_var])) %>%
		select(var, v, n_viremic_participants, n_viremic_sequenced)
}
v_par = bind_rows(v_par)

# pre-treatment hiv+ participants
p_par = list()
p_par[['round']] = hiv_dr_cat %>%
		filter(finalhiv == 'P' & viremic & pre_treatment) %>%
		rename(v = round) %>%
		mutate(v = as.character(v)) %>%
		group_by(v) %>%
		summarise(
			n_pretreat_participants = n(), 
			n_pretreat_sequenced = sum(!is.na(pi)), 
			.group='drop') %>%
		mutate(var = 'Survey round') %>%
		select(var, v, n_pretreat_participants, n_pretreat_sequenced)
for (i_var in table_vars){
	p_par[[i_var]] = hiv_dr_cat %>%
		filter(finalhiv == 'P' & viremic & pre_treatment) %>%
		rename(v = !!i_var) %>%
		group_by(v) %>%
		summarise(n_pretreat_participants=n(), 
			n_pretreat_sequenced = sum(!is.na(pi)), 
			.groups='drop') %>%
		mutate('var' = str_to_sentence(var_rename[i_var])) %>%
		select(var, v, n_pretreat_participants, n_pretreat_sequenced)
}
p_par = bind_rows(p_par)


# treatment hiv+ participants
t_par = list()
t_par[['round']] = hiv_dr_cat %>%
		filter(finalhiv == 'P' & viremic & !pre_treatment) %>%
		rename(v = round) %>%
		mutate(v = as.character(v)) %>%
		group_by(v) %>%
		summarise(
			n_treat_participants = n(), 
			n_treat_sequenced = sum(!is.na(pi)), 
			.group='drop') %>%
		mutate(var = 'Survey round') %>%
		select(var, v, n_treat_participants, n_treat_sequenced)


for (i_var in table_vars){
	t_par[[i_var]] = hiv_dr_cat %>%
		filter(finalhiv == 'P' & viremic & !pre_treatment) %>%
		rename(v = !!i_var) %>%
		group_by(v) %>%
		summarise(n_treat_participants=n(), 
			n_treat_sequenced = sum(!is.na(pi)), 
			.groups='drop') %>%
		mutate('var' = str_to_sentence(var_rename[i_var])) %>%
		select(var, v, n_treat_participants, n_treat_sequenced)
}
t_par = bind_rows(t_par)


# now column bind all of these
table = v_par %>% 
	inner_join(p_par,by=c("var", "v")) %>%
	left_join(t_par, by=c("var", "v")) %>%
	arrange(var, v) 

# get overall data for each round
oa = hiv_dr_cat %>% 
	summarize(
		n_viremic_participants = sum(viremic & (!is.na(finalhiv) & finalhiv == 'P')),
		n_viremic_sequenced = sum(viremic & (!is.na(finalhiv) & finalhiv == 'P') & (!is.na(pi))),
		n_pretreat_participants = sum(pre_treatment & viremic & (!is.na(finalhiv) & finalhiv == 'P')),
		n_pretreat_sequenced = sum(pre_treatment & viremic & (!is.na(finalhiv) & finalhiv == 'P') & (!is.na(pi))),
		n_treat_participants = sum(!pre_treatment & viremic & (!is.na(finalhiv) & finalhiv == 'P')),
		n_treat_sequenced = sum(!pre_treatment & viremic & (!is.na(finalhiv) & finalhiv == 'P')  & 
			!is.na(pi))) %>%
	mutate(var = 'Overall', v = NA, n_viremic_p = NA, n_pretreat_p = NA, n_treat_p = NA)



# need to add chi2 p_values comparing all categories to the number of participants
table  = bind_rows(
	table, 
	bind_rows(table %>% group_by(var) %>% group_map(~
		tibble(
			var = .y[[1]],
			v = NA,
			n_viremic_participants = NA,
			n_viremic_p = chisq.test(rbind(.$n_viremic_participants, .$n_viremic_sequenced))$p.value,
			n_preatreat_participants = NA,
			n_pretreat_p = chisq.test(rbind(.$n_pretreat_participants, .$n_pretreat_sequenced))$p.value,
			n_treat_participants = NA,
			n_treat_p = chisq.test(rbind(na.omit(.$n_treat_participants), na.omit(.$n_treat_sequenced)))$p.value))),
	oa) %>%
	arrange(!(var == "Overall"), !is.na(var), var, !is.na(v), v) %>%
	select(var, v, n_viremic_participants, n_viremic_sequenced, n_viremic_p, 
		n_pretreat_participants, n_pretreat_sequenced, n_pretreat_p, n_treat_participants, n_treat_sequenced, n_treat_p)

# format columns
table = table %>% 
	mutate(
		n_viremic_sequenced = 
			if_else(!is.na(n_viremic_sequenced), format_col(n_viremic_sequenced, n_viremic_participants), NA),
		n_pretreat_sequenced = 
			if_else(!is.na(n_pretreat_sequenced), format_col(n_pretreat_sequenced, n_pretreat_participants), NA),
		n_treat_sequenced = 
			if_else(!is.na(n_treat_sequenced), format_col(n_treat_sequenced, n_treat_participants), NA)) %>%
	mutate_at(vars(n_viremic_p, n_pretreat_p, n_treat_p), ~if_else(!is.na(.), if_else(. < 1E-4, "<1e-4", as.character(signif(., 2))), NA)) %>%
	mutate_all(~if_else(is.na(.), "", as.character(.))) %>%
	mutate(viremic_spacer = '', pretreat_spacer = '') %>%
	select(var, v, n_viremic_participants, n_viremic_sequenced, n_viremic_p, viremic_spacer, 
		n_pretreat_participants, n_pretreat_sequenced, n_pretreat_p, pretreat_spacer, 
		n_treat_participants, n_treat_sequenced, n_treat_p)





write_tsv(table, 'tables/table_s17.tsv')
