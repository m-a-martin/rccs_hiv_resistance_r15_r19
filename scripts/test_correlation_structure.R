library(tidyverse)
library(lmtest)
library(sandwich)
source('scripts/utils.R')
library(geepack)

dat_file = 'data/rakai_drug_resistance_categorized.tsv'
pretreat_w_file = 'models/pretreat_weights.tsv'
treat_w_file = 'models/treat_weights.tsv'
all_w_file = 'models/all_weights.tsv'
# read in data
hiv_dr_cat = read_tsv(dat_file, show_col_types = FALSE)
# sampling weights
pretreat_weights = read_tsv(pretreat_w_file, show_col_types = FALSE)
treat_weights = read_tsv(treat_w_file, show_col_types = FALSE)
all_weights = read_tsv(all_w_file, show_col_types = FALSE)
# merge em
pretreat_dat = hiv_dr_cat %>% left_join(pretreat_weights, by=c("study_id", "round")) %>% 
	filter(pre_treatment & finalhiv == "P" & viremic) %>% arrange(study_id, round)

# iterate over drugs and correlation structures to assess best fit
all_pretreat_q = list()
all_pretreat_m = list()
for (class in c("nnrti", "nrti", "pi")){
	class_dat = pretreat_dat %>% 
		# TODO, change paste to handle variable better
		rename(
			r = !!class, 
			w = !!paste("w_", class, sep='')) %>%
		filter(!is.na(r)) %>%
		mutate(r = 1 * (r == 'intermediate/high'), round = as.character(round))
	class_dat = class_dat %>% 
		left_join(
			class_dat %>% 
				select(study_id) %>% 
				unique() %>%
				mutate(idx = seq(1,n())), by="study_id") %>% arrange(idx, round)
	#class_dat = class_dat %>% group_by(study_id) %>% filter(n() > 1)
	for (corr in c("independence", "exchangeable", "ar1", "unstructured")){
		m = geeglm(r ~ round, data=class_dat, id = idx, family=poisson(link="log"), corstr=corr, weights=w)
		q = QIC(m)
		all_pretreat_q[[paste(class, corr, sep="_")]] = tibble(type = names(q), val=q, class = class, corr=corr)
		all_pretreat_m[[paste(class, corr, sep='-')]] = tibble(round = names(m$coefficients), 
			coefficient=m$coefficients, std.err = summary(m)$coefficients[,2], class = class, corr=corr)
	}
}

pretreat_q_table = bind_rows(all_pretreat_q) %>% pivot_wider(names_from=type, values_from=val)  %>% group_by(class) %>% arrange(QIC)
pretreat_m_table = bind_rows(all_pretreat_m)

write_tsv(pretreat_q_table, 'models/pretreat_q_table.tsv')
write_tsv(pretreat_m_table, 'models/pretreat_m_table.tsv')

treat_dat = hiv_dr_cat %>% left_join(treat_weights, by=c("study_id", "round")) %>% 
	filter(!pre_treatment & finalhiv == "P" & viremic & round > 16 & round < 19)  %>% arrange(study_id, round)

# iterate over drugs and correlation structures to assess best fit
all_treat_q = list()
all_treat_m = list()
for (class in c("nnrti", "nrti", "pi")){
	class_dat = treat_dat %>% 
		rename(
			r = !!class, 
			w = !!paste("w_", class, sep='')) %>%
		filter(!is.na(r)) %>%
		mutate(r = 1 * (r == 'intermediate/high'), round = as.character(round))
	class_dat = class_dat %>% 
		left_join(
			class_dat %>% 
				select(study_id) %>% 
				unique() %>%
				mutate(idx = seq(1,n())), by="study_id") %>% arrange(idx, round)
	#class_dat = class_dat %>% group_by(study_id) %>% filter(n() > 1)
	for (corr in c("independence", "exchangeable")){
		m = geeglm(r ~ round, data=class_dat, id = idx, family=poisson(link="log"), corstr=corr, weights=w)
		q = QIC(m)
		all_treat_q[[paste(class, corr, sep="_")]] = tibble(type = names(q), val=q, class = class, corr=corr)
		all_treat_m[[paste(class, corr, sep='-')]] = tibble(round = names(m$coefficients), 
			coefficient=m$coefficients, std.err = summary(m)$coefficients[,2], class = class, corr=corr)
	}
	if (class_dat %>% group_by(study_id) %>% summarise(n=n()) %>% select(n) %>% max() > 2){
		for (corr in c("ar1", "unstructured")){
			m = geeglm(r ~ round, data=class_dat, id = idx, family=poisson(link="log"), corstr=corr, weights=w)
			q = QIC(m)
				all_treat_q[[paste(class, corr, sep="_")]] = tibble(type = names(q), val=q, class = class, corr=corr)
				all_treat_m[[paste(class, corr, sep='-')]] = tibble(round = names(m$coefficients), 
					coefficient=m$coefficients, std.err = summary(m)$coefficients[,2], class = class, corr=corr)
		}
	}
}


treat_q_table = bind_rows(all_treat_q) %>% pivot_wider(names_from=type, values_from=val)  %>% group_by(class) %>% arrange(QIC)
treat_m_table = bind_rows(all_treat_m)

write_tsv(treat_q_table, 'models/treat_q_table.tsv')
write_tsv(treat_m_table, 'models/treat_m_table.tsv')


all_dat = hiv_dr_cat %>% left_join(all_weights, by=c("study_id", "round")) %>% 
	filter(finalhiv == "P" & viremic & round > 16 & round < 19)  %>% arrange(study_id, round)

# iterate over drugs and correlation structures to assess best fit
all_all_q = list()
all_all_m = list()
for (class in c("nnrti", "nrti", "pi")){
	class_dat = all_dat %>% 
		rename(
			r = !!class, 
			w = !!paste("w_", class, sep='')) %>%
		filter(!is.na(r)) %>%
		mutate(r = 1 * (r == 'intermediate/high'), round = as.character(round))
	class_dat = class_dat %>% 
		left_join(
			class_dat %>% 
				select(study_id) %>% 
				unique() %>%
				mutate(idx = seq(1,n())), by="study_id") %>% arrange(idx, round)
	#class_dat = class_dat %>% group_by(study_id) %>% filter(n() > 1)
	for (corr in c("independence", "exchangeable")){
		m = geeglm(r ~ round, data=class_dat, id = idx, family=poisson(link="log"), corstr=corr, weights=w)
		q = QIC(m)
		all_all_q[[paste(class, corr, sep="_")]] = tibble(type = names(q), val=q, class = class, corr=corr)
		all_all_m[[paste(class, corr, sep='-')]] = tibble(round = names(m$coefficients), 
			coefficient=m$coefficients, std.err = summary(m)$coefficients[,2], class = class, corr=corr)
	}
	if (class_dat %>% group_by(study_id) %>% summarise(n=n()) %>% select(n) %>% max() > 2){
		for (corr in c("ar1", "unstructured")){
			m = geeglm(r ~ round, data=class_dat, id = idx, family=poisson(link="log"), corstr=corr, weights=w)
			all_all_q[[paste(class, corr, sep="_")]] = tibble(type = names(q), val=q, class = class, corr=corr)
			all_all_m[[paste(class, corr, sep='-')]] = tibble(round = names(m$coefficients), 
				coefficient=m$coefficients, std.err = summary(m)$coefficients[,2], class = class, corr=corr)
		}
	}
}


all_q_table = bind_rows(all_all_q) %>% pivot_wider(names_from=type, values_from=val)  %>% group_by(class) %>% arrange(QIC)
all_m_table = bind_rows(all_all_m)

write_tsv(all_q_table, 'models/all_q_table.tsv')
write_tsv(all_m_table, 'models/all_m_table.tsv')



