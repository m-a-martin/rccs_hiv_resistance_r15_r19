library(tidyverse)
source('scripts/utils.R')


sort_cols = function(x){
	class_order = c("nnrti" = 1, "nrti" = 2, "pi" = 3)
	type_order = c("Prev. % (95% CI)" = 1, "Prev. ratio (95% CI)" = 2, "p-value" = 3, "spacer" = 10)
	non_round_x = x[x!= "round"]
	s1 = class_order[str_split(non_round_x, "_", simplify=T)[,1]]
	s2 = type_order[str_split(non_round_x, "_", simplify=T)[,2]]
	scores = tibble(val = non_round_x, s1 = s1, s2 = s2) %>%
		arrange(s1, s2)
	return(c("round", scores$val))
}


#### PRE-TREATMENT RESISTANCE BY ROUND AMONG PLWHIV ####

dat_file = 'data/rakai_drug_resistance_categorized.tsv'
table_vars = c("age_cat", "comm_type", "sex")



pretreat_dr_pred_file = 'models/pretreat_amongPar_prev_pred.tsv'
pretreat_dr_rr_file = 'models/pretreat_amongPar_rr.tsv'

pretreat_dr_pred = read_tsv(pretreat_dr_pred_file)
# format for table
pretreat_dr_pred = pretreat_dr_pred %>%
	select(-se.fit) %>%
	mutate(across(fit:upr, ~format_digit(.x*100))) %>%
	unite("val", fit:lwr, sep=' (') %>%
	unite("val", val:upr, sep=', ') %>%
	mutate(
		val = paste(val, ')', sep=''),
		type = paste(class, "Prev. % (95% CI)", sep='_')) %>%
	select(-class)

pretreat_dr_rr = read_tsv(pretreat_dr_rr_file) %>%
	mutate(across(RR:UCI, ~format_digit(.x))) %>%
	unite("val", RR:LCI, sep=' (') %>%
	unite("val", val: UCI, sep=', ') %>%
	mutate(
		val = paste(val, ')', sep=''),
		val = if_else(`var` == "(Intercept)", 'ref', val),
		P = if_else(P < 1E-4, "<0.0001", as.character(signif(P,2))),
		P = if_else(`var` == "(Intercept)", 'ref', P),
		`var` = as.numeric(str_split(`var`, 'round', simplify=T)[,2]),
		`var` = if_else(is.na(`var`), 
			min(`var`, na.rm=TRUE)-1,
			`var`),
		type = paste(class, "Prev. ratio (95% CI)", sep='_')) %>%
	rename(c("round"='var')) %>%
	select(-class)

pretreat_dr_rr = 
	bind_rows(
		pretreat_dr_rr %>% select(-P),
	pretreat_dr_rr %>% select(-val) %>% rename(c("val" = "P")) %>%
		mutate(type = paste(str_split(type, "_", simplify=T)[,1], "p-value", sep='_')))


ts = bind_rows(pretreat_dr_pred, pretreat_dr_rr)

ts = ts %>% pivot_wider(names_from=type, values_from=val) %>% mutate(
	'nnrti_spacer' = '', 'nrti_spacer' = '', 'pi_spacer' = '')
# sort
ts = ts %>% select(sort_cols(colnames(ts))) %>%
	arrange(round)

write_tsv(ts, 'tables/table_s21.tsv')
