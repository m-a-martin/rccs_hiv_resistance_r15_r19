library(tidyverse)
source('scripts/utils.R')


sort_cols = function(x){
	#round_order = c("nnrti" = 1, "nrti" = 2, "pi" = 3)
	type_order = c("Prev. (95% CI)" = 1, "Prev. ratio (95% CI)" = 2, "p-value" = 3, "spacer" = 4)
	non_round_x = x[x!= "s_var" & x != 's']
	s1 = as.numeric(str_split(non_round_x, "_", simplify=T)[,1])
	s2 = type_order[str_split(non_round_x, "_", simplify=T)[,2]]
	scores = tibble(val = non_round_x, s1 = s1, s2 = s2) %>%
		arrange(s1, s2)
	return(c("s", "s_var", scores$val))
}



#### PRE-TREATMENT rtK103N BY ROUND AND SEX AMONG PLWHIV ####

pretreat_dr_pred_file = 
'models/pretreat_dr_muts_amongPLWHIV_prev_pred_stratified.tsv'
pretreat_dr_rr_file = 'models/pretreat_dr_muts_amongPLWHIV_rr_pred_stratified.tsv'

# format for table
pretreat_dr_pred = read_tsv(pretreat_dr_pred_file, show_col_types=FALSE) %>%
	filter(mut == 'rtK103N') %>%
	select(-se.fit, -mut) %>%
	mutate(fit = round(fit*100, 2), 
		lwr = round(lwr*100, 2), upr = round(upr*100, 2)) %>%
	unite("val", fit:lwr, sep=' (') %>%
	unite("val", val:upr, sep=', ') %>%
	mutate(
		val = paste(val, ')', sep=''),
		type = paste(round, "Prev. (95% CI)", sep='_')) %>%
	rename(s_var = sex)


# todo infer this somehow?
ref_cat = c('sex'='F', 'age_cat'='(14,24]', 'comm_type' = 'Agrarian')

pretreat_dr_rr = read_tsv(pretreat_dr_rr_file) %>%
	filter(mut == 'rtK103N') %>%	
	select(-mut) %>%
	mutate(RR = round(RR, 2),
		LCI = round(LCI, 2),
		UCI = round(UCI, 2)) %>%
	unite("val", RR:LCI, sep=' (') %>%
	unite("val", val: UCI, sep=', ') %>%
	mutate(
		val = paste(val, ')', sep=''),
		val = if_else(`var` == "(Intercept)", 'ref', val),
		P = if_else(P < 1E-4, "<0.0001", as.character(signif(P,2))),
		P = if_else(`var` == "(Intercept)", 'ref', P),
		round = as.numeric(str_split(str_split(
			var, ':', simplify=T)[,1], 'round', simplify=T)[,2]),
		round = replace_na(round, min(round, na.rm=TRUE)-1),
		s_var = str_split(str_split(var, ":", simplify=T)[,2], 's_var', simplify=T)[,1],
		s_var = case_when(
			s_var == "" & !grepl("sex", var) ~ ref_cat[s],
			s_var == "" & grepl("sex", var) ~ str_split(var, s, simplify=T)[,2],
			s_var != '' ~ str_split(s_var, s, simplify=T)[,2]),
		type = paste(round, "Prev. ratio (95% CI)", sep='_')) %>%
	select(-var)


pretreat_dr_rr = 
	bind_rows(
		pretreat_dr_rr %>% select(-P),
	pretreat_dr_rr %>% select(-val) %>% rename(c("val" = "P")) %>%
		mutate(type = paste(str_split(type, "_", simplify=T)[,1], "p-value", sep='_')))



ts = bind_rows(pretreat_dr_pred, pretreat_dr_rr)


ts = ts %>% select(-round) %>% pivot_wider(names_from="type", values_from="val") %>%
	mutate(`15_spacer` = '', `16_spacer` = '', `17_spacer` = '', `18_spacer` = '')




ts = ts %>%
	select(sort_cols(colnames(ts)))

ts = bind_rows(ts, tibble(s = unique(ts$s))) %>% 
	arrange(s, !is.na(s_var), s_var)

write_tsv(ts, 'tables/table_s30.tsv', na='')


