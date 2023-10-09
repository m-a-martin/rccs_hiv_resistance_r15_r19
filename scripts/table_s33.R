library(tidyverse)
source('scripts/utils.R')


sort_cols = function(x, excl=c()){
	#round_order = c("nnrti" = 1, "nrti" = 2, "pi" = 3)
	type_order = c("Prev. (95% CI)" = 1, "Prev. ratio (95% CI)" = 2, "p-value" = 3, "spacer"=4)
	non_excl_x = x[!(x %in% excl)]
	s1 = as.numeric(str_split(non_excl_x, "_", simplify=T)[,1])
	s2 = type_order[str_split(non_excl_x, "_", simplify=T)[,2]]
	scores = tibble(val = non_excl_x, s1 = s1, s2 = s2) %>%
		arrange(s1, s2)
	return(c(excl, scores$val))
}



#### PRE-TREATMENT RESISTANCE BY ROUND AMONG PLWHIV ####


pretreat_dr_muts_prev_file = 'models/pretreat_int97a_amongPLWHIV_prev_pred.tsv'
pretreat_dr_muts_rr_file =  'models/pretreat_int97a_amongPLWHIV_rr.tsv'

# format for table
pretreat_dr_pred = read_tsv(pretreat_dr_muts_prev_file) %>%
	select(-se.fit) %>%
	mutate(fit = round(fit*100, 2), 
		lwr = round(lwr*100, 2), upr = round(upr*100, 2)) %>%
	unite("val", fit:lwr, sep=' (') %>%
	unite("val", val:upr, sep=', ') %>%
	mutate(
		val = paste(val, ')', sep=''),
		type = paste(round, "Prev. (95% CI)", sep='_')) 


# todo infer this somehow?
ref_cat = c('sex'='F', 'age_cat'='(14,24]', 'comm_type' = 'Agrarian', 'subtype' = "A1")

pretreat_dr_rr = read_tsv(pretreat_dr_muts_rr_file) %>%
	mutate(RR = round(RR, 3),
		LCI = round(LCI, 3),
		UCI = round(UCI, 3)) %>%
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
		s_var = str_split(str_split(var, ":", simplify=T)[,2], 's_var', simplify=T)[,2],
		s_var = case_when(
			s_var == "" & !grepl("s_var", var) ~ ref_cat[s],
			s_var == "" & grepl("s_var", var) ~ str_split(var, "s_var", simplify=T)[,2],
			s_var != '' ~ s_var),
		type = paste(round, "Prev. ratio (95% CI)", sep='_')) %>%
	select(-var)

pretreat_dr_rr = 
	bind_rows(
		pretreat_dr_rr %>% select(-P),
	pretreat_dr_rr %>% select(-val) %>% rename(c("val" = "P")) %>%
		mutate(type = paste(str_split(type, "_", simplify=T)[,1], "p-value", sep='_')))



ts11 = bind_rows(pretreat_dr_pred, pretreat_dr_rr)
ts11 = ts11 %>% select(-round) %>% pivot_wider(names_from="type", values_from="val")

# blank cols
for (r in unique(read_tsv(pretreat_dr_muts_prev_file)$round)){
	ts11[paste(r, "_spacer", sep='')] = NA
}

ts11 = ts11 %>%
	select(sort_cols(colnames(ts11), c('s', 's_var')))


#blank rows
ts11 = bind_rows(
	ts11, 
	 tibble(s = unique(ts11$s))) %>%
	arrange(s, !is.na(s_var), s_var)


write_tsv(ts11, 'tables/table_s33.tsv', na='')
