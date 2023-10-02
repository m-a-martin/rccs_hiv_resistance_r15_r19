library(tidyverse)
source('scripts/utils.R')

sort_cols = function(x){
		class_order = c("independence" = 1, "exchangeable" = 2, "ar1" = 3, "unstructured" = 4)
		type_order = c("Prev. % (95% CI)" = 1, "Prev. ratio (95% CI)" = 2, "p-value" = 3, "spacer" = 4)
		non_round_x = x[x!= "round"]
		s1 = class_order[str_split(non_round_x, "_", simplify=T)[,1]]
		s2 = type_order[str_split(non_round_x, "_", simplify=T)[,2]]
		scores = tibble(val = non_round_x, s1 = s1, s2 = s2) %>%
			arrange(s1, s2)
		return(c("round", scores$val))
}

pred_file = 'models/pretreat_int97a_amongPLWHIV_gee_prev_pred.tsv'
rr_file = 'models/pretreat_int97a_amongPLWHIV_gee_rr.tsv'

pred = read_tsv(pred_file)
# format for table
pred = pred %>%
	select(-se.fit) %>%
	mutate(across(fit:upr, ~format_digit(.x*100))) %>%
	unite("val", fit:lwr, sep=' (') %>%
	unite("val", val:upr, sep=', ') %>%
	mutate(
		val = paste(val, ')', sep=''), 
		type = paste(c, "Prev. % (95% CI)", sep='_')) %>%
		select(-c)

rr = read_tsv(rr_file) %>%
	mutate(across(RR:UCI, ~format_digit(.x))) %>%
	unite("val", RR:LCI, sep=' (') %>%
	unite("val", val: UCI, sep=', ') %>%
	mutate(
		val = paste(val, ')', sep=''),
		val = if_else(`var` == "(Intercept)", 'ref', val),
		P = if_else(P < 1E-4, "<0.0001", as.character(signif(P, 2))),
		P = if_else(`var` == "(Intercept)", 'ref', P),
		`var` = as.numeric(str_split(`var`, 'round', simplify=T)[,2]),
		`var` = if_else(is.na(`var`), 
			min(`var`, na.rm=TRUE)-1,
			`var`),
		type = paste(c, "Prev. ratio (95% CI)", sep='_')) %>%
	rename(c("round"='var')) %>%
	select(-c)

	rr = 
		bind_rows(
			rr %>% select(-P),
			rr %>% select(-val) %>% rename(c("val" = "P")) %>%
				mutate(type = paste(str_split(type, "_", simplify=T)[,1], "p-value", sep='_')))

	t = bind_rows(pred, rr)
	t = t %>% pivot_wider(names_from=type, values_from=val)
	# add blank columns
	for (c in unique(read_tsv(pred_file)$c)){
		t[paste(c, 'spacer', sep='_')] = NA
	}
	# sort
	t = t %>% select(sort_cols(colnames(t))) %>%
		arrange(round)

write_tsv(t, 'tables/table_s31.tsv')