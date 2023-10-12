library(tidyverse)
source('scripts/utils.R')



sort_cols = function(x, excl=c()){
	round_order = c("nnrti" = 1, "nrti" = 2, "pi" = 3, "all" = 4)
	type_order = c("Coefficient (95% CI)" = 1, "p-value" = 3, "spacer"=4)
	non_excl_x = x[!(x %in% excl)]
	s1 = round_order[(str_split(non_excl_x, "_", simplify=T)[,1])]
	s2 = type_order[str_split(non_excl_x, "_", simplify=T)[,2]]
	scores = tibble(val = non_excl_x, s1 = s1, s2 = s2) %>%
		arrange(s1, s2)
	return(c(excl, scores$val))
}


# todo infer this somehow?
ref_cat = c('Sex'='Female', 'Age category'='(14,24] years old', 
	'Survey round' = 'round 15')

ref_table = tibble(var_type = names(ref_cat), Variable = ref_cat, `Coefficient (95% CI)` = 'ref', 
	`p-value` = 'ref')

row_headers = tibble(var_type = c(names(ref_cat), 'Viral load'), Variable = c(names(ref_cat), 'Viral load'), 
	`Coefficient (95% CI)` = '', 
	`p-value` = '' )

ref_table = bind_rows(row_headers, ref_table)


all_col=list()
for (i_class in c('nnrti')){
	o_class = read_tsv(paste(c('models/pretreat_inc_', i_class, '_weights_raw.tsv'), collapse=''), show_col_types = FALSE) %>%
		mutate(RR = round(RR, 2),
			LCI = round(LCI, 2),
			UCI = round(UCI, 2)) %>%
	rename(var = `...1`) %>%
	unite("val", RR:LCI, sep=' (') %>%
	unite("val", val: UCI, sep=', ') %>%
	mutate(
		val = paste(val, ')', sep=''),
		#val = if_else(`var` == "(Intercept)", 'ref', val),
		P = if_else(P < 1E-4, "<0.0001", as.character(signif(P, 2))),
		var_type = case_when(
			str_detect(var, 'comm_type') ~ 'Community type',
			str_detect(var, 'sex') ~ 'Sex',
			str_detect(var, 'age_cat') ~ 'Age category',
			str_detect(var, 'vl') ~ 'Viral load',
			str_detect(var, 'round') ~ 'Survey round'),
		var = if_else(
			str_detect(var, ":"),
			paste(
				str_split(var, ":", simplify=T)[,2], 
				str_split(var, ":", simplify=T)[,1],
				sep=':'),
			var),
		var = gsub(':',' \u00D7 ',var),
		var = gsub('comm_type', '', var),
		var = gsub('Fishing', 'Fishing community', var),
		var = gsub('Trading', 'Trading community', var),
		var = gsub('age_cat', '', var),
		var = gsub(']', '] years old', var),
		var = gsub('round', 'round ', var),
		var = gsub('sexM', 'Male', var),
		var = gsub('log10vl', 'log10 VL', var)) %>%
	rename(Variable = var, `Coefficient (95% CI)` = val, `p-value` = P)
	o_class = bind_rows(o_class, ref_table) %>%
		mutate(class = i_class)
	all_col[[i_class]] = o_class
}

all_col = bind_rows(all_col)


all_col = all_col %>% 
	pivot_longer(-c(Variable, class), names_to='name', values_to='var') %>%
	pivot_wider(names_from=c('class', 'name'), values_from='var') %>%
	mutate('pi_spacer' = '', 'nnrti_spacer' = '', 'nrti_spacer' = '')


all_col = all_col %>% select(sort_cols(colnames(all_col), c('Variable')))

all_col = all_col %>% arrange(!is.na(nnrti_var_type), 
	!(nnrti_var_type == 'Survey round'), 
	nnrti_var_type, 
	!(`nnrti_p-value` == ''),
	!(`nnrti_p-value` == 'ref'))

# finally, drop unnecessary rows
all_col = all_col %>% select(names(all_col)[!str_detect(names(all_col), 'var_type')])

write_tsv(all_col, 'tables/table_s8.tsv')
