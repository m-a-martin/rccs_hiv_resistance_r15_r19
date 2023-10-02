library(tidyverse)
source('scripts/utils.R')



sort_cols = function(x, excl=c()){
	round_order = c("gen" = 1, "nnrti" = 2, "nrti" = 3, "pi" = 4)
	type_order = c("Coefficient (95% CI)" = 1, "p-value" = 3, "spacer"=4)
	non_excl_x = x[!(x %in% excl)]
	s1 = round_order[(str_split(non_excl_x, "_", simplify=T)[,1])]
	s2 = type_order[str_split(non_excl_x, "_", simplify=T)[,2]]
	scores = tibble(val = non_excl_x, s1 = s1, s2 = s2) %>%
		arrange(s1, s2)
	return(c(excl, scores$val))
}

all_col=list()
for (i_class in c('nnrti')){
	o_class = read_tsv(paste(c('models/pretreat_inc_', i_class, '_weights_raw.tsv'), collapse='')) %>%
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
		var = gsub(':',' & ',var),
		var = gsub('comm_type', '', var),
		var = gsub('Fishing', 'Fishing community', var),
		var = gsub('Trading', 'Trading community', var),
		var = gsub('age_cat', '', var),
		var = gsub(']', '] years old', var),
		var = gsub('round', 'round ', var),
		var = gsub('sexM', 'Male', var),
		var = gsub('missing_vlTRUE', 'missing VL', var),
		var = gsub('log10vl', 'log10 VL', var),
		class = i_class) %>%
	rename(Variable = var, `Coefficient (95% CI)` = val, `p-value` = P)
	all_col[[i_class]] = o_class
}

all_col = bind_rows(all_col)


all_col = all_col %>% 
	pivot_longer(-c(Variable, class), names_to='name', values_to='var') %>%
	pivot_wider(names_from=c('class', 'name'), values_from='var')
all_col = all_col %>% select(sort_cols(colnames(all_col), c('Variable')))

write_tsv(all_col, 'tables/table_s8.tsv')
