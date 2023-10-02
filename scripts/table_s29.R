library(tidyverse)
source('scripts/utils.R')


sort_cols = function(x){
		type_order = c("Prev. (95% CI)" = 1, "Prev. ratio (95% CI)" = 2, "p-value" = 3, "spacer" = 4)
		non_mut_x = x[x!= "mut"]
		s1 = as.numeric(str_split(non_mut_x, "_", simplify=T)[,1])
		s2 = type_order[str_split(non_mut_x, "_", simplify=T)[,2]]
		scores = tibble(val = non_mut_x, s1 = s1, s2 = s2) %>%
			arrange(s1, s2)
		return(c("mut", scores$val))
}


# risk ratio of 15 most frequent mutations in pretreated data compared to R15
pretreat_dr_muts_prev_file = 'models/pretreat_dr_muts_amongPLWHIV_prev_pred.tsv'
pretreat_dr_muts_prev = read_tsv(pretreat_dr_muts_prev_file)
plot_muts = (pretreat_dr_muts_prev %>% filter(round == 18) %>% arrange(-fit) %>%
	slice(1:15) %>% select(mut) %>% mutate(idx = seq(1,n())))

pretreat_dr_muts_prev = pretreat_dr_muts_prev %>% filter(mut %in% plot_muts$mut) %>%
	select(-se.fit) %>%
	mutate_at(vars(fit, lwr, upr), ~round(.x*100, 2)) %>%
	unite("val", fit:lwr, sep=' (') %>%
	unite("val", val:upr, sep=', ') %>%
	mutate(
		val = paste(val, ')', sep=''),
		round = paste(round, "_Prev. (95% CI)", sep='')) %>%
	pivot_wider(names_from=round, values_from=val)


pretreat_dr_muts_rr_file = 'models/pretreat_dr_muts_amongPLWHIV_rr_pred.tsv'
pretreat_dr_muts_rr = read_tsv(pretreat_dr_muts_rr_file)
pretreat_dr_muts_rr = pretreat_dr_muts_rr %>% filter(mut %in% plot_muts$mut) %>%
	mutate_at(vars(RR, LCI, UCI), ~round(., 2)) %>%
	unite("val", RR:LCI, sep=' (') %>%
	unite("val", val:UCI, sep=', ') %>%
	mutate(
		val = paste(val, ')', sep=''),
		val = if_else(`var` == "(Intercept)", 'ref', val),
		var = as.numeric(str_split(`var`, 'round', simplify=T)[,2]),
		var = if_else(is.na(`var`), 
			min(`var`, na.rm=TRUE)-1,
			`var`), 
		P = if_else(P < 1e-3, "<0.001", as.character(signif(P, 2)))) %>%
	rename(`Prev. ratio (95% CI)` = val, `p-value` = P) %>%
	pivot_longer(-c(mut, var)) %>%
	mutate(name = paste(var, name, sep="_")) %>%
	select(-var) %>%
	pivot_wider(names_from=name, values_from=value)


t5 = pretreat_dr_muts_prev %>% left_join(pretreat_dr_muts_rr, by="mut") %>%
	# sort rows
	left_join(plot_muts, by="mut") %>%
	arrange(idx) %>%
	select(-idx)

for (r in unique(read_tsv(pretreat_dr_muts_prev_file)$round)){
	t5[paste(r, "_spacer", sep='')] = NA
}

t5 = t5 %>% select(sort_cols(colnames(t5)))
write_tsv(t5, 'tables/table_s29.tsv', na='')

