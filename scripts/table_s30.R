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
treat_dr_muts_prev_file = 'models/treat_dr_muts_amongPLWHIV_prev_pred.tsv'
treat_dr_muts_prev = read_tsv(treat_dr_muts_prev_file)


plot_muts = (treat_dr_muts_prev %>% filter(round == 18) %>% arrange(-fit) %>%
	slice(1:15) %>% select(mut) %>% mutate(idx = seq(1,n())))

treat_dr_muts_prev = treat_dr_muts_prev %>% filter(mut %in% plot_muts$mut) %>%
	select(-se.fit) %>%
	mutate_at(vars(fit, lwr, upr), ~round(.x*100, 2)) %>%
	unite("val", fit:lwr, sep=' (') %>%
	unite("val", val:upr, sep=', ') %>%
	mutate(
		val = paste(val, ')', sep=''),
		round = paste(round, "_Prev. (95% CI)", sep='')) %>%
	pivot_wider(names_from=round, values_from=val)


t = treat_dr_muts_prev %>% 
	# sort rows
	left_join(plot_muts, by="mut") %>%
	arrange(idx) %>%
	select(-idx)

for (r in unique(read_tsv(treat_dr_muts_prev_file)$round)){
	t[paste(r, "_spacer", sep='')] = NA
}

t = t %>% select(sort_cols(colnames(t)))
write_tsv(t, 'tables/table_s30.tsv', na='')

