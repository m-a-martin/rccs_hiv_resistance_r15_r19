library(tidyverse)


d = read_tsv('models/all_multi_amongPLWHIV_prev_pred.tsv')

t = d %>% select(-se.fit) %>%
	mutate(across(fit:upr, ~format_digit(.x*100))) %>%
	unite("val", fit:lwr, sep=' (') %>%
	unite("val", val:upr, sep=', ') %>%
	mutate(
		val = paste(val, ')', sep='')) %>%
	rename(`Prev. % (95% CI)` = val)

write_tsv(t, 'tables/table_s27.tsv')