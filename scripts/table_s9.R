library(tidyverse)


order=c('independence'=1, 'exchangeable'=2, 'ar1'=3, 'unstructured'=4)
pretreat_q_dat = read_tsv('models/pretreat_q_table.tsv') %>%
	mutate(data = 'viremic pre-treatment PLHIV')
treat_q_dat = read_tsv('models/treat_q_table.tsv')%>%
	mutate(data = 'viremic treatment-experienced PLHIV')
all_q_dat = read_tsv('models/all_q_table.tsv')%>%
	mutate(data = 'viremic PLHIV')

pretreat_q_dat = pretreat_q_dat %>% 
	mutate(data = 'viremic pre-treatment PLHIV') %>%
	select(data, class, `Quasi Lik`, corr, params, QIC, CIC) %>% 
		pivot_longer(-c(data, class, corr, params)) %>%
		mutate(value = round(value, 3)) %>%
		arrange(class) %>%
		pivot_wider(names_from=c(name, class), values_from=value) %>%
	arrange(by=order[corr])

treat_q_dat = treat_q_dat %>% 
	mutate(data = 'viremic treatment-experienced PLHIV') %>%
	select(data, class, `Quasi Lik`, corr, params,  QIC, CIC) %>% 
		pivot_longer(-c(data, class, corr, params)) %>%
		mutate(value = round(value, 3)) %>%
		arrange(class) %>%
		pivot_wider(names_from=c(name, class), values_from=value) %>%
	arrange(by=order[corr])

all_q_dat = all_q_dat %>%
	mutate(data = 'viremic PLHIV') %>%
		select(data, class, `Quasi Lik`, corr, params,  QIC, CIC) %>% 
			pivot_longer(-c(data, class, corr, params)) %>%
			mutate(value = round(value, 3)) %>%
			arrange(class) %>%
			pivot_wider(names_from=c(name, class), values_from=value) %>%
		arrange(by=order[corr])

q_dat = bind_rows(pretreat_q_dat, treat_q_dat, all_q_dat) 

write_tsv(q_dat, 'tables/table_s9.tsv')