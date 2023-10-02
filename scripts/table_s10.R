library(tidyverse)


order=c('independence'=1, 'exchangeable'=2, 'ar1'=3, 'unstructured'=4)
pretreat_m_dat = read_tsv('models/pretreat_m_table.tsv') 

m = pretreat_m_dat %>%
	select(class, round, corr, coefficient, std.err) %>%
	pivot_longer(-c(class, round, corr)) %>%
	mutate(value = round(value, 3)) %>%
	pivot_wider(names_from=c(corr, name), values_from=value)

write_tsv(m, 'tables/table_s10.tsv')