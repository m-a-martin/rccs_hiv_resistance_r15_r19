library(tidyverse)
source('scripts/utils.R')



dat_file = 'data/rakai_drug_resistance_categorized.tsv'
table_vars = c("age_cat", "comm_type", "sex", "vl_cat")

# read in data and filter out bad rounds and age categories
hiv_dr_cat = read_tsv(dat_file) %>% 
	filter(
		age_cat != '(50, 100]' & 
		round > 14 &
		viremic &
		finalhiv == 'P')

insti = hiv_dr_cat %>% filter(!is.na(insti)) %>%
	summarise(class = 'insti',
		n = n(), 
		n_resistant = sum(insti == 'intermediate/high'))

nnrti = hiv_dr_cat %>% filter(!is.na(nnrti)) %>%
	summarise(class = 'nnrti',
		n = n(), 
		n_resistant = sum(nnrti == 'intermediate/high'))


nrti = hiv_dr_cat %>% filter(!is.na(nrti)) %>%
	summarise(class = 'nrti',
		n = n(), 
		n_resistant = sum(nrti == 'intermediate/high'))

pi = hiv_dr_cat %>% filter(!is.na(pi)) %>%
	summarise(class = 'pi',
		n = n(), 
		n_resistant = sum(pi == 'intermediate/high'))

t = bind_rows(insti, nnrti, nrti, pi)
t$n_resistant = format_col(t$n_resistant, t$n)


write_tsv(t, 'tables/table_s18.tsv')
