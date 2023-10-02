library(tidyverse)
source('scripts/utils.R')


table_s2 = bind_rows(
	as_tibble(cbind(nnrti_drugs, nnrti_cols)) %>% 
		rename(
			`generic name` = nnrti_drugs,
			`abbreviation` = nnrti_cols) %>%
		mutate(class = 'nnrti'),
	as_tibble(cbind(nrti_drugs, nrti_cols)) %>% 
		rename(
			`generic name` = nrti_drugs,
			`abbreviation` = nrti_cols) %>%
		mutate(class = 'nrti'),
	as_tibble(cbind(pi_drugs, pi_cols)) %>% 
		rename(
			`generic name` = pi_drugs,
			`abbreviation` = pi_cols) %>%
		mutate(class = 'pi'),
	as_tibble(cbind(insti_drugs, insti_cols)) %>% 
		rename(
			`generic name` = insti_drugs,
			`abbreviation` = insti_cols) %>%
		mutate(class = 'insti'))

write_tsv(table_s2, 'tables/table_s3.tsv')