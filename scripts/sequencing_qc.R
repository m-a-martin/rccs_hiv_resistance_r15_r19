suppressMessages(library(tidyverse))


get_cols = c('readnum_hiv', 'mapped_num', 'duprate', 'insertsize_05', 'insertsize_median', 'insertsize_95',
	'length_strict')

valid_dr_dat = read_tsv("data/rakai_drug_resistance_formatted.tsv", show_col_types=FALSE) %>%
	select(all_of(c('study_id', 'round', 'sequence_id', 'valid_dr_dat', get_cols))) %>%
	filter(valid_dr_dat == TRUE) %>%
	pivot_longer(-c('study_id', 'round', 'sequence_id', 'valid_dr_dat'))

# now summarise
probs =	c(0.025, 0.25, 0.5, 0.75, 0.975)
dr_cat_sum = bind_rows(
	valid_dr_dat %>% 
		filter(value != "NULL") %>%
		mutate(
			value = as.numeric(value),
			value = if_else(name == 'duprate', value*100, value)) %>%
		group_by(name) %>%
		summarise(as_tibble_row(quantile(value, probs=probs))) %>%
		mutate_at(paste(probs*100, '%', sep=''), ~round(.x, 2)),
	tibble(name='insert_size')) %>%
	arrange(
		!(name == 'readnum_hiv'),
		!(name == 'mapped_num'),
		!(name == 'insert_size'),
		!(name == 'insertsize_05'),
		!(name == 'insertsize_median'),
		!(name == 'insertsize_95'),
		!(name == 'duprate'),
		!(name == 'length_strict'))

write_tsv(dr_cat_sum, 'tables/seq_qc_among_valid.tsv', na='')