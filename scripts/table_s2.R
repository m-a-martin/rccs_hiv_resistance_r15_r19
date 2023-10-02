library(tidyverse)


format_col = function(x1, x2){
		return(paste(paste(paste(x1, " (", sep=''),
			round((x1 / x2)*100, 2), sep=''),
			'%)', sep=''))
}


hiv_dr_cat = read_tsv('data/rakai_drug_resistance_categorized.tsv') 

round_tech = hiv_dr_cat %>%
	mutate(round = as.character(round)) %>%
	filter(valid_dr_dat) %>%
	group_by(round) %>%
	summarise(
		N = n(),
		PANGEA1 = format_col(sum(sequenced == 'PANGEA1'), N),
		PANGEA2 = format_col(sum(sequenced == 'PANGEA2'), N))

oa_tech = hiv_dr_cat %>%
	filter(valid_dr_dat) %>%
	summarise(
		round = 'overall',
		N = n(),
		PANGEA1 = format_col(sum(sequenced == 'PANGEA1'), N),
		PANGEA2 = format_col(sum(sequenced == 'PANGEA2'), N))

t = bind_rows(
	oa_tech, tibble(round = 'round', PANGEA1 = '', PANGEA2 = ''),
	round_tech)

write_tsv(t, 'tables/table_s2.tsv')