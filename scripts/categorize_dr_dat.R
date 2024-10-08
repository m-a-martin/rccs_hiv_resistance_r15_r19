suppressMessages(library(tidyverse))
suppressMessages(source('scripts/utils.R'))



#get_resistance = function(.data, cols, label){
#	.data = .data %>% mutate(!!label := case_when(
#		# if any of the columns are a 3 and 4 and 
#		# all of the columns are in S, 1, 2, 3, or 4 (i.e. not X) then call it intermediate/high
#		if_any(all_of(cols), function(x) x %in% c("3", "4")) & 
#			if_all(all_of(cols), function(x) x %in% c("P", "S", "0", "1", "2", "3", "4")) ~ 'intermediate/high',
#		# if all of the columns are an 'S', '1', or '2'
#		# then call it susceptible
#		if_all(all_of(cols), function(x) x %in% c("P", "S", "0", "1", "2")) ~ 'susceptible'))
#		# otherwise will return NA
#	return(.data)
#}

get_resistance = function(.data, cols, label){
	.data = .data %>% mutate(!!label := case_when(
		# if any of the columns are a 3 and 4 and 
		# all of the columns are in S, 1, 2, 3, or 4 (i.e. not X) then call it intermediate/high
		if_any(all_of(cols), function(x) x %in% c("4")) & 
			if_all(all_of(cols), function(x) x %in% c("P", "S", "0", "1", "2", "3", "4")) ~ 'high',
		if_any(all_of(cols), function(x) x %in% c("3")) & 
			if_all(all_of(cols), function(x) x %in% c("P", "S", "0", "1", "2", "3", "4")) ~ 'intermediate',
		if_any(all_of(cols), function(x) x %in% c("2")) & 
			if_all(all_of(cols), function(x) x %in% c("P", "S", "0", "1", "2", "3", "4")) ~ 'low',
		# if all of the columns are an 'S', '1', or '2'
		# then call it susceptible
		if_all(all_of(cols), function(x) x %in% c("P", "S", "0", "1", "2")) ~ 'susceptible'))
		# otherwise will return NA
	return(.data)
}


hiv_dr = read_tsv('data/rakai_drug_resistance_formatted.tsv', show_col_types=FALSE)


hiv_dr_cat = hiv_dr %>%
	get_resistance(insti_cols, "insti") %>%
	get_resistance(nnrti_cols, "nnrti") %>%
	get_resistance(nrti_cols, "nrti") %>%
	get_resistance(pi_cols, "pi") %>%
	select(-all_of(c(insti_cols, nnrti_cols, nrti_cols, pi_cols, insti_cols)))

write_tsv(hiv_dr_cat, 'data/rakai_drug_resistance_categorized.tsv')
