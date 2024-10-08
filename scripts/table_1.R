suppressMessages(library(tidyverse))
source('scripts/utils.R')

format_col = function(x1, x2){
		return(paste(paste(paste(x1, " (", sep=''),
			round((x1 / x2)*100, 2), sep=''),
			'%)', sep=''))
}


dat_file = 'data/rakai_drug_resistance_categorized.tsv'
table_vars = c("age_cat", "comm_type", "sex")

# read in data and filter out bad rounds and age categories
hiv_dr_cat = read_tsv(dat_file, show_col_types=FALSE) %>% 
	mutate(round = as.character(round),
		viremic = if_else(round == 15 & finalhiv == 'P' & pre_treatment == FALSE, NA, viremic))

t1 = bind_rows(
	bind_rows(
			hiv_dr_cat %>% 
				mutate(round = 'overall'), 
			hiv_dr_cat) %>%
		group_by(round) %>%
		summarise(
			all = n(),
			all_plhiv = sum(finalhiv == 'P', na.rm=TRUE), 
			viremic_plhiv = sum(finalhiv == 'P' & viremic == TRUE, na.rm=TRUE), 
			pt_viremic_plhiv = sum(finalhiv == 'P' & viremic == TRUE & pre_treatment == TRUE, na.rm=TRUE),
			t_viremic_plhiv = sum(finalhiv == 'P' & viremic == TRUE & pre_treatment == FALSE, na.rm=TRUE)),
	tibble(round = "Survey round")) %>%
	arrange(round != 'overall', round != 'Survey round') %>%
	mutate(viremic_plhiv = if_else(round != 15, viremic_plhiv, NA),
		t_viremic_plhiv = if_else(round != 15, t_viremic_plhiv, NA))


t1[ncol(t1)] = 
		c(as.character(t1[ncol(t1)][1,]),
		if_else(
			!is.na(t1[ncol(t1)][2:nrow(t1),]),
			format_col(t1[[ncol(t1)]][2:nrow(t1)], t1[[ncol(t1)-2]][2:nrow(t1)]),
			NA))

for (i in (ncol(t1)-1):4){
	t1[colnames(t1)[i]] = 
		gsub(' \\(NA%\\)', '', 
			c(as.character(t1[[colnames(t1)[i]]][1]),
				if_else(
					!is.na(t1[[colnames(t1)[i]]][2:nrow(t1)]),
					format_col(t1[[colnames(t1)[i]]][2:nrow(t1)], t1[[colnames(t1)[i-1]]][2:nrow(t1)]),
					NA)))
}

for (i in (3:3)){
	t1[colnames(t1)[i]] = 
		gsub(' \\(NA%\\)', '', if_else(
			!is.na(t1[[colnames(t1)[i]]]),
			format_col(t1[[colnames(t1)[i]]], t1[[colnames(t1)[(i-1)]]]),
			NA))

}


if (!dir.exists('tables')){dir.create('tables')}

write_tsv(t1, "tables/table_1.tsv", na='')



