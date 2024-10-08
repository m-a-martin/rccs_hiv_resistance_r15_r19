suppressMessages(require(tidyverse))
suppressMessages(source('scripts/utils.R'))

# todo read in as arguments
dr_cat_file = 'data/rakai_drug_resistance_categorized.tsv'
dat_file = 'data/rakai_drug_resistance_mut_formatted.tsv'
w_file = 'models/all_sequencing_probs.tsv'

# read in mutation data
hiv_muts = read_tsv(dat_file, show_col_types=FALSE)
# read in data
hiv_dr_cat = read_tsv(dr_cat_file, show_col_types=FALSE)

# calc viremic participant-visits
totals = hiv_dr_cat %>% filter(viremic == TRUE & finalhiv == 'P' & !is.na(sequenced) & 
	!is.na(nnrti) & !is.na(nrti) & !is.na(pi) & !is.na(insti)) %>%
	group_by(pre_treatment) %>% summarise(n=n()) %>%
	mutate(
		type = if_else(pre_treatment == FALSE, 'treatment_experienced', 'pre_treatment'),
		label=paste(paste(type, '_n_', sep=''), n, sep=''),
		pre_treatment = as.character(pre_treatment))

totals = bind_rows(totals, 
	tibble(
		pre_treatment = 'all', 
		n = sum(totals$n),
		label = paste('total_n_', sum(totals$n), sep='')))

labels = setNames(totals$label, totals$pre_treatment)

# get all mut counts
mut_counts = hiv_muts %>% 
	filter(viremic == TRUE & finalhiv == 'P') %>%
	group_by(pre_treatment, mut) %>%
	summarise(n=n()) %>%
	mutate(pre_treatment = as.character(pre_treatment)) %>%
	pivot_wider(names_from='pre_treatment', values_from='n', values_fill=0) %>%
	mutate(all = `FALSE` + 	`TRUE`) %>%
	select(c('mut', 'all', 'TRUE', 'FALSE')) %>%
	arrange(-all, -`TRUE`, -`FALSE`)

colnames(mut_counts)[2:ncol(mut_counts)] = labels[colnames(mut_counts)[2:ncol(mut_counts)]]

write_tsv(mut_counts, 'tables/mut_counts.tsv')
