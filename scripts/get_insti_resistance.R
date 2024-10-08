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

insti_muts = hiv_muts %>% filter(substr(mut, 1, 2) == 'in')


insti_dr = hiv_dr_cat %>% select(study_id, round, pre_treatment, insti) %>%
	filter(insti == 'intermediate' | insti == 'high')

insti_dr = insti_dr %>% left_join(
	insti_muts %>% 
		select(mut, freq, reads, study_id, round), 
	by=c('study_id', 'round'))

insti_dr %>% arrange(pre_treatment, study_id, round, mut) %>% print(n=100)

nrow(insti_dr %>% filter(pre_treatment == TRUE) %>% select(study_id) %>% unique())
nrow(insti_dr %>% filter(pre_treatment == TRUE))
insti_dr %>% filter(pre_treatment == TRUE) %>% group_by(mut) %>% summarise(n=n())


nrow(insti_dr %>% filter(pre_treatment == FALSE) %>% select(study_id) %>% unique())
nrow(insti_dr %>% filter(pre_treatment == FALSE))
insti_dr %>% filter(pre_treatment == FALSE) %>% group_by(mut) %>% summarise(n=n())
