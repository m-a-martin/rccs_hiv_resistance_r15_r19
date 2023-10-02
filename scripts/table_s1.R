library(tidyverse)
library(haven)
source('scripts/utils.R')

rccs_metadata_backup = read_stata('data/RCCSdata_R001_R019_VOIs.dta')
rccs_metadata = read_format_rccs_dat('data/RCCSdata_R001_R019_VOIs.dta')
hiv_dr_cat = read_tsv('data/rakai_drug_resistance_categorized.tsv')

# get 99th percentile bounds
round_bounds = rccs_metadata %>% 
	filter(round != 15.1) %>%
	group_by(round) %>%
	summarise(
		`survey start date` = quantile(int_date, c(5E-4), type=1, na.rm=TRUE),
		`survey median date` = quantile(int_date, c(0.5), type=1, na.rm=TRUE),
		`survey end date` = quantile(int_date, c(1-5E-4), type=1, na.rm=TRUE))

write_tsv(round_bounds, 'tables/table_s1.tsv')
