suppressMessages(require(tidyverse))
suppressMessages(require(haven))
suppressMessages(require(readxl))
source('scripts/utils.R')

# impute round 15 viral loads
set.seed(1)
hiv_dr = hiv_dr %>%
	left_join(
		hiv_dr %>% filter(round == 15 & finalhiv == 'P' & vl_cat != 'missing' & pre_treatment == TRUE) %>%
			group_by(round, finalhiv, pre_treatment, valid_dr_dat) %>%
			summarise(p_viremic = sum(viremic_raw)/n(), .groups='drop') %>%
			mutate(vl_cat = 'missing'),
		by=c('round', 'finalhiv', 'pre_treatment', 'valid_dr_dat', 'vl_cat')) %>% 
	mutate(viremic = case_when(
		!is.na(viremic_raw) ~ viremic_raw,
		round== 15 & pre_treatment == TRUE & is.na(viremic_raw) ~ as.logical(rbinom(n(), 1, replace_na(p_viremic, 0))),
		round == 15 & pre_treatment == FALSE & is.na(viremic_raw) ~ NA))
