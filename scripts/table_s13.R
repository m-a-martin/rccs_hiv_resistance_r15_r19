library(tidyverse)

source('scripts/utils.R')

dat_file = 'data/rakai_drug_resistance_categorized.tsv'
table_vars = c("age_cat", "comm_type", "sex")


hiv_prev_pred_file = 'models/HIV_prev_pred.tsv'
hiv_pred_file = 'models/HIV_format.tsv'

viremic_file = 'models/HIV_viremic_format.tsv'
viremic_pred_file = 'models/HIV_viremic_prev_pred.tsv'
pretreat_pred_file = 'models/HIV_pretreat_viremic_prev_pred.tsv'
pretreat_file = 'models/HIV_pretreat_format.tsv'
treat_file = 'models/HIV_treat_format.tsv'
treat_pred_file = 'models/HIV_treat_viremic_prev_pred.tsv'

	
t_c1 = read_tsv(hiv_prev_pred_file) %>%
	select(-se.fit) %>%
	mutate_at(vars(fit, lwr, upr), ~round(.x*100, 2)) %>%
	unite("Prev. (95% CI)", fit:lwr, sep=' (') %>%
	unite("Prev. (95% CI)", `Prev. (95% CI)`:upr, sep=', ') %>%
	mutate(`Prev. (95% CI)` = paste(`Prev. (95% CI)`, ')', sep='')) %>%
	rename(
		'HIV prev. (95% CI)' = `Prev. (95% CI)`,
		'survey round' = round)

t_c2 = read_tsv(hiv_pred_file) %>%
	rename(`survey round` = V1) %>%
	unite("Prev. ratio (95% CI)", RR:LCI, sep=' (') %>%
	unite("Prev. ratio (95% CI)", `Prev. ratio (95% CI)`: UCI, sep=', ') %>%
	mutate(
		`Prev. ratio (95% CI)` = paste(`Prev. ratio (95% CI)`, ')', sep=''),
		`Prev. ratio (95% CI)` = if_else(`survey round` == "(Intercept)", 'ref', `Prev. ratio (95% CI)`),
		P = if_else(`survey round` == "(Intercept)", 'ref', P),
		`survey round` = as.numeric(str_split(`survey round`, 'round', simplify=T)[,2]),
		`survey round` = if_else(is.na(`survey round`), 
			min(`survey round`, na.rm=TRUE)-1,
			`survey round`)) %>%
	rename('HIV prev. ratio (95% CI)' = `Prev. ratio (95% CI)`)

t_c3 = read_tsv(viremic_pred_file) %>%
	select(-se.fit) %>%
	mutate_at(vars(fit, lwr, upr), ~round(.x*100, 2)) %>%
	unite("Prev. (95% CI)", fit:lwr, sep=' (') %>%
	unite("Prev. (95% CI)", `Prev. (95% CI)`:upr, sep=', ') %>%
	mutate(`Prev. (95% CI)` = paste(`Prev. (95% CI)`, ')', sep='')) %>%
	rename(
		'Viremic HIV prev. (95% CI)' = `Prev. (95% CI)`,
		'survey round' = round)

t_c4 = read_tsv(viremic_file) %>%
	rename(`survey round` = V1) %>%
	unite("Prev. ratio (95% CI)", RR:LCI, sep=' (') %>%
	unite("Prev. ratio (95% CI)", `Prev. ratio (95% CI)`: UCI, sep=', ') %>%
	mutate(
		`Prev. ratio (95% CI)` = paste(`Prev. ratio (95% CI)`, ')', sep=''),
		`Prev. ratio (95% CI)` = if_else(`survey round` == "(Intercept)", 'ref', `Prev. ratio (95% CI)`),
		P = if_else(`survey round` == "(Intercept)", 'ref', P),
		`survey round` = as.numeric(str_split(`survey round`, 'round', simplify=T)[,2]),
		`survey round` = if_else(is.na(`survey round`), 
			min(`survey round`, na.rm=TRUE)-1,
			`survey round`)) %>%
	rename('HIV prev. ratio (95% CI)' = `Prev. ratio (95% CI)`)

t_c5 = read_tsv(pretreat_pred_file) %>%
	select(-se.fit) %>%
	mutate_at(vars(fit, lwr, upr), ~round(.x*100, 2)) %>%
	unite("Prev. (95% CI)", fit:lwr, sep=' (') %>%
	unite("Prev. (95% CI)", `Prev. (95% CI)`:upr, sep=', ') %>%
	mutate(`Prev. (95% CI)` = paste(`Prev. (95% CI)`, ')', sep='')) %>%
	rename(
		'Treatment-naive HIV prev. (95% CI)' = `Prev. (95% CI)`,
		'survey round' = round)

t_c6 = read_tsv(pretreat_file) %>%
	rename(`survey round` = V1) %>%
	unite("Prev. ratio (95% CI)", RR:LCI, sep=' (') %>%
	unite("Prev. ratio (95% CI)", `Prev. ratio (95% CI)`: UCI, sep=', ') %>%
	mutate(
		`Prev. ratio (95% CI)` = paste(`Prev. ratio (95% CI)`, ')', sep=''),
		`Prev. ratio (95% CI)` = if_else(`survey round` == "(Intercept)", 'ref', `Prev. ratio (95% CI)`),
		P = if_else(`survey round` == "(Intercept)", 'ref', P),
		`survey round` = as.numeric(str_split(`survey round`, 'round', simplify=T)[,2]),
		`survey round` = if_else(is.na(`survey round`), 
			min(`survey round`, na.rm=TRUE)-1,
			`survey round`)) %>%
	rename('Treatment-naive HIV prev. ratio (95% CI)' = `Prev. ratio (95% CI)`)

t_c7 = read_tsv(treat_pred_file) %>%
	select(-se.fit) %>%
	mutate_at(vars(fit, lwr, upr), ~round(.x*100, 2)) %>%
	unite("Prev. (95% CI)", fit:lwr, sep=' (') %>%
	unite("Prev. (95% CI)", `Prev. (95% CI)`:upr, sep=', ') %>%
	mutate(`Prev. (95% CI)` = paste(`Prev. (95% CI)`, ')', sep='')) %>%
	rename(
		'Treatment-experienced HIV prev. (95% CI)' = `Prev. (95% CI)`,
		'survey round' = round)

t_c8 = read_tsv(treat_file) %>%
	rename(`survey round` = V1) %>%
	unite("Prev. ratio (95% CI)", RR:LCI, sep=' (') %>%
	unite("Prev. ratio (95% CI)", `Prev. ratio (95% CI)`: UCI, sep=', ') %>%
	mutate(
		`Prev. ratio (95% CI)` = paste(`Prev. ratio (95% CI)`, ')', sep=''),
		`Prev. ratio (95% CI)` = if_else(`survey round` == "(Intercept)", 'ref', `Prev. ratio (95% CI)`),
		P = if_else(`survey round` == "(Intercept)", 'ref', P),
		`survey round` = as.numeric(str_split(`survey round`, 'round', simplify=T)[,2]),
		`survey round` = if_else(is.na(`survey round`), 
			min(`survey round`, na.rm=TRUE)-1,
			`survey round`)) %>%
	rename('Treatment-experienced HIV prev. ratio (95% CI)' = `Prev. ratio (95% CI)`)

blank_col = tibble(`survey round` = t_c2$`survey round`, ` `= rep('', nrow(t_c1)))

t = t_c1 %>%
	left_join(t_c2, by='survey round') %>% left_join(blank_col, by='survey round') %>%
	left_join(t_c3, by='survey round') %>% left_join(t_c4, by='survey round') %>% 
	left_join(blank_col, by='survey round') %>%
	left_join(t_c5, by='survey round') %>% left_join(t_c6, by='survey round') %>% 
	left_join(blank_col, by='survey round') %>%
	left_join(t_c7, by='survey round') %>%
	left_join(t_c8, by='survey round') 


write_tsv(t, 'tables/table_s13.tsv', na='..')



