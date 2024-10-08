suppressMessages(library(tidyverse))
suppressMessages(library(argparser))
source('scripts/utils.R')


format_p = function(x){
	return(if_else(x < 1E-4, "<0.0001", as.character(round(x,4))))
}


time_stratified_table = function(hiv_dr_cat, table_vars){
	sum_dat = list()
	for (i_var in table_vars){
		i_dat = hiv_dr_cat %>%
					rename(s_var = !!i_var) %>%
					filter(!is.na(s_var)) %>%
					mutate(s_var = if_else(
						s_var == 'intermediate' | s_var == 'high', 
						'intermediate/high', 
						s_var))
		oa_dat = bind_rows(
				i_dat %>%
					summarise(n=n()) %>%
					mutate(round = 'all', s_var='all'),
				i_dat %>%
					group_by(round) %>%
					summarise(n=n(), .groups='drop') %>%
					mutate(round = as.character(round),
						s_var = 'all')) 
		var_sum_dat = bind_rows(
			i_dat %>%
				group_by(s_var) %>%
				summarise(n=n(), .groups='drop') %>%
				mutate(round = 'all'),
			i_dat %>%
				group_by(round, s_var) %>%
				summarise(n=n(), .groups='drop') %>%
				mutate(round = as.character(round)))
		# calculate p_values
		#p_val_dat = var_sum_dat %>% filter(round != 'all') %>%
		#	pivot_wider(names_from=round, values_from=n, values_fill=0)
		#p_val = format_p(
		#	chisq.test(apply(as.matrix(p_val_dat)[,2:ncol(p_val_dat)], 2,as.numeric))$p.value)
		# calculate proportions
		var_sum_fmt = var_sum_dat %>% 
			left_join(oa_dat %>% rename(n_all=n) %>% select(-s_var), by='round') %>%
			mutate(p = n/n_all,
				val = paste(paste(paste(n, ' (', sep=''), round(100*p,2), sep=''), '%)', sep='')) %>%
			select(-n, -n_all, -p)
		t = bind_rows(
			oa_dat %>% rename(val=n) %>% 
				mutate(val = as.character(val)),
			var_sum_fmt %>%
				mutate(val = as.character(val))) %>%
			pivot_wider(names_from=round, values_from=val, values_fill = '0 (0%)') 
		#t[['p-val']] = c(p_val, rep(NA, nrow(t)-1))
		t = t %>% arrange(
				s_var != 'all',
				s_var != 'susceptible',
				s_var != 'low',
				s_var != 'intermediate/high') %>%
			mutate(s = i_var)
		sum_dat[[i_var]] = t
	}
	out = bind_rows(sum_dat)
	return(out)
}




#### --------------- ####
#### MAIN CODE BLOCK ####
#### --------------- ####
p <- arg_parser("round prev rr table")
p <- add_argument(p, "--out", help="output file name", 
	nargs=1)
p <- add_argument(p, "--dat", help="dat file path", 
	nargs=1)
args <- parse_args(p)


#args$dat = 'data/rakai_drug_resistance_categorized.tsv'
args$tableVars = c("insti", "nnrti", "nrti", "pi")

hiv_dr_cat = read_tsv(args$dat, show_col_types=FALSE) %>%
	filter(finalhiv == 'P' & valid_dr_dat == TRUE & viremic == TRUE)

t = time_stratified_table(hiv_dr_cat, args$tableVars) %>% 
	mutate(type = 'all')

pretreat_t = time_stratified_table(hiv_dr_cat %>% filter(pre_treatment == TRUE), args$tableVars) %>% 
	mutate(type = 'pretreat')

treat_t = time_stratified_table(hiv_dr_cat %>% filter(pre_treatment == FALSE), args$tableVars) %>% 
	mutate(type = 'treat')

t = bind_rows(t, pretreat_t, treat_t)

if (!dir.exists('tables')){dir.create('tables')}

write_tsv(t, 
	paste(c('tables/', args$out, '.tsv'), collapse=''), 
	na='')


