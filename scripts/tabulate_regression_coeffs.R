suppressMessages(library(tidyverse))
suppressMessages(library(argparser))
suppressMessages(source('scripts/utils.R'))


sort_cols = function(x, excl=c()){
	round_order = c("nnrti" = 1, "nrti" = 2, "pi" = 3, "NNRTI-NRTI-PI" = 4)
	type_order = c("Coefficient (95% CI)" = 1, "p-value" = 3, "spacer"=4)
	non_excl_x = x[!(x %in% excl)]
	s1 = round_order[(str_split(non_excl_x, "_", simplify=T)[,1])]
	s2 = type_order[str_split(non_excl_x, "_", simplify=T)[,2]]
	scores = tibble(val = non_excl_x, s1 = s1, s2 = s2) %>%
		arrange(s1, s2)
	return(c(excl, scores$val))
}


tabulate_regression_coeffs = function(m){
	d = read_tsv(m, show_col_types=FALSE)
	d_table = d %>%
		mutate(RR = round(RR, 2),
			LCI = round(LCI, 2),
			UCI = round(UCI, 2)) %>%
	rename(var = `...1`) %>%
	unite("val", RR:LCI, sep=' (') %>%
	unite("val", val: UCI, sep=', ') %>%
	mutate(
		val = paste(val, ')', sep=''),
		#val = if_else(`var` == "(Intercept)", 'ref', val),
		P = if_else(P < 1E-4, "<0.0001", as.character(signif(P, 2))),
		var_type = case_when(
			str_detect(var, 'comm_type') ~ 'Community type',
			str_detect(var, 'sex') ~ 'Sex',
			str_detect(var, 'age_cat') ~ 'Age category',
			str_detect(var, 'log10vl') ~ 'Viral load',
			str_detect(var, 'missing_vl') ~ 'Viral load',
			str_detect(var, 'round') ~ 'Survey round',
			str_detect(var, 'pre_treatment') ~ 'Treatment')) %>%
	rename(Variable = var, `Coefficient (95% CI)` = val, `p-value` = P)
	# todo infer this somehow?
	ref_cat = c('Sex'='sexF', 'Age category'='age_cat(14,24]', 
		'Community type' = 'comm_typeAgrarian')
	if (grepl('pre_treatment', d[,1])){
		ref_cat = c(ref_cat, 'Treatment'='pre_treatmentFALSE')
	}
	ref_table = tibble(var_type = names(ref_cat), Variable = ref_cat, `Coefficient (95% CI)` = 'ref', 
		`p-value` = 'ref')
	header_table = tibble(var_type = names(ref_cat), Variable = NA, `Coefficient (95% CI)` = NA,
		`p-value` = NA)
	m_label = str_split(m, "/", simplify=TRUE)
	m_label = str_split(m_label[,ncol(m_label)], "\\.", simplify=TRUE)[,1]
	d_table = bind_rows(d_table, ref_table, header_table) %>% 
		arrange((Variable != '(Intercept)' | is.na(Variable)), 
			var_type, !is.na(Variable), `Coefficient (95% CI)` != 'ref') %>%
		select(var_type, Variable, `Coefficient (95% CI)`, `p-value`) %>%
		rename(
			!!paste(m_label, "Coefficient (95% CI)", sep='_') := `Coefficient (95% CI)`,
			!!paste(m_label, "p-value", sep='_') := `p-value`)
	return(d_table)
}


#### --------------- ####
#### MAIN CODE BLOCK ####
#### --------------- ####
p <- arg_parser("round prev rr table")
p <- add_argument(p, "--out", help="output file name", 
	nargs=1)
p <- add_argument(p, "--dat", help="dat file path", 
	nargs=1)
p <- add_argument(p, "--models", help="regression models to summarize",
	nargs=Inf)
args <- parse_args(p)

#args$models = c(
#	'models/all_R15_nnrti_weights.tsv',
#	'models/all_R15_nrti_weights.tsv',
#	'models/all_R15_pi_weights.tsv',
#	'models/all_R15_nnrti_nrti_pi_weights.tsv',
#	'models/all_R15_insti_nnrti_nrti_pi_weights.tsv')

out = list()
for (m in args$models){
	out[[m]] = tabulate_regression_coeffs(m)
	out[[m]] = bind_cols(out[[m]], tibble(spacer = rep(NA, nrow(out[[m]]))))
}

out = reduce(out, inner_join, c('var_type', 'Variable'))

write_tsv(out, paste(c('tables/', args$out, '.tsv'), collapse=''), 
	na='')
