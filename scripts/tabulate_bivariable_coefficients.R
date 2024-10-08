suppressMessages(library(tidyverse))
suppressMessages(library(argparser))
source('scripts/utils.R')

sort_cols = function(x, class_order){
		type_order = c("n" = 1, "y" = 2, "Prev. % (95% CI)" = 3, "Prev. ratio (95% CI)" = 4, "p-value" = 5, "spacer" = 6)
		non_round_x = x[x!= "round" & x != 's' & x != 's_var']
		non_round_x_split = str_split(non_round_x, "_")
		s1 = class_order[unlist(lapply(non_round_x_split, function(x){paste(x[1:(length(x)-1)], collapse='_')}))]
		s2 = type_order[unlist(lapply(non_round_x_split, function(x){paste(x[(length(x))], collapse='_')}))]
		scores = tibble(val = non_round_x, s1 = s1, s2 = s2) %>%
			arrange(s1, s2)
		return(c("round", 's', 's_var', scores$val))
	}
	
#### --------------- ####
#### MAIN CODE BLOCK ####
#### --------------- ####
p <- arg_parser("round prev rr table")
p <- add_argument(p, "--out", help="output file name", 
	nargs=1)
p <- add_argument(p, "--coeffs", help="estimated prevalence file path", 
	nargs=1)
p <- add_argument(p, "--classOrder", help="estimated rr file path", 
	nargs=Inf)
args <- parse_args(p)

#args$coeffs = 'models/dr_amongPretreat_bivar_coeffs.tsv'
#args$rr = 'models/dr_amongPretreat_rr_stratified.tsv'
#args$prev = 'models/dr_amongPretreat_prev_pred_stratified.tsvv'
#args$class_order =c('plhiv', 'v_plhiv', 'pt_v_plhiv', 't_v_plhiv')
#args$classOrder = c('nnrti', 'nrti', 'pi', 'nnrti_nrti', 'nnrti_pi', 'nrti_pi', 'nnrti_nrti_pi')
#pred_file = args$prev
#rr_file = args$rr
class_order = setNames(seq(1,length(args$classOrder)), args$classOrder)

coeffs = read_tsv(args$coeffs, show_col_types=FALSE) %>%
	filter(!grepl("round", var) & var != "(Intercept)")

coeffs = coeffs %>%
		select(-any_of(c("SE", 'corr', ''))) %>%
		mutate(across(RR:UCI, ~format_digit(.x))) %>%
		unite("val", RR:LCI, sep=' (') %>%
		unite("val", val: UCI, sep=', ') %>%
		mutate(
			val = paste(val, ')', sep=''),
			P = if_else(P < 1E-4, "<0.0001", as.character(formatC(signif(P, 2), format='f', digits=4)))) %>%
		pivot_longer(c(val, P)) %>%
		mutate(
			type = if_else(name == 'val',
				paste(class, "Prev. ratio (95% CI)", sep='_'),
				paste(class, "p-value", sep='_'))) %>%
		select(-any_of(c('name', 'class'))) %>%
		pivot_wider(names_from='type', values_from='value')

t = class_strat_round_prev_rr_table(args$prev, args$rr, class_order)

if (!dir.exists('tables')){dir.create('tables')}

write_tsv(t, 
	paste(c('tables/', args$out, '.tsv'), collapse=''), 
	na='..')