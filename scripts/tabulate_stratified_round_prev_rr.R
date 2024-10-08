suppressMessages(library(tidyverse))
suppressMessages(library(argparser))
source('scripts/utils.R')


class_strat_round_prev_rr_table = function(pred_file, rr_file, class_order){
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
	all_pred = read_tsv(pred_file, show_col_types=FALSE) %>%
		mutate(round = as.character(round), pre_treatment = as.character(pre_treatment))
	# format for table
	pred = all_pred %>%
		select(-any_of(c("n_outcome_obs", "se.fit", "df", "n", "y"))) %>%
		filter(!is.na(fit)) %>%
		mutate(across(fit:upr, ~format_digit(.x*100))) %>%
		unite("val", fit:lwr, sep=' (') %>%
		unite("val", val:upr, sep=', ') %>%
		mutate(
			val = paste(val, ')', sep=''),
			type = paste(class, "Prev. % (95% CI)", sep='_')) %>%
		select(-class,  -corr)
	if (any(colnames(all_pred) == 'n')){
		pred = bind_rows(
			pred,
			all_pred %>% 
				mutate(
					type = paste(class, "n", sep='_'), val=as.character(n)) %>% 
				select(-any_of(c("n_outcome_obs", "se.fit", "df", "n", "y", "corr", "fit", "lwr", "upr", "class"))),
			all_pred %>% 
				mutate(type = paste(class, "y", sep='_'), val=format_col(y,n)) %>% 
				select(-any_of(c("n_outcome_obs", "se.fit", "df", "n", "y", "corr", "fit", "lwr", "upr", "class"))))
	}
	#pred["s"] = colnames(pred)[!(colnames(pred) %in% c('round', 'val', 'type'))][
	#	which(!is.na(pred[!(colnames(pred) %in% c('round', 'val', 'type'))]), 
	#		arr.ind=TRUE)[,2]]
	pred["s_var"] = apply(pred[,!(colnames(pred) %in% c('round', 'val', 'type', 's'))] %>% 
		replace(is.na(.), ""), 1, paste, collapse='')
	pred = pred %>% 
		mutate(
			round = if_else(is.na(round), s_var, round),
			s_var = if_else(s_var == round, "", s_var)) %>%
		select(-any_of(colnames(pred)[!(colnames(pred) %in% 
			c('round', 'val', 'type', 's', 's_var'))]))
	all_cats = pred %>% select(round, s, s_var) %>% unique()
	rr = read_tsv(rr_file, show_col_types=FALSE) %>%
		select(-any_of(c("SE", 'corr', ''))) %>%
		mutate(across(RR:UCI, ~format_digit(.x))) %>%
		unite("val", RR:LCI, sep=' (') %>%
		unite("val", val: UCI, sep=', ') %>%
		filter(var != "(Intercept)") %>%
		mutate(
			val = paste(val, ')', sep=''),
			P = if_else(P < 1E-4, "<0.0001", as.character(formatC(signif(P, 2), format='f', digits=4))),
			round = gsub("round", "", var)) %>%
		pivot_longer(c(val, P)) %>%
		mutate(
			type = if_else(name == 'val',
				paste(class, "Prev. ratio (95% CI)", sep='_'),
				paste(class, "p-value", sep='_'))) %>%
		rename(val=value) %>%
		select(-any_of(c('class', 'var', 'name')))
	rr["s_var"] = apply(rr[,!(colnames(rr) %in% c('round', 'val', 'type', 's'))] %>% 
		mutate_all(~as.character(.)) %>%
		replace(is.na(.), ""), 1, paste, collapse='')
	rr = rr	%>%
		select(all_of(c('round', 'val', 'type', 's', 's_var')))
	t = bind_rows(pred, rr)
	t %>% print(n=250)
	t = t  %>% pivot_wider(names_from=type, values_from=val) %>% 	
		arrange(s, !is.na(suppressWarnings(as.numeric(round))), round, s_var) 
	# gotta be a tidy way to do this 
	for (i in colnames(t)[grepl('_Prev.',colnames(t)) | grepl('_RR.',colnames(t)) | grepl('p-value', colnames(t))]){
		t[[i]] = replace_na(t[[i]], "ref")
	}
	# add blank columns
	for (class in unique(read_tsv(pred_file, show_col_types=FALSE)$class)){
		t[paste(class, 'spacer', sep='_')] = " "
	}

	t = bind_rows(t, all_cats %>% select(-round) %>% unique())
	t[is.na(t$round),colnames(t)[4:ncol(t)]] = ""
	# sort
	t = t %>% select(sort_cols(colnames(t), class_order)) %>%
		arrange(s, s_var, !is.na(round), round)
	#%>%
	#	arrange(s, round, !is.na(s_var), s_var)
	return(t)
}

#### --------------- ####
#### MAIN CODE BLOCK ####
#### --------------- ####
p <- arg_parser("round prev rr table")
p <- add_argument(p, "--out", help="output file name", 
	nargs=1)
p <- add_argument(p, "--prev", help="estimated prevalence file path", 
	nargs=1)
p <- add_argument(p, "--rr", help="estimated rr file path", 
	nargs=1)
p <- add_argument(p, "--classOrder", help="estimated rr file path", 
	nargs=Inf)
args <- parse_args(p)

#args$rr = 'models/dr_amongPar_rr_stratified.tsv'
#args$prev = 'models/dr_amongPar_prev_pred_stratified.tsv'
#args$class_order =c('nnrti nrti pi')
#pred_file = args$prev
#rr_file = args$rr
class_order = setNames(seq(1,length(args$classOrder)), args$classOrder)

t = class_strat_round_prev_rr_table(args$prev, args$rr, class_order)

if (!dir.exists('tables')){dir.create('tables')}

write_tsv(t, 
	paste(c('tables/', args$out, '.tsv'), collapse=''), 
	na='..')


