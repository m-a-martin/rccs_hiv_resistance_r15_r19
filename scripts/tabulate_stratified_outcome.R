suppressMessages(library(tidyverse))
suppressMessages(library(argparser))
source('scripts/utils.R')


format_col = function(x1, x2){
		return(paste(paste(paste(x1, " (", sep=''),
			round((x1 / x2)*100, 2), sep=''),
			'%)', sep=''))
}

chisq_p = function(x, a, b){
	return(chisq.test(rbind(na.omit(x[[a]]), na.omit(x[[b]])))$p.value)
}


stratified_outcome_table = function(hiv_dr_cat, outcome_label, table_cats, table_vars){
	# add data columns for compound variables
	for (i in table_vars[grepl(':', table_vars)]){
		split_i = str_split(i, ":", simplify=TRUE)[1,]
		hiv_dr_cat["i"] = Reduce(paste, hiv_dr_cat[,split_i])
	}
	oa_dat = list()
	tab_dat = list()
	for (i in seq(1,length(table_cats$names))){
		par_col = paste(c('n_', table_cats$names[[i]], '_participants'), collapse='')
		outcome_col = paste(c('n_', table_cats$names[[i]], '_outcome'), collapse='')
		p_col = paste(c('n_', table_cats$names[[i]], '_p'), collapse='')
		spacer_col = paste(c('n_', table_cats$names[[i]], '_spacer'), collapse='')
		oa_dat[[table_cats$names[i]]] = hiv_dr_cat %>%
			filter(eval(parse(text=table_cats$get[[i]]))) %>%
			summarise(
				!!par_col := n(), 
				!!outcome_col := sum(outcome), 
				!!p_col := NA,
				!!spacer_col := NA,
				.groups='drop')
		for (i_var in table_vars){
			var_sum_dat = hiv_dr_cat %>%
						filter(eval(parse(text=table_cats$get[[i]]))) %>%
						rename(v = !!i_var) %>%
						group_by(v) %>%
						summarise(
							!!par_col := n(), 
							!!outcome_col := sum(outcome), 
							!!spacer_col := NA,
							.groups='drop')
			tab_dat[[table_cats$names[i]]] = 
				rbind(
					tab_dat[[table_cats$names[i]]],
					tibble(
						var = str_to_sentence(var_rename[i_var]),
						v = NA, 
						!!par_col := NA, 
						!!outcome_col := NA,
						!!p_col := chisq.test(var_sum_dat[[3]], p = var_sum_dat[[2]]/sum(var_sum_dat[[2]]), 
							simulate.p.value=TRUE)$p.value,
						#!!p_col := chisq.test(rbind(
						#		var_sum_dat[[2]] - var_sum_dat[[3]],
						#		var_sum_dat[[3]]))$p.value,
						!!spacer_col := NA),
					var_sum_dat %>%
						mutate('var' = str_to_sentence(var_rename[i_var]), !!p_col:=NA) %>%
						select(all_of(c('var', 'v', par_col, outcome_col, p_col, spacer_col))))
		}
	}
	# now column bind all of these
	sub_table =  Reduce(full_join, tab_dat) 
	table = rbind(
		Reduce(cbind, oa_dat) %>%
				mutate(var = 'Overall', v=NA),
			sub_table)
	table = table[colnames(sub_table)]
	# format columns
	for (i in 1:ncol(table)){
		if (str_detect(colnames(table)[i], 'outcome')){
			table[colnames(table)[i]] = 
				if_else(
					!is.na(table[[colnames(table)[i]]]),
					format_col(table[[colnames(table)[i]]], table[[colnames(table)[(i-1)]]]),
					NA)
			colnames(table)[i] = gsub('outcome', outcome_label, colnames(table)[i])
		}
		if (rev(str_split(colnames(table)[i], '_', simplify=TRUE))[1] == 'p'){
			table[colnames(table)[i]] = 
				if_else(!is.na(table[[colnames(table)[i]]]), 
					if_else(table[[colnames(table)[i]]] < 1E-4, 
						"<0.0001", 
						as.character(formatC(signif(table[[colnames(table)[i]]], 2), format='f', digits=4))), 
					NA)
		}
	}
	return(table)
}




#### --------------- ####
#### MAIN CODE BLOCK ####
#### --------------- ####
p <- arg_parser("round prev rr table")
p <- add_argument(p, "--out", help="output file name", 
	nargs=1)
p <- add_argument(p, "--dat", help="dat file path", 
	nargs=1)
p <- add_argument(p, "--filter", help="filter string for data", 
	nargs=1, default="TRUE")
p <- add_argument(p, "--outcome", help="data outcome", 
	nargs=1)
p <- add_argument(p, "--outcomeLabel", help="data outcomelabel", 
	nargs=1)
p <- add_argument(p, "--tableVars", help="table variables", 
	nargs=Inf)
p <- add_argument(p, "--tableCats", help="table categories", 
	nargs=Inf)
p <- add_argument(p, "--tableCatLabels", help="table categories labels", 
	nargs=Inf)
args <- parse_args(p)

#args$dat = 'data/rakai_drug_resistance_categorized.tsv'
#args$tableVars = c("age_cat", "comm_type", "sex", "vl_cat")
#args$outcome = 'valid_dr_dat == TRUE'
#args$outcomeLabel = 'sequenced'

# read in data and filter out bad rounds and age categories
hiv_dr_cat = read_tsv(args$dat, show_col_types=FALSE)

hiv_dr_cat = hiv_dr_cat %>%
	filter(eval(parse(text=args$filter)))
# add outcome column
hiv_dr_cat = hiv_dr_cat %>% mutate(outcome = eval(parse(text=args$outcome)))


table_cats = list(
	names=args$tableCatLabels, 
	get=args$tableCats)

t = stratified_outcome_table(hiv_dr_cat, args$outcomeLabel, table_cats, args$tableVars)

if (!dir.exists('tables')){dir.create('tables')}

write_tsv(t, 
	paste(c('tables/', args$out, '.tsv'), collapse=''), 
	na='')