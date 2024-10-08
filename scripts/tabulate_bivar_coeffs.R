suppressMessages(library(tidyverse))
suppressMessages(library(argparser))
source('scripts/utils.R')


sort_cols = function(x, class_order){
	type_order = c("coeff" = 1, "P" = 2, "spacer" = 3)
	non_round_x = x[x!= "s" & x != 'var']
	non_round_x_split = str_split(non_round_x, "_")
	s1 = class_order[unlist(lapply(non_round_x_split, function(x){x[1]}))]
	s2 = type_order[unlist(lapply(non_round_x_split, function(x){paste(x[(length(x))], collapse='_')}))]
	scores = tibble(val = non_round_x, s1 = s1, s2 = s2) %>%
		arrange(s1, s2)
	return(c("s", "var", scores$val))
}


#### --------------- ####
#### MAIN CODE BLOCK ####
#### --------------- ####
p <- arg_parser("bivariate coeffs table")
p <- add_argument(p, "--out", help="output file name", 
	nargs=1)
p <- add_argument(p, "--coeffs", help="estimated prevalence file path", 
	nargs=1)
p <- add_argument(p, "--classOrder", help="estimated rr file path", 
	nargs=Inf)
p <- add_argument(p, "--nOut", help="estimated rr file path", 
	nargs=1, default="all")
args <- parse_args(p)


if (all(is.na(args$classOrder))){
	prev = read_tsv(args$prev, show_col_types=FALSE) %>% filter(round == args$endRound) %>%
		arrange(-fit)
	class_order = setNames(seq(1,nrow(prev)), prev$class)
}else{
	class_order = setNames(seq(1,length(args$classOrder)), args$classOrder)
}


#args$classOrder = c('nnrti', 'nrti', 'pi')
#args$coeffs = 'models/dr_amongPar_bivar_coeffs.tsv'

coeffs = read_tsv(args$coeffs, show_col_types=FALSE)
tabulated = coeffs %>% 
	select(-corr) %>%
	mutate(
		across(RR:UCI, ~format_digit(.x)),
		P = if_else(P < 1E-4, "<0.0001", as.character(signif(P, 2))),
				P = if_else(`var` == "(Intercept)", 'ref', P),) %>%
	unite("coeff", RR:LCI, sep=' (') %>%
	unite("coeff", coeff:UCI, sep=', ')  %>%
	mutate(coeff = paste(coeff, ")", sep='')) %>%
	pivot_longer(-c(s, var, class)) %>%
	unite("name", c(class, name), sep="_") %>%
	pivot_wider(names_from='name', values_from='value')


# add blank columns
for (class in unique(read_tsv(args$coeffs, show_col_types=FALSE)$class)){
	tabulated[paste(class, 'spacer', sep='_')] = " "
}

# sort
tabulated = tabulated %>% select(sort_cols(colnames(tabulated), class_order)) %>%
	arrange(s)


write_tsv(tabulated, 
	paste(c('tables/', args$out, '.tsv'), collapse=''), 
	na='')

