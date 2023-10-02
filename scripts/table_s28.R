library(tidyverse)
source('scripts/utils.R')


sort_cols = function(x){
		class_order = c("viremia" = 1, "nnrti" = 2, "nrti" = 3, "pi" = 4)
		type_order = c("expect" = 2, "spacer" = 3)
		non_comm_x = x[x!= "comm_num"]
		s1 = class_order[str_split(non_comm_x, "_", simplify=T)[,1]]
		s2 = type_order[str_split(non_comm_x, "_", simplify=T)[,2]]
		scores = tibble(
			val = non_comm_x, s1 = s1, s2 = s2) %>%
			mutate(s2 = replace_na(s2, 1)) %>%
			arrange(s1, s2)
		return(c("comm_num", scores$val))
	}

dat_file = 'data/rakai_drug_resistance_categorized.tsv'
comm_num_prev_file = 'models/all_amongPar_comm_num_R18_prev_pred.tsv'
comm_num_v_prev_file = 'models/HIV_comm_num_viremic_prev_pred.tsv'
class_prev_pred_file = 'models/all_amongPLWHIV_prev_pred_link.tsv'


comm_num_type = read_tsv(dat_file) %>% select(comm_type, comm_num) %>% unique() %>% 
	mutate(comm_num = as.character(comm_num))

comm_num_prev = read_tsv(comm_num_prev_file) %>% filter(round == 18) %>% 
	mutate(comm_num = as.character(comm_num)) %>%
	select(class, comm_num, fit, lwr, upr)


comm_num_v_prev = read_tsv(comm_num_v_prev_file) %>% 
	filter(round == 18) %>% 
	mutate(comm_num = as.character(comm_num),
	class = 'viremia') %>%
	select(class, comm_num, fit, lwr, upr)


class_prev = read_tsv(class_prev_pred_file) %>% filter(round == 18)

ts11 = bind_rows(comm_num_prev, comm_num_v_prev)


# set up class_prev_list
prev = list()
for (c in unique(class_prev$class)){
	e = as.character(class_prev %>% filter(class == c))
	names(e) = colnames(class_prev)
	prev[[c]] = e
}


v_prev = (ts11 %>% filter(class == 'viremia'))
for (c in unique(class_prev$class)){
	ts11 = bind_rows(
		ts11, 
		tibble(
			comm_num = v_prev$comm_num,
			fit = v_prev$fit * exp(as.numeric(prev[[c]][["fit"]])),
			lwr = v_prev$fit * exp(as.numeric(prev[[c]][["fit"]]) + qnorm(0.05/2)* as.numeric(prev[[c]][["se.fit"]])),
			upr =  v_prev$fit * exp(as.numeric(prev[[c]][["fit"]]) - qnorm(0.05/2)* as.numeric(prev[[c]][["se.fit"]])),
			class = paste(c, "_expect", sep="")))
}

# merge fit lwr and upr
# then pivot wider
ts11 = ts11 %>% 
	mutate_at(vars(fit, lwr, upr), ~round(.x*100, 2)) %>%
	unite("val", fit:lwr, sep=' (') %>%
	unite("val", val:upr, sep=', ') %>%
	mutate(
		val = paste(val, ')', sep='')) %>%
	pivot_wider(names_from=class, values_from=val)

# add spacer cols
for (c in unique(class_prev$class)){
	ts11[paste(c, "_spacer", sep='')] = NA
}

# sort cols 
ts11 = ts11 %>% 
	select(sort_cols(colnames(ts11))) 

# finally, add indicator of community type 
ts11 = comm_num_type %>% 
	right_join(ts11) %>%
	mutate(
		comm_type = if_else(is.na(comm_type), comm_num, comm_type)) %>% 
	arrange(!(comm_num == 'all'), desc(viremia))


# and save
write_tsv(ts11, 'tables/table_s28.tsv', na='')
