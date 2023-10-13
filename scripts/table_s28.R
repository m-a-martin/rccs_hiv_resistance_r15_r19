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

comm_type_prev_file = 'models/all_amongPar_prev_pred_stratified.tsv'
comm_type_v_prev_file = 'models/HIV_viremic_prev_pred_stratified.tsv'
class_prev_pred_file = 'models/all_amongPLWHIV_prev_pred_link.tsv'


comm_type_prev = read_tsv(comm_type_prev_file) %>% filter(round == 18) %>%
	mutate(comm_num = s_var) %>%
	select(class, comm_num, fit, lwr, upr)

comm_type_v_prev = read_tsv(comm_type_v_prev_file) %>% 
	filter(round == 18) %>%
	mutate(
		comm_num = comm_type,
		class = 'viremia') %>%
	select(class, comm_num, fit, lwr, upr)


class_prev = read_tsv(class_prev_pred_file) %>% filter(round == 18)

t4 = bind_rows(comm_type_prev, comm_type_v_prev)


# set up class_prev_list
prev = list()
for (c in unique(class_prev$class)){
	e = as.character(class_prev %>% filter(class == c))
	names(e) = colnames(class_prev)
	prev[[c]] = e
}


v_prev = (t4 %>% filter(class == 'viremia'))
for (c in unique(class_prev$class)){
	t4 = bind_rows(
		t4, 
		tibble(
			comm_num = v_prev$comm_num,
			fit = v_prev$fit * exp(as.numeric(prev[[c]][["fit"]])),
			lwr = v_prev$fit * exp(as.numeric(prev[[c]][["fit"]]) + qnorm(0.05/2)* as.numeric(prev[[c]][["se.fit"]])),
			upr =  v_prev$fit * exp(as.numeric(prev[[c]][["fit"]]) - qnorm(0.05/2)* as.numeric(prev[[c]][["se.fit"]])),
			class = paste(c, "_expect", sep="")))
}

# merge fit lwr and upr
# then pivot wider
t4 = t4 %>% 
	mutate_at(vars(fit, lwr, upr), ~round(.x*100, 2)) %>%
	unite("val", fit:lwr, sep=' (') %>%
	unite("val", val:upr, sep=', ') %>%
	mutate(
		val = paste(val, ')', sep='')) %>%
	pivot_wider(names_from=class, values_from=val)

# add spacer cols
for (c in unique(class_prev$class)){
	t4[paste(c, "_spacer", sep='')] = NA
}

# sort cols 
t4 = t4 %>% 
	select(sort_cols(colnames(t4))) 


# and save
write_tsv(t4, 'tables/table_s28.tsv', na='')
