
gtheme = theme_bw() + 
	theme(text=element_text(color='#333333'),
			plot.title = element_text(hjust = 0.5, size=14),
			axis.text=element_text(color='#707070', size=12),
			legend.text = element_text(color='#333333', size=12),
			legend.title = element_text(color='#333333', size=12),
			axis.title=element_text(size=16),
			axis.ticks = element_line(color = "#707070"),
			panel.grid.major = element_line('#eaeaea', 0.5),
		  panel.grid.minor = element_line('#eaeaea', 0.5),
		  legend.background = element_rect(fill='transparent'),
		  strip.background = element_blank(), 
		  strip.text = element_text(size=14))


format_digit = function(x, digits=2){
	return(round(x, digits))
}


format_col = function(x1, x2){
		return(paste(paste(paste(x1, " (", sep=''),
			round((x1 / x2)*100, 2), sep=''),
			'%)', sep=''))
}
# now make plots
# read in color file as argument
# colors
colors = c(
	'insti'="#CD5C5C", 
	'INSTI' = "#CD5C5C", 
	'nnrti'='#4682B4', 
	'NNRTI' = '#4682B4',
	'treatment-naive nnrti'="#4682B4",
	'treatment-experienced nnrti'='#a2c0d9',
	'nrti'="#56731E",
	'NRTI' = "#56731E",
	'pi'="#B46B93",
	'PI' = "#B46B93",
	'M'="#4682B4",
	'F'="#89446B",
	"(34,50]"="#8693CA",
	"(24,34]"="#4682B4",
	"(14,24]"="#267B7F",
	"A1" = "#5C927C",
	"B" = "#7c5c92",
	"C" = "#927c5c",
	"AG" = "#1c7953",
	"BF" = "#72925c",
	"06A1" = "#5c8d92",
	"A2" = "#3e6354",
	"02G" = "#8e214d",
	"DF" = "#92615c",
	"H" = "#8D925C",
	"BF1" = "#7fa068",
	"AE" = "#518B89",
	"D" = "#925c72",
	"AD" = "#5c7292",
	"cpx" = "#615c92",
	"A1D" = "#1c4279",
	"BC" = "#79531c",
	"CD" = "#1c7179",
	"G" = "#791c42",
	"A1BD" = "#791c42",
	"A1CD" = "#615c92",
	"A1BCD" = "#7c5c92",
	"A1C" = "#f7f0c7",
	"A1DF1" = "#84c8ed",
	"01A1D" = "#1c2b79",
	"other" = "#eaeaea",
	"other/missing" = "#eaeaea",
	"other (recombinant)" = "#eaeaea",
	"Agrarian" = "#8EB9A6",
	"Fishing" = "#96B3D2",
	"Trading" = "#ADA5D5",
	"agrarian" = "#8EB9A6",
	"fishing" = "#96B3D2",
	"trading" = "#ADA5D5",
	'all' = '#888888',
	'all PLHIV' = '#888888',
	'pre-treatment viremic' = '#a05851',
	'treatment-experienced viremic' = '#5172a0',
	'viremic' = '#5199a0',
	'in' = "#DBC665",
	'pr' = "#B46B93",
	'rt' = "#4682B4"
	)

# drug names
nnrti_cols = c("DOR", "EFV", "ETR", "NVP", "RPV")
nnrti_drugs = c("doravirine", "efavirenz", "Etravirine", "nevirapine", "rilpivirine")

nrti_cols = c("ABC", "AZT", "D4T", "DDI", "FTC", "3TC", "TDF")
nrti_drugs = c("abacavir", "zidovudine", "stavudine", "didanosine", "emtricitabine", "lamivudine", "tenofovir")

pi_cols = c("ATV", "DRV", "FPV", "IDV", "LPV", "NFV", "SQV", "TPV")
pi_drugs = c("atazanavir", "darunavir", "fosamprenavir", "indinavir", "loprinavir", "nelfinavir", "saquinavir", "tipranavir")

insti_cols = c("BIC", "DTG", "EVG", "RAL", "CAB")
insti_drugs = c("bictegravir", "dolutegravir", "elvitegravir", "raltegravir", "cabotegravir ")

all_cols = c(nnrti_cols, nrti_cols, pi_cols, insti_cols)

# rename stratification variables
# to do read from file
var_rename = c(
	"age_cat" = "age category",
	"subtype" = "subtype",
	"comm_type" = "community type",
	"sex" = "sex",
	'F' = 'female',
	'M' = 'male',
	'in' = 'integrase',
	'rt' = 'reverse transcriptase',
	'pr' = 'protease',
	'ref_subtype' = 'subtype\n(best reference)',
	'genome_subtype' = 'subtype\n(genome)',
	'pol_subtype' = 'subtype\n(pol)',
	'vl_cat' = "viral load",
	'round' = ' survey round'
)

patterns = c(
	'NNRTI' = 'none',
	'NRTI' = 'none',
	'PI' = 'none',
	'NNRTI, NRTI' = 'stripe',
	'NNRTI, NRTI, PI' = 'stripe',
	'NNRTI, PI' = 'stripe',
	'NRTI, PI' = 'stripe')

shapes=c(
	'insti' = 23, 'nnrti' = 24, 'nrti' = 25, 'pi' = 22, 'all' = 21,
	'NNRTI' = 24, 'NRTI' = 25, 'PI' = 22, 
	'rt' = 24, 'in'=23, 'pr'=22,
	'viremic pre-treatment' = 21,
	'viremic treatment-experienced' = 21,
	'viremic' = 21)

strata = c('sex', 'age_cat', 'comm_type', 'subtype')
#drop_stratum = c('subtype:B')


class_round_prev_rr_table = function(pred_file, rr_file){
	sort_cols = function(x){
		class_order = c("nnrti" = 1, "nrti" = 2, "pi" = 3)
		type_order = c("Prev. % (95% CI)" = 1, "Prev. ratio (95% CI)" = 2, "p-value" = 3, "spacer" = 4)
		non_round_x = x[x!= "round"]
		s1 = class_order[str_split(non_round_x, "_", simplify=T)[,1]]
		s2 = type_order[str_split(non_round_x, "_", simplify=T)[,2]]
		scores = tibble(val = non_round_x, s1 = s1, s2 = s2) %>%
			arrange(s1, s2)
		return(c("round", scores$val))
	}
	pred = read_tsv(pred_file)
	# format for table
	pred = pred %>%
		select(-se.fit) %>%
		mutate(across(fit:upr, ~format_digit(.x*100))) %>%
		unite("val", fit:lwr, sep=' (') %>%
		unite("val", val:upr, sep=', ') %>%
		mutate(
			val = paste(val, ')', sep=''),
			type = paste(class, "Prev. % (95% CI)", sep='_')) %>%
		select(-class)

	rr = read_tsv(rr_file) %>%
		mutate(across(RR:UCI, ~format_digit(.x))) %>%
		unite("val", RR:LCI, sep=' (') %>%
		unite("val", val: UCI, sep=', ') %>%
		mutate(
			val = paste(val, ')', sep=''),
			val = if_else(`var` == "(Intercept)", 'ref', val),
			P = if_else(P < 1E-4, "<0.0001", as.character(signif(P, 2))),
			P = if_else(`var` == "(Intercept)", 'ref', P),
			`var` = as.numeric(str_split(`var`, 'round', simplify=T)[,2]),
			`var` = if_else(is.na(`var`), 
				min(`var`, na.rm=TRUE)-1,
				`var`),
			type = paste(class, "Prev. ratio (95% CI)", sep='_')) %>%
		rename(c("round"='var')) %>%
		select(-class)

	rr = 
		bind_rows(
			rr %>% select(-P),
			rr %>% select(-val) %>% rename(c("val" = "P")) %>%
				mutate(type = paste(str_split(type, "_", simplify=T)[,1], "p-value", sep='_')))

	t = bind_rows(pred, rr)
	t = t %>% pivot_wider(names_from=type, values_from=val)
	# add blank columns
	for (class in unique(read_tsv(pred_file)$class)){
		t[paste(class, 'spacer', sep='_')] = NA
	}
	# sort
	t = t %>% select(sort_cols(colnames(t))) %>%
		arrange(round)
	return(t)
}


run_process_model = function(d, m_formula, out_file, corstr='fit', id_col='study_id', order_col='round'){
	suppressMessages(require(geepack))
	suppressMessages(require(tidyverse))
	output = list()
	# for GEE need to sort by individual and time
	# make order column numeric to enable sorting
	d[["order_col_numeric"]] = as.numeric(d[[order_col]])
	# add id column for geepack
	d = d %>% group_by_at(id_col) %>%
		mutate(idx = cur_group_id()) %>%
		ungroup() %>% 
		arrange(!!as.symbol(id_col), 
			order_col_numeric)
	# add check to ensure sorted for GEE
	stopifnot((d %>% mutate(delta_idx = idx - lag(idx)) %>%
		summarise(m = min(delta_idx, na.rm=TRUE)))$m >= 0)
	# only do this if more than one round
	if (length(unique(d[["order_col_numeric"]])) > 1){
		stopifnot((d %>% group_by(idx) %>% mutate(delta_order = 
				order_col_numeric - 
					lag(order_col_numeric)) %>% 
			ungroup() %>%
			select(study_id, round, delta_order) %>%
			summarise(m = min(delta_order, na.rm=TRUE)))$m > 0)
	}
	if (corstr == 'fit'){
		all_models_q = list()		
		all_models_m = list()		
		for (corr in c("independence", "exchangeable")){
			m = geeglm(m_formula, data=d, id = idx, family=poisson(link="log"), corstr=corr, 
				weights=w, waves=order_col_numeric)
			q = QIC(m)
			all_models_q[[corr]] = tibble(type = names(q), val=q, corr=corr)
			all_models_m[[corr]] = m
		}
		# ar1 is only different if >2 time points
		if (length(unique(d[[order_col]])) > 2){
			for (corr in c("ar1")){
				m = geeglm(m_formula, data=d, id = idx, family=poisson(link="log"), corstr=corr, 
					weights=w, waves=order_col_numeric)
				q = QIC(m)
				all_models_q[[corr]] = tibble(type = names(q), val=q, corr=corr)
				all_models_m[[corr]] = m
			}
		}
		all_models_q = bind_rows(all_models_q) %>% pivot_wider(names_from=type, values_from=val)
		# identify correlation structure with minimum qic
		min_qic_corr = (all_models_q %>% filter(QIC == min(QIC)))$corr[[1]]
		m = all_models_m[[min_qic_corr]]
	}
	else{
		m = geeglm(m_formula, data=d, id = idx, family=poisson(link="log"), corstr=corstr, 
				weights=w, waves=order_col_numeric)
		min_qic_corr = NA
		all_models_q = NA
	}
	# save formatted and unformatted table
	r = summary(m)$coefficients
	o = cbind(exp(cbind(
		RR = r[,1], 
	                   LCI = r[,1] + qnorm(0.05/2)*r[,2],
	                   UCI = r[,1] - qnorm(0.05/2)*r[,2])),
	                   P = r[,4])
	rownames(o) = rownames(r)
	saveRDS(m, paste(c(out_file, '.rds'), collapse=''))
	output[['o']] = o
	output[['m']] = m
	output[['corr']] = all_models_q
	output[['min_qic_corr']] = min_qic_corr
	return(output)
}


get_pg_rccs_mapping = function(rccs_dat, neuro_dat){
	pg_rccs_id = rccs_dat %>% drop_na(Pangea.id) %>% 
		select(RCCS_studyid, visit, Pangea.id) %>%
		mutate(
			Pangea.id = str_split(Pangea.id, "-", simplify=T)[,2]) 
	# add in Neuro study metadata
	pg_rccs_id = rbind(
		pg_rccs_id, 
			neuro_dat %>%
			select(studyid, Pangea.id) %>%
			mutate(
				Pangea.id = str_split(Pangea.id, "-", simplify=T)[,2],
				visit="NEU") %>%
			rename(RCCS_studyid = studyid)) %>% unique()
	# need to ensure that pangea IDs are unique
	# first check that none are repeated
	print(nrow(pg_rccs_id %>% unique() %>% group_by(Pangea.id) %>% filter(n() > 1)))
	# now check that when we have longitudinal data the samples
	# from the same person have unique pangea IDs
	print(nrow(pg_rccs_id %>% unique() %>% group_by(RCCS_studyid) %>%
		filter(length(unique(Pangea.id)) != n())))
	# turned into a named vector for easy lookup
	pg_rccs_mapping = paste(pg_rccs_id$RCCS_studyid, as.character(pg_rccs_id$visit), sep='R')
	names(pg_rccs_mapping) = pg_rccs_id$Pangea.id 
	return(pg_rccs_mapping)
}

read_format_rccs_dat = function(rccs_df){
	rccs_metadata = read_dta(rccs_df)
	# foramt round
	rccs_metadata = rccs_metadata %>% 
		mutate(
			round = as.numeric(str_replace(str_remove(round, "^R0+"), 'S', '.1')),
			int_date = 
				as.Date(int_date, 
					format = 
						ifelse(grepl("/", int_date, fixed = TRUE), 
							ifelse(nchar(str_split(int_date, "/", simplify=TRUE)[,3]) == 2, 
								"%d/%m/%y", "%d/%m/%Y"),
							"%Y-%m-%d")))
	# bad examples: 
	# H013632  21/06/1995 R001
	# G113713  06/11/2020
	# todo add date formatting!
	return(rccs_metadata)
}


#get_incident = function(rccs_metadata, tspan){
#	last_neg = rccs_metadata %>% filter(finalhiv == 'N') %>%
#		group_by(study_id) %>% 
#			summarise(last_neg = max(round), last_neg_date = max(int_date))
#	first_pos = rccs_metadata %>% filter(finalhiv == 'P') %>%
#		group_by(study_id) %>% 
#		summarise(first_pos = min(round), first_pos_date = min(int_date))
#	incident = last_neg %>% inner_join(first_pos) %>% 
#		mutate(
#			inf_date = (first_pos_date - last_neg_date)/2 + last_neg_date,
#			incident = (first_pos > last_neg) & (first_pos_date - last_neg_date) <= tspan,
#			incident_all = TRUE)
#	return(incident)
#}
		

get_incident = function(rccs_metadata){
	last_neg = rccs_metadata %>% filter(finalhiv == 'N') %>%
		group_by(study_id) %>% 
			summarise(last_neg = max(round))
	first_pos = rccs_metadata %>% filter(finalhiv == 'P') %>%
		group_by(study_id) %>% 
		summarise(first_pos = min(round))
	incident = last_neg %>% inner_join(first_pos, by = 'study_id') %>% 
		mutate(
			incident = (first_pos > last_neg) & (first_pos - last_neg) == 1)
	return(incident)
}

get_first_arv = function(rccs_metadata){
	first_arv = rccs_metadata %>% select(study_id) %>% unique() %>% 
		left_join(
			rccs_metadata %>% select(study_id, round, int_date, arvmed, cuarvmed) %>%
				mutate(any_arv = replace_na(arvmed, 0) == 1 | replace_na(cuarvmed, 0) == 1) %>%
				filter(any_arv) %>%
				group_by(study_id) %>%
				filter(round == min(round)) %>%
				select(-c(arvmed, cuarvmed, any_arv)) %>%
				ungroup() %>%
				rename(first_arv_round = round, 
					first_arv_date = int_date)) %>% 
				mutate(first_arv_round = replace_na(first_arv_round, Inf),
					first_arv_date = replace_na(first_arv_date, as.Date(Inf)))
	return(first_arv)
}


calc_wilson_ci = function(p,n, alpha=0.05){
	#https://www.statisticshowto.com/wilson-ci/
	#doi:10.1080/01621459.1927
	z = qnorm(1-alpha/2)
	q = 1-p
	# break the formulat into three components
	comp_1 = p + z^2 / (2*n)
	comp_2 = z * ((p*q)/n + z^2 / (4*n^2))^0.5
	comp_3 = 1 + z^2/n
	ul = (comp_1 + comp_2) / comp_3
	ll = (comp_1 - comp_2)/comp_3
	return(cbind(ll, ul))
}


get_inf_date = function(rccs_metadata){
	last_neg = rccs_metadata %>% 
		select(study_id, round, int_date, finalhiv) %>%
		filter(finalhiv=='N') %>%
		group_by(study_id) %>%
		filter(round == max(round)) %>%
		ungroup() %>%
		select(-finalhiv) %>%
		rename(last_neg_round = round,
			last_neg_date = int_date)
	first_pos = rccs_metadata %>% 
		select(study_id, round, int_date, finalhiv) %>%
		filter(finalhiv=='P') %>%
		group_by(study_id) %>%
		filter(round == min(round)) %>%
		ungroup() %>%
		select(-finalhiv) %>%
		rename(first_pos_round = round,
			first_pos_date = int_date)
	rccs_metadata = rccs_metadata %>% 
		left_join(last_neg, by='study_id') %>%
		left_join(first_pos, by='study_id') %>%
		mutate(
			inf_date = case_when(
				!is.na(last_neg_date) & !is.na(first_pos_date) ~ pmin(last_neg_date, first_pos_date) + abs(first_pos_date-last_neg_date)/2),
			inf_uncertainty = case_when(
				!is.na(last_neg_date) & !is.na(first_pos_date) ~ abs(first_pos_date-last_neg_date))) 
	return(rccs_metadata)
}

