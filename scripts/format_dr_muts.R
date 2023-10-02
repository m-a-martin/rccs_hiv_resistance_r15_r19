library(tidyverse)
#source('scripts/utils.R')

# formatted data
hiv_dr = read_tsv('data/rakai_drug_resistance_formatted.tsv')
drug = c(
		"DOR", "EFV", "ETR", "NVP", "RPV",
		"ABC", "AZT", "D4T", "DDI", "FTC", "3TC", "TDF",
		"ATV", "DRV", "FPV", "IDV", "LPV", "NFV", "SQV", "TPV", 
		"BIC", "DTG", "EVG", "RAL", "CAB")
drug_cols = c(drug, "Mut", "Haplo")
group_cols = names(hiv_dr)[!names(hiv_dr) %in% drug_cols]

# add number of mutations
hiv_dr = hiv_dr %>% mutate(n_muts = replace_na(str_count(Mut, "\\|")+1, 0)) %>%
	filter(age_cat != '(50, 100]' & round > 14)

extract = function(x, idx){
	if (all(is.na(x))){
		return(NA)
	}else{
		return(x[idx])
	}
}

hiv_muts = bind_cols(
	tibble(mut = unlist(lapply(as.vector(t(str_split(hiv_dr$Mut, "\\|", simplify=T))),
			function(x) str_split(x, ";", simplify=T)[,1])),
		freq=unlist(lapply(as.vector(t(str_split(hiv_dr$Mut, "\\|", simplify=T))),
			function(x) rev(str_split(x, ";", simplify=T))[1])),
		reads=unlist(lapply(as.vector(t(str_split(hiv_dr$Mut, "\\|", simplify=T))),
			function(x) extract(str_split(x, ";", simplify=T), 2)))) %>%
		filter(!is.na(mut) & mut != ''),
	(hiv_dr %>% 
			select(study_id, round, finalhiv, viremic, pre_treatment))[
		rep(seq(1,nrow(hiv_dr)), hiv_dr$n_muts),]) %>% unique()


write_tsv(hiv_muts, 'data/rakai_drug_resistance_mut_formatted.tsv')
