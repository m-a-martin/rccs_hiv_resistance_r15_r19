suppressMessages(library(tidyverse))

# formatted data
hiv_dr = read_tsv('data/rakai_drug_resistance_formatted.tsv', show_col_types = FALSE)
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


extract_col = function(x, idx){
	if (all(is.na(x))){
		return(NA)
	}else{
		return(x[,idx])
	}
}


hiv_muts = bind_cols(
	tibble(
		mut = 
			replace_na(
				unlist(lapply(
					str_split(hiv_dr$Mut, "\\|"),
					function(x) str_split(x, ";", simplify=TRUE)[,1])),
				'WT'),
		freq = 
			replace_na(
				unlist(
					lapply(str_split(hiv_dr$Mut, "\\|"), function(x) extract_col(str_split(x, ";", simplify=TRUE), 3))),
				""),
		reads = 
			replace_na(
				unlist(
					lapply(str_split(hiv_dr$Mut, "\\|"), function(x) extract_col(str_split(x, ";", simplify=TRUE), 2))),
				"")),	
	(hiv_dr %>% 
				select(study_id, round, finalhiv, viremic, pre_treatment))[
			rep(seq(1,nrow(hiv_dr)), if_else(hiv_dr$n_muts==0, 1, hiv_dr$n_muts)),]) %>%
	filter(mut != 'WT')


write_tsv(hiv_muts, 'data/rakai_drug_resistance_mut_formatted.tsv')
