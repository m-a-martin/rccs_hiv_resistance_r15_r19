library(tidyverse)
library(readxl)
library(haven)
source('scripts/utils.R')


format_rccs_metadata = function(rccs_metadata, v=c(), rounds=c()){
	# list of covariates we want to include
	# make new variables
	if ("ageyrs" %in% v){
		breaks=c(14, 24, 34, 50, 100)
		rccs_metadata = rccs_metadata %>% 
			mutate(
				ageyrs = ifelse(ageyrs == 98, NA, ageyrs),
				age_cat = cut(ageyrs, breaks=breaks))
		v = c("age_cat", v[!(v %in% c("ageyrs"))])
	}

	if ("comm_type" %in% v){
		rccs_metadata = rccs_metadata %>% mutate(comm_type = as_factor(comm_type))
	}

	if ("copies" %in% v){
		vl_cat = rccs_metadata$copies
		vl_cat[vl_cat == ''] = 'missing'
		vl_cat[vl_cat == 'BD'] = 'BD'
		where_numeric = vl_cat != 'missing' & vl_cat != 'BD'
		vl_cat[where_numeric] = gsub(',','',vl_cat[where_numeric])
		vl_cat[where_numeric] = 
			as.character(cut(log10(as.numeric(vl_cat[where_numeric])), breaks=c(0,1,3,4,5,Inf)))
		rccs_metadata = 
			rccs_metadata %>% 
				mutate(
					vl_cat = vl_cat, 
					missing_vl = copies == '',
					bd_vl = copies == 'BD',
					log10vl = case_when(
						copies != '' & copies != 'missing' & copies != 'BD' ~ log10(as.numeric(gsub(',','',copies)))))
		v = c(v[!(v %in% c("copies"))], 'copies', "vl_cat",  'log10vl')
	}
	# add in treatment info
	# if viral load is BD we assume not pre-treatment
	rccs_metadata = rccs_metadata %>% 
		left_join(get_first_arv(rccs_metadata), by='study_id') %>%
		left_join(get_incident(rccs_metadata), by='study_id') %>%
		mutate(
			pre_treatment = round < first_arv_round,
			viremic = case_when(
				finalhiv == 'N' | finalhiv == '' ~ FALSE,
				finalhiv == 'P' & pre_treatment & copies == '' ~ TRUE,
				!pre_treatment & copies == '' ~ FALSE,
				finalhiv == 'P' & !is.na(copies) & as.numeric(gsub(',','',copies)) >= 1000 ~ TRUE,
				!is.na(copies) & as.numeric(gsub(',','',copies)) <1000 ~ FALSE,
				!is.na(copies) & copies == 'BD' ~ FALSE))
	# filter rccs_metadata based on inclusion criteria
	rccs_metadata = rccs_metadata %>% 
		select(study_id, round, finalhiv, pre_treatment, last_neg, first_pos,  incident, viremic, all_of(v)) %>%
		filter(round %in% rounds & age_cat != '(50,100]')
	# further filter missing viral loads
	missing_vl = rccs_metadata %>% filter(round > 15 & finalhiv == 'P' & vl_cat == "missing") %>% 
		select(study_id, round, finalhiv, copies, vl_cat)
	write_tsv(missing_vl, 'data/rakai_drug_resistance_missing_vl.tsv')
	print(paste(nrow(missing_vl), "hiv+ participants-visits from round 16+ are missing viral load measurements"))
	print(nrow(rccs_metadata))
	rccs_metadata = rccs_metadata %>%
		filter(round == 15 | round > 15 & finalhiv != 'P' | round > 15 & finalhiv == 'P' & vl_cat != 'missing')
	print(nrow(rccs_metadata))
	return(rccs_metadata)
}

# todo loop over sheets
hiv_dr = read_excel("data/Rakai_Drug_Resistance_20220921.xlsx", sheet="5%_10reads")  %>%
	mutate(
		sampleID = str_split(sampleID, '5pct', simplify=TRUE)[,1])
names(hiv_dr) = gsub('/r', '', names(hiv_dr))

# read in subtype data
subtype_dat = 
	read_tsv('data/PANGEA_Genomes_WithSubtype.tsv') %>% 
	select(sequence_id, subtype_bestref, Genome_rip4_sub, Pol_rip4_sub) %>%
	rename(pangea_id = sequence_id, 
		ref_subtype = subtype_bestref, 
		genome_subtype = Genome_rip4_sub, 
		pol_subtype = Pol_rip4_sub) %>%
	mutate(pangea_id = str_split(pangea_id, "_", simplify=T)[,1])

# add to resistance dat
hiv_dr = hiv_dr %>%
	mutate(
		pangea_id = str_split(sampleID, "_", simplify=T)[,1]) %>%
	left_join(subtype_dat)

# read in mapping to rakai ID and date
metadata = read_csv('data/rakai_sequence_id_mappings_cohort.csv') %>%
	mutate(drmseq_prefix = str_split(drmseq_prefix, '20pct', simplify=TRUE)[,1])

# now merge
hiv_dr = hiv_dr %>% 
	left_join(metadata %>% mutate(sampleID = drmseq_prefix),
	by=c('sampleID')) %>%
	mutate(
		study_id = str_split(pt_id, "-", simplify=TRUE)[,2]) %>%
	group_by(ref_subtype) %>%
	mutate(
		subtype = 
			case_when(
				is.na(ref_subtype) ~ 'missing',
				n() >= 300 & !is.na(ref_subtype) ~ ref_subtype,
				n() < 300 & grepl('_', ref_subtype) & !is.na(ref_subtype)~ 'other (recombinant)',
				n() < 300 & !grepl('_', ref_subtype) & !is.na(ref_subtype) ~ 'other (non-recombinant)')) %>%
	ungroup()

# filter our data from Neuro study
hiv_dr = hiv_dr %>% filter(cohort_id == 'U-JHU-RCCS')

# read in rccs metadata
rccs_metadata_backup = read_format_rccs_dat('data/RCCSdata_R001_R019_VOIs.dta')
rccs_metadata = rccs_metadata_backup

# format rccs metadata
v = c('ageyrs', 'sex', 'comm_type', 'comm_num', 'copies', 'int_date')
rounds = seq(15,19)
rccs_metadata = format_rccs_metadata(rccs_metadata, v, rounds)

# and merge with dr data
hiv_dr = rccs_metadata %>% 
	rename(visit_dt = int_date) %>%
	left_join(
		hiv_dr %>% 
			mutate(
				dr_dat = TRUE,
				valid_dr_dat = !if_all(all_of(all_cols), function(x) x %in% c("X"))), 
		by=c('study_id', 'visit_dt')) %>%
	mutate(dr_dat=replace_na(dr_dat, FALSE),
		valid_dr_dat=replace_na(valid_dr_dat, FALSE))

# add survey dates
survey_dates = rccs_metadata_backup %>% group_by(round) %>% 
	summarise(round_mid_date = quantile(int_date, c(0.5), type=1, na.rm=TRUE))

hiv_dr = hiv_dr %>% left_join(survey_dates %>% select(round, round_mid_date), by="round")

# select columns
hiv_dr = hiv_dr %>% select(
	-c('main_cohort_id', 'cohort_id', 'pt_id',
		'drmseq_prefix', 'sequence_id'))

# and save
write_tsv(hiv_dr, "data/rakai_drug_resistance_formatted.tsv")

# save anonymized data for GitHub
# turn community numbers into NA
study_ids = unique(hiv_dr$study_id)
# generate random IDs
set.seed(111) 
random_ids = sample(1:length(study_ids), length(study_ids), replace=F)
random_ids = sprintf(paste(c('%0', 
		nchar(as.character(length(study_ids))), "d"), collapse=''), 
	random_ids)
# create mapping
random_ids = tibble(study_id=study_ids, anonymized_id=random_ids)
write_tsv(random_ids, "data/anonymized/anonymized_study_id_map.tsv")

anon_hiv_dr = hiv_dr %>% left_join(random_ids, by='study_id') %>% 
	select(-c(study_id, pangea_id, comm_num, sampleID)) %>%
	rename(study_id = anonymized_id) 

# reorder columns
anon_hiv_dr = anon_hiv_dr %>%
	select(c('study_id', names(anon_hiv_dr)[seq(1,length(names(anon_hiv_dr))-1)]))

# and save
write_tsv(anon_hiv_dr, 'data/anonymized/rakai_drug_resistance_formatted.tsv')



names(hiv_dr)
