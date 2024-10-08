suppressMessages(require(tidyverse))
suppressMessages(require(haven))
suppressMessages(require(readxl))
source('scripts/utils.R')


format_rccs_metadata = function(rccs_metadata, v=c(), rounds=c()){
	# list of covariates we want to include
	# make new variables
	# categorize age
	if ("ageyrs" %in% v){
		breaks = c(14,24,34,49)
		rccs_metadata = rccs_metadata %>% 
			mutate(
				age_cat = cut(ageyrs, breaks=breaks))
		v = c("age_cat", v)
	}
	# converte community type to factor instead of labelled factor
	if ("comm_type" %in% v){
		rccs_metadata = rccs_metadata %>% mutate(comm_type = as_factor(comm_type))
	}
	# categorize viral load 
	# encodes bd as NA
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
	rccs_metadata = rccs_metadata %>% 
		left_join(get_first_arv(rccs_metadata), by='study_id') %>%
		mutate(
			pre_treatment = round < first_arv_round,
			viremic_raw = case_when(
				finalhiv == 'N' | finalhiv == '' ~ FALSE,
				finalhiv == 'P' &  copies == '' ~ NA,
				finalhiv == 'P' & !is.na(copies) & as.numeric(gsub(',','',copies)) >= 1000 ~ TRUE,
				finalhiv == 'P' & !is.na(copies) & as.numeric(gsub(',','',copies)) < 1000 ~ FALSE,
				finalhiv == 'P' & !is.na(copies) & copies == 'BD' ~ FALSE))
	# filter rccs_metadata based on inclusion criteria
	rccs_metadata = rccs_metadata %>% 
		select(study_id, round, finalhiv, pre_treatment, viremic_raw, all_of(v)) %>%
		filter(round %in% rounds & !is.na(age_cat))
	# further filter missing viral loads
	missing_vl = rccs_metadata %>% filter(round > 15 & finalhiv == 'P' & vl_cat == "missing") %>% 
		select(study_id, round, finalhiv, copies, vl_cat)
	write_tsv(missing_vl, 'data/rakai_drug_resistance_missing_vl.tsv')
	print(paste(nrow(missing_vl), "hiv+ participants-visits from round 16+ are missing viral load measurements"))
	rccs_metadata = rccs_metadata %>%
		filter(round == 15 | (round > 15 & finalhiv != 'P') | (round > 15 & finalhiv == 'P' & vl_cat != 'missing'))
	return(rccs_metadata)
}

#### MAIN CODE STARTS HERE ####
# read in drug resistance data
#hiv_dr = read_excel("data/Rakai_Drug_Resistance_20220921.xlsx", sheet="5%_10reads")  %>%
#	mutate(
#		sampleID = str_split(sampleID, '5pct', simplify=TRUE)[,1])
hiv_dr = read_tsv("data/Rakai_Drug_Resistance_20220921_filter.tsv", show_col_types=FALSE) %>%
	filter(sheet == "5%_10reads") %>%
	mutate(
		sampleID = str_split(sampleID, '5pct', simplify=TRUE)[,1])

names(hiv_dr) = gsub('/r', '', names(hiv_dr))

# read in subtype data from pangea metadata
subtype_dat = 
	read_tsv('data/PANGEA_Genomes_WithSubtype.tsv', show_col_types=FALSE) %>% 
	select(sequence_id, subtype_bestref, Genome_rip4_sub, Pol_rip4_sub) %>%
	rename(pangea_id = sequence_id, 
		ref_subtype = subtype_bestref, 
		genome_subtype = Genome_rip4_sub, 
		pol_subtype = Pol_rip4_sub) %>%
	mutate(pangea_id = str_split(pangea_id, "_", simplify=T)[,1])

# merge subtype and resistance data
hiv_dr = hiv_dr %>%
	mutate(
		pangea_id = str_split(sampleID, "_", simplify=T)[,1]) %>%
	left_join(subtype_dat)

# read in mapping to rakai ID and date
metadata = read_csv('data/rakai_sequence_id_mappings_cohort.csv', show_col_types=FALSE) %>%
	mutate(drmseq_prefix = str_split(drmseq_prefix, '20pct', simplify=TRUE)[,1])

# now merge
hiv_dr = hiv_dr %>% 
	left_join(metadata %>% mutate(sampleID = drmseq_prefix),
	by=c('sampleID')) %>%
	mutate(
		study_id = str_split(pt_id, "-", simplify=TRUE)[,2]) %>%
	group_by(ref_subtype) %>%
	# consolidate subtype labels
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

# format rccs metadata
v = c('ageyrs', 'sex', 'comm_type', 'copies', 'int_date')
rounds = seq(15,19)
# read in rccs metadata
rccs_metadata_backup = read_format_rccs_dat('data/RCCSdata_R001_R019_VOIs.dta')
rccs_metadata = rccs_metadata_backup
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

# impute round 15 viral loads
set.seed(1)
hiv_dr = hiv_dr %>%
	left_join(
		hiv_dr %>% filter(round == 15 & finalhiv == 'P' & vl_cat != 'missing' & pre_treatment == TRUE) %>%
			group_by(round, finalhiv, pre_treatment, valid_dr_dat) %>%
			summarise(p_viremic = sum(viremic_raw)/n(), .groups='drop') %>%
			mutate(vl_cat = 'missing'),
		by=c('round', 'finalhiv', 'pre_treatment', 'valid_dr_dat', 'vl_cat')) %>% 
	mutate(viremic = case_when(
		!is.na(viremic_raw) ~ viremic_raw,
		round== 15 & pre_treatment == TRUE & is.na(viremic_raw) ~ as.logical(rbinom(n(), 1, replace_na(p_viremic, 0))),
		round == 15 & pre_treatment == FALSE & is.na(viremic_raw) ~ NA))

# add survey dates
hiv_dr = hiv_dr %>%
	left_join(
		rccs_metadata_backup %>% group_by(round) %>% 
			summarise(round_mid_date = quantile(int_date, c(0.5), type=1, na.rm=TRUE)))

# add anonymized ID
hiv_dr = hiv_dr %>% mutate(rccs_study_id = study_id) %>%
	group_by(study_id) %>%
	mutate(study_id = as.character(cur_group_id())) %>%
	ungroup() %>%
	mutate(
		study_id = str_pad(study_id, max(nchar(study_id)), pad='0', side='left'))


# finally add in sequencing qc data
seq_dat = read_csv('data/2024-10-02_pangea2_mike_sequence_id_metadata.csv', show_col_types=FALSE) %>%
	select(pt_id, visit_dt, sequence_id, readnum_hiv, mapped_num, duprate, insertsize_05, 
		insertsize_median, insertsize_95, length_strict) %>%
	unique()

hiv_dr = hiv_dr %>% left_join(seq_dat %>% select(-pt_id, -visit_dt), by=c('sequence_id'))

# select columns
hiv_dr = hiv_dr %>% select(
	-c('main_cohort_id', 'cohort_id', 'pt_id',
		'drmseq_prefix'))

# and save
write_tsv(hiv_dr %>% select(-rccs_study_id, -pangea_id, -visit_dt), "data/rakai_drug_resistance_formatted.tsv")
write_tsv(hiv_dr, "data/rakai_drug_resistance_formatted_orig_ids.tsv")



#write_csv(
#	rccs_metadata %>% 
#		filter(finalhiv == 'P' & round >= 15 & round != 15.1) %>%
#		# first time reporting art
#		left_join(
#			rccs_metadata %>% filter(arvmed == 1 | cuarvmed == 1) %>%
#				group_by(study_id) %>%
#				summarise(first_arv_round = min(round)) %>%
#				select(study_id, first_arv_round),
#			by=c('study_id')) %>%
#		mutate(first_arv_round = replace_na(first_arv_round, Inf)) %>%
#		filter(round < first_arv_round) %>%
#		# if not "no" to long term medication
#		filter((ltemmed != 2 | is.na(ltemmed))) %>%
#		filter((
#			cuarvmed != 2 | is.na(cuarvmed)) & 
#			arvmed != 2 | is.na(arvmed)) %>%
#		select(study_id, round, finalhiv, ltemmed, cuarvmed, arvmed),
#	'data/missing_art.csv')