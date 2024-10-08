suppressMessages(library(tidyverse))
suppressMessages(library(stats))
suppressMessages(require(lmtest))
suppressMessages(library(sandwich))

# read in categorized dr data
hiv_dr_cat = read_tsv('data/rakai_drug_resistance_categorized.tsv', show_col_types=FALSE) 


#### SAMPLING MODEL ####
# per round and per drug class
# outcome: drug resistance prediction
# predictors: missing_vl, log10vl, age_category, sex, community type
# population: viremic PLHIV

# subset to just virmeic HIV-infected individuals
# exclude treatment-experienced PLHIV from R15 and R19
dat = hiv_dr_cat %>% 
	filter(
		finalhiv == 'P' & 
		viremic &
		(pre_treatment == TRUE | 
			(pre_treatment == FALSE & 
				(round == 17 | round == 18)))) %>%
	mutate(
		round = as.character(round),
		nnrti_dat = !is.na(nnrti), 
		nrti_dat = !is.na(nrti),
		pi_dat = !is.na(pi),
		nnrti_nrti_pi_dat = !is.na(nnrti) & !is.na(nrti) & !is.na(pi),
		insti_nnrti_nrti_pi_dat = !is.na(insti) & !is.na(nnrti) & !is.na(nrti) & !is.na(pi),
		missing_vl = vl_cat == 'missing',
		log10vl = replace_na(log10vl, 0)) %>%
	filter(age_cat != '(50, 100]' & round >= 15)


# for each drug class get sampling model
if (!dir.exists('models')) {dir.create('models')}
if (!dir.exists('models/seq_probs')) {dir.create('models/seq_probs')}

all_probs = list()
for (round in unique(dat$round)){
	round_dat = dat[dat$round == round,]
	for (drug in c('nnrti', 'nrti', 'pi', 'nnrti_nrti_pi', 'insti_nnrti_nrti_pi')){
		if (round %in% c("17", "18")){
		f = as.formula(paste(drug, '_dat ~ pre_treatment + missing_vl + log10vl + 
					comm_type + age_cat + sex', sep=''))
		}else{
			# only pre-treatment in R15,R16, and R19
			f = as.formula(paste(drug, '_dat ~ missing_vl + log10vl + 
						comm_type + age_cat + sex', sep=''))
		}
		m = glm(f,
			data=round_dat,
			family=poisson(link='log'))
		r = coeftest(m, vcov=sandwich)
		# create output files
		o <- cbind(exp(cbind(RR = r[,1], 
			       LCI = r[,1] + qnorm(0.05/2)*r[,2],
			       UCI = r[,1] - qnorm(0.05/2)*r[,2])),
			       P = r[,4])
		write.table(o, 
			paste(c('models/seq_probs/all_R', round, '_',  drug, '_probs.tsv'), collapse=''), 
			sep='\t', col.names=NA)
		saveRDS(m, paste(c('models/seq_probs/all_R', round, '_', drug, '_probs.rds'), collapse=''))
		# now predict probability of sequencing success
		# among those which were successfully sequenced
		p = predict(m, round_dat, type='response')
		temp_prob_dat = 
			round_dat[,c('study_id', 'round', paste(drug, '_dat', sep=''))] %>%
				mutate(drug = drug, p = p)
		all_probs[[paste(round, drug, sep='_')]] = 
			temp_prob_dat[temp_prob_dat[paste(drug, '_dat', sep='')] == TRUE,]
	}
}

# format probs
all_probs = bind_rows(all_probs)
all_probs = all_probs %>% 
	select(c(study_id, round, drug, p)) %>% 
	pivot_wider(names_from=drug, values_from=p, names_prefix='p_') 

# and finally save to file
write_tsv(all_probs, "models/all_sequencing_probs.tsv")


