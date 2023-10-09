library(tidyverse)
library(stats)
library(lmtest)
library(sandwich)


# read in categorized dr data
hiv_dr_cat = read_tsv('data/rakai_drug_resistance_categorized.tsv') 

#### SAMPLING MODEL ####
# per round and per drug
# outcome: drug resistance prediction
# predictors: missing_vl, log10vl, age_category, sex, community type
# population: pretreatment HIV+ individuals

# subset to just pre-treatment, HIV-infected individuals
# remove weird ages and early round data
treat_dat = hiv_dr_cat %>% filter(!pre_treatment & finalhiv == 'P' & viremic) %>%
	mutate(
		nnrti_dat = !is.na(nnrti), 
		nrti_dat = !is.na(nrti),
		pi_dat = !is.na(pi),
		missing_vl = vl_cat == 'missing',
		log10vl = replace_na(log10vl, 0)) %>%
	filter(round > 16 & round < 19)


# for each round for each drunk get sampling model
all_weights = list()
# first just general sampling model 
#temp_dat = treat_dat %>% 
#			mutate(round = as.character(round)) %>%
#			mutate(d = valid_dr_dat)
# include all interaction terms so this is essentially a per-round model
#m = glm(d ~ round*missing_vl + round*log10vl + round*comm_type + round*age_cat + round*sex, 
#	data = temp_dat, 
#	family=poisson(link="log"))
# don't actually use the CI for anything but good to see our predictors are significant
#r = coeftest(m, vcov=sandwich)
# create output files
#o <- cbind(exp(cbind(RR = r[,1], 
#	       LCI = r[,1] + qnorm(0.05/2)*r[,2],
#	       UCI = r[,1] - qnorm(0.05/2)*r[,2])),
#	       P = r[,4])
#write.table(o, 
#	paste(c('models/treat_gen_weights_raw.tsv'), collapse=''), 
#	sep='\t', col.names=NA)
#saveRDS(m, paste(c('models/treat_gen_weights.rds'), collapse=''))
# now predict probability of sequencing success
# among those which were successfully sequenced
#p = predict(m, temp_dat, type='response')
#temp_weight_dat = 
#	temp_dat %>% select(study_id, round, d) %>% 
#		mutate(drug = 'gen', p = p) %>% 
#		group_by(round) %>%
#		mutate(x = n()) %>%
#		filter(d == 1) %>%
#		mutate(w = x * p / sum(p))
#all_weights[['gen']] = temp_weight_dat

for (drug in c('nnrti', 'nrti', 'pi')){
		# temporary dataset renaming outcome variable
		# can't use dynamic names in glm formula (as far as I can tell)
		temp_dat = treat_dat %>% 
			mutate(round = as.character(round)) %>%
			rename(d = paste(drug, "_dat", sep=''))
		# include all interaction terms so this is essentially a per-round model
		m = glm(d ~ round*missing_vl + round*log10vl + round*comm_type + round*age_cat + round*sex, 
			data = temp_dat, 
			family=poisson(link="log"))
		# don't actually use the CI for anything but good to see our predictors are significant
		r = coeftest(m, vcov=sandwich)
		# create output files
		o <- cbind(exp(cbind(RR = r[,1], 
			       LCI = r[,1] + qnorm(0.05/2)*r[,2],
			       UCI = r[,1] - qnorm(0.05/2)*r[,2])),
			       P = r[,4])
		write.table(o, 
			paste(c('models/treat_', drug, '_weights_raw.tsv'), collapse=''), 
			sep='\t', col.names=NA)
		saveRDS(m, paste(c('models/treat_', drug, '_weights.rds'), collapse=''))
		# now predict probability of sequencing success
		# among those which were successfully sequenced
		p = predict(m, temp_dat, type='response')
		temp_weight_dat = 
			temp_dat %>% select(study_id, round, d) %>% 
				mutate(drug = drug, p = p) %>% 
				group_by(round) %>%
				mutate(x = n()) %>%
				filter(d == 1) %>%
				mutate(w = x * p / sum(p))
		all_weights[[drug]] = temp_weight_dat
}

# weights for sequences with data for all drugs
temp_dat = treat_dat %>% 
			mutate(round = as.character(round)) %>%
			mutate(d = !is.na(nnrti) & !is.na(pi) & !is.na(insti) & !is.na(nrti))
# include all interaction terms so this is essentially a per-round model
m = glm(d ~ round*missing_vl + round*log10vl + round*comm_type + round*age_cat + round*sex, 
	data = temp_dat, 
	family=poisson(link="log"))
# don't actually use the CI for anything but good to see our predictors are significant
r = coeftest(m, vcov=sandwich)
# create output files
o <- cbind(exp(cbind(RR = r[,1], 
	       LCI = r[,1] + qnorm(0.05/2)*r[,2],
	       UCI = r[,1] - qnorm(0.05/2)*r[,2])),
	       P = r[,4])
write.table(o, 
	paste(c('models/treat_all_weights_raw.tsv'), collapse=''), 
	sep='\t', col.names=NA)
saveRDS(m, paste(c('models/treat_all_weights.rds'), collapse=''))
# now predict probability of sequencing success
# among those which were successfully sequenced
p = predict(m, temp_dat, type='response')
temp_weight_dat = 
	temp_dat %>% select(study_id, round, d) %>% 
		mutate(drug = 'all', p = p) %>% 
		group_by(round) %>%
		mutate(x = n()) %>%
		filter(d == 1) %>%
		mutate(w = x * p / sum(p))
all_weights[['all']] = temp_weight_dat



# format weights
all_weights = bind_rows(all_weights)
all_weights = all_weights %>% 
	select(c(study_id, round, drug, p)) %>% 
	pivot_wider(names_from=drug, values_from=p, names_prefix='p_') %>%
	left_join(all_weights %>% 
		select(c(study_id, round, drug, w)) %>% 
		pivot_wider(names_from=drug, values_from=w, names_prefix='w_'),
		by=c('study_id', 'round'))

# and finally save to file
write_tsv(all_weights, "models/treat_weights.tsv")



