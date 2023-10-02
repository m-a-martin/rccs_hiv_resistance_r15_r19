library(tidyverse)
library(stats)
library(lmtest)
library(sandwich)

# data files
dat_file = 'data/rakai_drug_resistance_categorized.tsv'

# read in data
hiv_dr_cat = read_tsv(dat_file)

# just HIV+ pretreatment incident individuals
# only data from incident individuals
# only data from their first post-infection RCCS visit
pretreat_dat = hiv_dr_cat %>% 
	filter(
		viremic &
		incident & 
		pre_treatment & 
		finalhiv == "P" &
		round == first_pos) %>%
	mutate(nnrti_dat = !is.na(nnrti),
		missing_vl = vl_cat == 'missing',
		log10vl = replace_na(log10vl, 0))


all_weights = list()
# for each participant visit, get probability of having nnrti resistance data
drug = 'nnrti'
temp_dat = pretreat_dat %>% 
	mutate(round = as.character(round)) %>%
	rename(d = paste(drug, "_dat", sep='')) %>%
	filter(round == 15 | (round != 15 & !missing_vl))

# include all interaction terms so this is essentially a per-round model
m = glm(d ~ round*missing_vl + round*log10vl + round*age_cat + round*sex, 
	data = temp_dat, 
	family=poisson(link="log"))

# don't actually use the CI for anything 
r = coeftest(m, vcov=sandwich)
# create output files
o <- cbind(exp(cbind(RR = r[,1], 
	       LCI = r[,1] + qnorm(0.05/2)*r[,2],
	       UCI = r[,1] - qnorm(0.05/2)*r[,2])),
	       P = r[,4])
write.table(o, 
	paste(c('models/pretreat_inc_', drug, '_weights_raw.tsv'), collapse=''), 
	sep='\t', col.names=NA)
saveRDS(m, paste(c('models/pretreat_inc_', drug, '_weights.rds'), collapse=''))
# now predict probability of sequencing success
# among those which were successfully sequenced
# note: rank deficiency warning is because there is 
# only missing viral loads in round 15
# not an issue because prediction data is same 
# as modeled data
p = predict(m, temp_dat, type='response')

# add p to temp data
# for each calculate a probability of NOT having sequence data
temp_weight_dat = 
	temp_dat %>% select(study_id, round, d) %>% 
		mutate(drug = drug, p = p, p_not = 1-p) 

# group by study id and get 
# probability of having at least 1 sequenced sample
# P(at least 1) = 1 - P(none)
# weight is inverse probability of p
grouped_temp_weight_dat = temp_weight_dat %>% group_by(study_id, drug, d) %>%
	summarise(p = 1 - prod(p_not), w = 1/p, .groups="drop") %>%
	filter(d) %>%
	select(-d) %>%
	# normalize weights 
	mutate(w = length(unique(temp_weight_dat$study_id)) * w / sum(w))

all_weights[[drug]] = grouped_temp_weight_dat

# format weights
all_weights = bind_rows(all_weights)
all_weights = all_weights %>% 
	select(c(study_id, drug, p)) %>% 
	pivot_wider(names_from=drug, values_from=p, names_prefix='p_') %>%
	left_join(
		all_weights %>% 
			select(c(study_id, drug, w)) %>% 
			pivot_wider(names_from=drug, values_from=w, names_prefix='w_'),
		by=c('study_id'))


# and finally save to file
write_tsv(all_weights, "models/pretreat_inc_weights.tsv")


