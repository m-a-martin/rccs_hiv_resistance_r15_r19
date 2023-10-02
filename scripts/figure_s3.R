library(tidyverse)
source('scripts/utils.R')

hiv_dr_cat = read_tsv('data/rakai_drug_resistance_categorized.tsv')
weights = read_tsv('models/all_weights.tsv')

weights = weights %>% select(study_id, round, 'w_nnrti', 'w_nrti', 'w_pi') %>%
	pivot_longer(-c('study_id', 'round'), names_to='class', values_to='w') %>%
	filter(round == 18 & !is.na(w)) %>%
	mutate(class = str_split(class, "_", simplify=T)[,2])

res_prev = hiv_dr_cat %>% 
	filter(round == 18 & finalhiv == 'P' & viremic) %>% 
	select(study_id, round, finalhiv, viremic, pre_treatment, nnrti, nrti, pi) %>%
	pivot_longer(-c('study_id', 'round', 'finalhiv', 'viremic', 'pre_treatment'), names_to='class', values_to='resistance') %>%
	filter(!is.na(resistance)) %>%
	left_join(weights, by=c('study_id', 'round', 'class')) 

res_prev = res_prev %>% mutate(cat = case_when(
		resistance == 'susceptible' ~ 'susceptible',
		resistance == 'intermediate/high' & pre_treatment ~ 'pre-treatment resistant',
		resistance == 'intermediate/high' & !pre_treatment ~ 'treatment-experienced resistant')) %>%
	group_by(class, cat) %>%
	summarise(n=sum(w)) %>%
	group_by(class) %>%
	mutate(p=n/sum(n),
	cat=ordered(cat, 
		levels=c('susceptible', 'pre-treatment resistant', 'treatment-experienced resistant')))


p = ggplot(res_prev, aes(x=class, y=p*100, fill=cat)) +
	geom_bar(stat='identity', position='stack', color='#333333') +
	scale_fill_manual(values=c("#eaeaea",  "#89446B", "#4682B4"), name='') +
	ylab('prevalence among viremic PLWHIV (%)') +
	xlab('drug class') +
	gtheme +
	theme(legend.text = element_text(size=10),
		axis.text = element_text(size=16))

ggsave('figures/final/figure_s3.pdf', p, width=6.4, height=4.8)