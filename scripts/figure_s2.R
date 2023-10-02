library(tidyverse)
library(cowplot)
source('scripts/utils.R')


# todo read in as arguments
dat_file = 'data/rakai_drug_resistance_categorized.tsv'

# read in data
hiv_dr_cat = read_tsv(dat_file) 

pretreat_dat = hiv_dr_cat %>% 
	filter(pre_treatment & finalhiv == 'P' & viremic) %>%
	mutate(nnrti_dat = !is.na(nnrti), nrti_dat = !is.na(nrti), pi_dat = !is.na(pi), insti_dat = !is.na(insti)) %>%
	group_by(study_id) %>%
	summarise(
		NNRTI = sum(nnrti_dat),
		NRTI = sum(nrti_dat),
		PI = sum(pi_dat)) %>%
	pivot_longer(-c(study_id)) %>%
	group_by(name, value) %>%
	summarise(n=n()) %>%
	filter(value > 0)

treat_dat = hiv_dr_cat %>% 
	filter(!pre_treatment & finalhiv == 'P' & viremic & round > 16 & round < 19) %>%
	mutate(nnrti_dat = !is.na(nnrti), nrti_dat = !is.na(nrti), pi_dat = !is.na(pi), insti_dat = !is.na(insti)) %>%
	group_by(study_id) %>%
	summarise(
		NNRTI = sum(nnrti_dat),
		NRTI = sum(nrti_dat),
		PI = sum(pi_dat)) %>%
	pivot_longer(-c(study_id))%>%
	group_by(name, value) %>%
	summarise(n=n()) %>%
	filter(value > 0)

all_dat =  hiv_dr_cat %>% 
	filter(finalhiv == 'P' & viremic & round > 16 & round < 19) %>%
	mutate(nnrti_dat = !is.na(nnrti), nrti_dat = !is.na(nrti), pi_dat = !is.na(pi), insti_dat = !is.na(insti)) %>%
	group_by(study_id) %>%
	summarise(
		NNRTI = sum(nnrti_dat),
		NRTI = sum(nrti_dat),
		PI = sum(pi_dat)) %>%
	pivot_longer(-c(study_id))%>%
	group_by(name, value) %>%
	summarise(n=n()) %>%
	filter(value > 0)


a = ggplot(pretreat_dat, aes(x=value, y=n)) + 
	geom_bar(stat="identity", fill='#eaeaea', color='#333333') +
	facet_grid(~name) +
	xlab('\n') +
	ylab('number of participants') +
	ggtitle('viremic pre-treatment') +
		xlim(0.5, 4.5) +
	gtheme +
	theme(panel.grid.major.x = element_blank(),
		  panel.grid.minor.x = element_blank(),
		  plot.title = element_text(hjust = 0.5, size=16))


b = ggplot(treat_dat, aes(x=value, y=n)) + 
	geom_bar(stat="identity", fill='#eaeaea', color='#333333') +
	facet_grid(~name) +
	xlab('\n') +
	ylab('number of participants') +
		xlim(0.5, 4.5) +
	ggtitle('viremic treatment-experienced') +
	gtheme +
	theme(panel.grid.major.x = element_blank(),
		  panel.grid.minor.x = element_blank(),
		  plot.title = element_text(hjust = 0.5, size=16))


c = ggplot(all_dat, aes(x=value, y=n)) + 
	geom_bar(stat="identity", fill='#eaeaea', color='#333333') +
	facet_grid(~name) +
	xlab('RCCS visits with resistance prediction') +
	ylab('number of participants') +
		xlim(0.5, 4.5) +
	ggtitle('all viremic') +
	gtheme +
	theme(panel.grid.major.x = element_blank(),
		  panel.grid.minor.x = element_blank(),
		  plot.title = element_text(hjust = 0.5, size=16))

p = plot_grid(a,b,c, nrow=3)
ggsave('figures/final/figure_s2.pdf', p, width=6.4*2, height=4.8*2)

