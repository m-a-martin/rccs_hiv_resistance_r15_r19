library(tidyverse)
source('scripts/utils.R')
library(cowplot)

t97a_file = 'models/pretreat_int97a_amongPLWHIV_prev_pred.tsv'
t97a = read_tsv(t97a_file)

dat_file = 'data/rakai_drug_resistance_categorized.tsv'
mut_file = 'data/rakai_drug_resistance_mut_formatted.tsv'

hiv_mut = read_tsv(mut_file)
hiv_dr_cat = read_tsv(dat_file)


# get samples with inT97A
# and merge with all data
hiv_dr_cat = hiv_dr_cat %>% left_join(hiv_mut %>% group_by(study_id, round) %>% summarise(t97a = 'inT97A' %in% mut) %>%
	filter(t97a)) %>%
	mutate(t97a = replace_na(t97a, FALSE))

# now filter
hiv_dr_cat = hiv_dr_cat %>%
	filter(
		valid_dr_dat &
		pre_treatment &
		finalhiv == 'P' & 
		viremic)


xmin = floor(log10(min(as.numeric((hiv_dr_cat %>% filter(!t97a & !is.na(copies)))$copies))))
xmax = ceiling(log10(max(as.numeric((hiv_dr_cat %>% filter(!t97a & !is.na(copies)))$copies))))


a = ggplot(hiv_dr_cat %>% filter(!t97a & !is.na(copies) & copies != '') %>% 
		mutate(label='ref samples'), aes(x=log10vl, fill=label)) +
	geom_histogram(breaks = seq(xmin,xmax,0.2), color='#333333') +
	scale_fill_manual(values=c('#eaeaea'), name='') +
	xlab('') +
	ylab('number of samples') +
	xlim(xmin, xmax) +
	gtheme +
	theme(legend.position=c(1,1), legend.justification=c(1,1))


#### PANEL B ####
b = ggplot(hiv_dr_cat %>% filter(t97a & !is.na(copies) & copies != '') %>% mutate(label='samples with inT97A'), aes(x=log10vl, fill=label)) +
	geom_histogram(breaks = seq(xmin,xmax,0.2), color='#333333') +
	scale_fill_manual(values=c("#8fbc8f"), name='') +
  	xlab(expression(log[10]~copies/mL)) + 
	ylab('number of samples') +
	xlim(xmin, xmax) +
	gtheme +
	theme(legend.position=c(1,1), legend.justification=c(1,1))


#### PANEL C ####
mut_file = 'data/rakai_drug_resistance_mut_formatted.tsv'
mut_dat = read_tsv(mut_file)


c = ggplot(mut_dat %>% filter(mut != 'inT97A') %>% mutate(label = 'non-inT97A resistance mutations'), aes(x=freq, fill=label)) +
	geom_histogram(breaks=seq(0, 1, 0.05), color='#333333') +
	scale_fill_manual(values=c('#eaeaea'), name='') +
	ylab('count') +
	xlab('') +
	gtheme+
	theme(legend.position=c(0,1), legend.justification=c(0,1))

d = ggplot(mut_dat %>% filter(mut == 'inT97A') %>% mutate(label = 'inT97A'), aes(x=freq, fill=label)) +
	geom_histogram(breaks=seq(0, 1, 0.05), color='#333333') +
	scale_fill_manual(values=c(as.character(colors['in'])), name='') + 
	xlab('mutation frequency') +
	ylab('count') +
	gtheme +
	theme(legend.position=c(0,1), legend.justification=c(0,1))

### panel e
max_x = 100*ceiling(max(as.numeric((mut_dat %>% mutate(var_reads = str_split(reads, '/', simplify=T)[,1]))$var_reads))/100)

e = ggplot(mut_dat %>% filter(mut != 'inT97A') %>% 
		mutate(
			label = 'non-inT97A resistance mutations', 
			var_reads = as.numeric(str_split(reads, '/', simplify=T)[,1])), 
		aes(x=var_reads, fill=label)) +
	geom_histogram(breaks=seq(0, max_x, 1000), color='#333333') +
	scale_fill_manual(values=c('#eaeaea'), name='', guide='none') +
	ylab('count') +
	xlab('') +
	scale_y_log10() + 
	gtheme+
	theme(legend.position=c(0,1), legend.justification=c(0,1))

f = ggplot(mut_dat %>% 
		filter(mut == 'inT97A') %>% 
		mutate(
			label = 'inT97A',
			var_reads = as.numeric(str_split(reads, '/', simplify=T)[,1])), 
		aes(x=var_reads, fill=label)) +
	geom_histogram(breaks=seq(0, max_x, 1000), color='#333333') +
	scale_fill_manual(values=c(as.character(colors['in'])), name='', guide='none') + 
	xlab('number of supporting reads') +
	ylab('count') +
	scale_y_log10() + 
	gtheme +
	theme(legend.position=c(0,1), legend.justification=c(0,1))


p2 = plot_grid(a, c, e, b, d, f, labels=c('a', 'c', 'e', 'b', 'd', 'f'), nrow=2)
ggsave('figures/final/figure_s8.pdf', p2, width=15, height=8)


