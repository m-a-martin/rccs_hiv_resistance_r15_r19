library(tidyverse)
source('scripts/utils.R')
library(cowplot)
library(ggpattern)



dat_file = 'data/rakai_drug_resistance_categorized.tsv'
prev_dat_file = 'models/all_multi_amongPLWHIV_prev_pred.tsv'

hiv_dr_cat = read_tsv(dat_file)

nrow(hiv_dr_cat %>% 
	filter(pre_treatment & !is.na(nnrti) & !is.na(nrti) & !is.na(pi)))

hiv_dr_cat %>% 
	filter(pre_treatment & !is.na(nnrti) & !is.na(nrti) & !is.na(pi)) %>%
	mutate(
		nnrti_r = nnrti == 'intermediate/high',
		nrti_r = nrti == 'intermediate/high',
		pi_r = pi == 'intermediate/high',
		r = case_when(
			nnrti_r & nrti_r & pi_r ~ 'NNRTI, NRTI, PI',
			nnrti_r & nrti_r & !pi_r ~ 'NNRTI, NRTI',
			nnrti_r & !nrti_r & pi_r ~ 'NNRTI, PI',
			nnrti_r & !nrti_r & !pi_r ~ 'NNRTI',
			!nnrti_r & nrti_r & pi_r ~ 'NRTI, PI',
			!nnrti_r & !nrti_r & pi_r ~ 'PI',
			!nnrti_r & nrti_r & !pi_r ~ 'NRTI')) %>% print(width=Inf)


pretreat_dat =  hiv_dr_cat %>% 
	filter(pre_treatment & !is.na(nnrti) & !is.na(nrti) & !is.na(pi)) %>%
	mutate(
		nnrti_r = nnrti == 'intermediate/high',
		nrti_r = nrti == 'intermediate/high',
		pi_r = pi == 'intermediate/high',
		r = case_when(
			nnrti_r & nrti_r & pi_r ~ 'NNRTI, NRTI, PI',
			nnrti_r & nrti_r & !pi_r ~ 'NNRTI, NRTI',
			nnrti_r & !nrti_r & pi_r ~ 'NNRTI, PI',
			nnrti_r & !nrti_r & !pi_r ~ 'NNRTI',
			!nnrti_r & nrti_r & pi_r ~ 'NRTI, PI',
			!nnrti_r & !nrti_r & pi_r ~ 'PI',
			!nnrti_r & nrti_r & !pi_r ~ 'NRTI')) %>% 
	filter(!is.na(r)) %>%
	group_by(r) %>%
	summarise(n=n()) %>%
	arrange(-n) %>%
	mutate(idx = seq(1,n()),
		class_1 = str_split(r, ',', simplify=T)[,1],
		class_2 = str_split(r, ', ', simplify=T)[,2],
		class_1 = if_else(r == 'NNRTI, NRTI, PI', 'all', class_1))


treat_dat =  hiv_dr_cat %>% 
	filter(!pre_treatment & !is.na(nnrti) & !is.na(nrti) & !is.na(pi)) %>%
	mutate(
		nnrti_r = nnrti == 'intermediate/high',
		nrti_r = nrti == 'intermediate/high',
		pi_r = pi == 'intermediate/high',
		r = case_when(
			nnrti_r & nrti_r & pi_r ~ 'NNRTI, NRTI, PI',
			nnrti_r & nrti_r & !pi_r ~ 'NNRTI, NRTI',
			nnrti_r & !nrti_r & pi_r ~ 'NNRTI, PI',
			nnrti_r & !nrti_r & !pi_r ~ 'NNRTI',
			!nnrti_r & nrti_r & pi_r ~ 'NRTI, PI',
			!nnrti_r & !nrti_r & pi_r ~ 'PI',
			!nnrti_r & nrti_r & !pi_r ~ 'NRTI')) %>% 
	filter(!is.na(r)) %>%
	group_by(r) %>%
	summarise(n=n()) %>%
	arrange(-n) %>%
	mutate(idx = seq(1,n()),
		class_1 = str_split(r, ',', simplify=T)[,1],
		class_2 = str_split(r, ', ', simplify=T)[,2],
		class_1 = if_else(r == 'NNRTI, NRTI, PI', 'all', class_1))

max_y = 10*ceiling(max(c(treat_dat$n, pretreat_dat$n))/10)

a = ggplot(pretreat_dat, aes(x=idx, y = n, group=r, pattern=r, fill=class_1, pattern_fill=class_2, label = n)) +
	geom_col_pattern(pattern_spacing=0.025, pattern_density=0.5, color='#333333') +
	geom_text(color='#333333', vjust = -max_y*0.005) +
	scale_x_continuous(breaks=pretreat_dat$idx, labels=pretreat_dat$r) +
	scale_pattern_manual(values=patterns, guide='none') +
	scale_fill_manual(values=colors, guide='none') +
	scale_pattern_fill_manual(values=colors, guide='none') +
	gtheme +
	xlab('resistance class') +
	ylab('number of samples') +
	ylim(c(0, max_y)) +
	ggtitle('viremic pre-treatment') +
	theme(
		panel.grid.major.x = element_blank(),
		  panel.grid.minor.x = element_blank(),
		  axis.text.x = element_text(angle=90, size=10, vjust = 0.5, hjust=1))

b = ggplot(treat_dat, aes(x=idx, y = n, group=r, pattern=r, fill=class_1, pattern_fill=class_2, label=n)) +
	geom_col_pattern(pattern_spacing=0.025, pattern_density=0.5, color='#333333') +
	geom_text(color='#333333', vjust = -max_y*0.005) +
	scale_x_continuous(breaks=treat_dat$idx, labels=treat_dat$r) +
	scale_pattern_manual(values=patterns, guide='none') +
	scale_fill_manual(values=colors, guide='none') +
	scale_pattern_fill_manual(values=colors, guide='none') +
	gtheme +
	xlab('resistance class') +
	ylab('') +
	ylim(c(0, max_y)) +
	ggtitle('viremic treatment-experienced') +
	theme(
		panel.grid.major.x = element_blank(),
		  panel.grid.minor.x = element_blank(), 
		  axis.text.x = element_text(angle=90, size=10, vjust = 0.5, hjust=1))

prev_dat = read_tsv(prev_dat_file) %>%
	mutate(
		type = case_when(
			type == 'multi' ~ 'multi-class',
			type == 'single' ~ 'single-class'),
		type = ordered(type, levels=c('single-class', 'multi-class')))

c = ggplot(prev_dat, aes(x=type, y=fit*100)) + 
	geom_errorbar(aes(ymin=lwr*100, ymax=upr*100), color='#333333', size=0.5, width=0) +
	geom_point(fill='#eaeaea', color='#333333', shape=21, size=5) +
	ggtitle('round 18 (2017)') +
	xlab('NNRTI') + 
	ylab('prevalence among viremic PLHIV (%)') +
	ylim(0,NA) +
	gtheme + 
	theme(
		panel.grid.major.x = element_blank(),
		  panel.grid.minor.x = element_blank())

p = plot_grid(a, b, c, nrow=1, labels=c('a', 'b', 'c'), rel_widths=c(0.4, 0.4, 0.25))
ggsave('figures/final/figure_2.pdf', p, width=15, height=6)
