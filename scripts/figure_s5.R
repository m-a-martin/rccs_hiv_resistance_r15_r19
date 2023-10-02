library(tidyverse)
source('scripts/utils.R')
library(cowplot)
library(ggpattern)

dat_file = 'data/rakai_drug_resistance_categorized.tsv'

hiv_dr_cat = read_tsv(dat_file) %>% filter(round == 17 | round == 18)

r18_dat =  hiv_dr_cat %>% 
	filter(!is.na(nnrti) & !is.na(nrti) & !is.na(pi)) %>%
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
	filter(round == 18 & (nnrti_r | nrti_r | pi_r)) %>%
	group_by(r) %>%
	summarise(n=n()) %>%
	arrange(-n) %>%
	mutate(idx = seq(1,n()),
		class_1 = str_split(r, ',', simplify=T)[,1],
		class_2 = str_split(r, ', ', simplify=T)[,2],
		class_1 = if_else(r == 'NNRTI, NRTI, PI', 'all', class_1))

max_y = 10*ceiling(max(c(r18_dat$n))/10)

a = ggplot(r18_dat, aes(x=idx, y = n, group=r, pattern=r, fill=class_1, pattern_fill=class_2, label = n)) +
	geom_col_pattern(pattern_spacing=0.025, pattern_density=0.5, color='#333333') +
	geom_text(color='#333333', vjust = -max_y*0.005) +
	scale_x_continuous(breaks=r18_dat$idx, labels=r18_dat$r) +
	scale_pattern_manual(values=patterns, guide='none') +
	scale_fill_manual(values=colors, guide='none') +
	scale_pattern_fill_manual(values=colors, guide='none') +
	gtheme +
	xlab('resistance class') +
	ylab('number of samples') +
	ggtitle('viremic round 18') +
	theme(
		panel.grid.major.x = element_blank(),
		  panel.grid.minor.x = element_blank(),
		  axis.text.x = element_text(size=10))


ggsave('figures/final/figure_s5.pdf', a, width=6.8, height=4.8)
