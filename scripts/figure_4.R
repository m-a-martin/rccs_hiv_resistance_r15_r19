library(tidyverse)
source('scripts/utils.R')
library(cowplot)
library(ggrepel)

dat_file = 'data/rakai_drug_resistance_mut_formatted.tsv'
dr_cat_file = 'data/rakai_drug_resistance_categorized.tsv'

hiv_dr_cat = read_tsv(dr_cat_file)
dat = read_tsv(dat_file)

#### PANEL A ####
# number of samples each mutation is observed in
dr_muts = dat %>% group_by(mut) %>% summarise(n=n()) %>%
	mutate(protein = str_replace_all(mut, "[[A-Z][0-9]]", ""),
			sub = str_replace_all(mut, "[[a-z]]", "")) %>%
	arrange(protein, -n) %>%
	mutate(idx = seq(1,n()))

blank_dat = tibble(protein=unique(dr_muts$protein)) %>%
	mutate(
		idx=as.numeric(NA), 
		x=as.numeric(NA), 
		y=as.numeric(NA),
		n=as.numeric(NA),
		idx=as.numeric(NA))

a = ggplot(dr_muts, aes(x=idx, y=n, fill=protein, shape=protein)) +
	geom_bar(stat="identity", position="identity", color='#333333') +
	geom_point(data = blank_dat, color='#333333', aes(fill=protein), size=4) +
	scale_x_continuous(breaks = dr_muts$idx,labels = dr_muts$sub, expand = c(0.01,0.01)) +
	xlab('substitution') +
	ylab('number of samples') + 
	scale_fill_manual(values=colors,
		labels=var_rename, name='') +
	scale_shape_manual(values=shapes, labels=var_rename, name='') +
	gtheme +
	theme(axis.text.x = element_text(angle=90, size=7, vjust = 0.5, hjust=1),
		panel.grid.minor.x = element_blank(),
		panel.grid.major.x = element_blank(),
		legend.text = element_text(size=14),
		legend.justification = c(1, 1), legend.position = c(1, 1)) 

#### PANEL B ####
# risk ratio of 15 most frequent mutations in pretreated data compared to R15
pretreat_dr_muts_prev_file = 'models/pretreat_dr_muts_amongPLWHIV_prev_pred.tsv'
pretreat_dr_muts_prev = read_tsv(pretreat_dr_muts_prev_file)
plot_muts = pretreat_dr_muts_prev %>% filter(round == 18) %>% arrange(-fit) %>%
	slice(1:15) %>% select(mut) %>% mutate(idx = seq(1,n()))

pretreat_dr_muts_rr_file = 'models/pretreat_dr_muts_amongPLWHIV_rr_pred.tsv'
pretreat_dr_muts_rr = read_tsv(pretreat_dr_muts_rr_file)
print(pretreat_dr_muts_rr)
pretreat_dr_muts_rr = pretreat_dr_muts_rr %>% inner_join(plot_muts, by="mut") %>%
	filter(var == 'round18') %>%
	mutate(protein = str_replace_all(mut, "[[A-Z][0-9]]", ""),
			sub = str_replace_all(mut, "[[a-z]]", "")) %>%
	arrange(-RR)


y_breaks = 2^seq(0, ceiling(log2(max(pretreat_dr_muts_rr$UCI))),1)


b = ggplot(pretreat_dr_muts_rr, aes(x=idx, y=RR, fill=protein, color=protein, shape=protein)) +
	geom_hline(yintercept=1, linetype="dashed", 
                color = "#333333", size=1) +
	geom_errorbar(aes(ymin=LCI, ymax=UCI), color='#333333', linewidth=1, width=0) +
	geom_errorbar(aes(ymin=LCI, ymax=UCI), width=0) +
	geom_point(color='#333333', size=4) +
	scale_x_continuous(breaks = pretreat_dr_muts_rr$idx,labels = rep('', nrow(pretreat_dr_muts_rr)), expand = c(0.02,0.02)) +
	scale_y_continuous(trans='log2', breaks=y_breaks) +
	scale_fill_manual(values=colors, guide='none') +
	scale_color_manual(values=colors, guide='none') +
	scale_shape_manual(values=shapes, guide='none') +
	xlab('') +
	ylab('\nprevalence ratio among viremic\npre-treatment PLHIV (v. R15)') +
	gtheme +
	theme(
		axis.title.y = element_text(size=12),
		axis.text.x = element_text(angle=90, size=10, vjust = 0.5, hjust=1),
		panel.grid.minor.x = element_blank(),
		panel.grid.major.x = element_blank(),
		legend.justification = c(0, 1), legend.position = c(0, 1)) 


#### PANEL C ####
# prevalence of 15 most frequent mutations in pretreated data compared to R15
pretreat_dr_muts_prev = pretreat_dr_muts_prev %>% inner_join(plot_muts) %>%
	filter(round == 18) %>%
	mutate(protein = str_replace_all(mut, "[[A-Z][0-9]]", ""),
			sub = str_replace_all(mut, "[[a-z]]", ""))


c = ggplot(pretreat_dr_muts_prev, aes(x=idx, y=fit*100, fill=protein, color=protein, shape=protein)) +
	geom_errorbar(aes(ymin=lwr*100, ymax=upr*100), color='#333333', linewidth=1, width=0) +
	geom_errorbar(aes(ymin=lwr*100, ymax=upr*100), width=0) +
	geom_point(color='#333333', size=4) +
	scale_x_continuous(breaks = pretreat_dr_muts_rr$idx,labels = pretreat_dr_muts_rr$sub, expand = c(0.02,0.02)) +
	scale_fill_manual(values=colors, guide='none') +
	scale_color_manual(values=colors, guide='none') +
	scale_shape_manual(values=shapes, guide='none') +
	xlab('substitution') +
	ylab('\nprevalence among viremic\npre-treatment PLHIV (%)') +
	gtheme +
	theme(
		axis.title.y = element_text(size=12),
		axis.text.x = element_text(angle=90, size=10, vjust = 0.5, hjust=1),
		panel.grid.minor.x = element_blank(),
		panel.grid.major.x = element_blank(),
		legend.justification = c(0, 1), legend.position = c(0, 1)) 


#### PANEL D ####
# prevalence of 15 most frequent mutations in treated data compared to R15
treat_dr_muts_prev_file = 'models/treat_dr_muts_amongPLWHIV_prev_pred.tsv'
treat_dr_muts_prev = read_tsv(treat_dr_muts_prev_file)
plot_muts = treat_dr_muts_prev %>% filter(round == 18) %>% arrange(-fit) %>%
	slice(1:15) %>% select(mut) %>% mutate(idx = seq(1,n()))
treat_dr_muts_prev = treat_dr_muts_prev %>% inner_join(plot_muts, by='mut') %>%
	filter(round == 18) %>%
	mutate(protein = str_replace_all(mut, "[[A-Z][0-9]]", ""),
			sub = str_replace_all(mut, "[[a-z]]", ""))

d = ggplot(treat_dr_muts_prev, aes(x=idx, y=fit*100, fill=protein, color=protein, shape=protein)) +
	geom_errorbar(aes(ymin=lwr*100, ymax=upr*100), color='#333333', linewidth=1, width=0) +
	geom_errorbar(aes(ymin=lwr*100, ymax=upr*100), width=0) +
	geom_point(color='#333333', size=4) +
	scale_x_continuous(breaks = treat_dr_muts_prev$idx,labels = treat_dr_muts_prev$sub, expand = c(0.02,0.02)) +
	scale_fill_manual(values=colors, guide='none') +
	scale_color_manual(values=colors, guide='none') +
	scale_shape_manual(values=shapes, guide='none') +
	xlab('substitution') +
	ylab('prevalence among viremic\ntreatment-experienced PLHIV (%)') +
	gtheme +
	theme(axis.text.x = element_text(angle=90, size=10, vjust = 0.5, hjust=1),
		panel.grid.minor.x = element_blank(),
		panel.grid.major.x = element_blank(),
		legend.justification = c(0, 1), legend.position = c(0, 1)) 




#### PANEL E ####
# prevalence correlation between treated and untreated individuals
treat_dr_muts_prev = read_tsv(treat_dr_muts_prev_file) %>% filter(round == 18)
pretreat_dr_muts_prev = read_tsv(pretreat_dr_muts_prev_file) %>% filter(round == 18)
dr_muts_prev = treat_dr_muts_prev %>% full_join(pretreat_dr_muts_prev, 
	by=c("mut", "round"), suffix=c("_treat", "_pretreat")) %>%
	mutate(
		fit_treat = replace_na(fit_treat, 0),
		lwr_treat = replace_na(lwr_treat, 0),
		upr_treat = replace_na(upr_treat, 0),
		fit_pretreat = replace_na(fit_pretreat, 0),
		lwr_pretreat = replace_na(lwr_pretreat, 0),
		upr_pretreat = replace_na(upr_pretreat, 0),
		protein = str_replace_all(mut, "[[A-Z][0-9]]", ""),
			sub = str_replace_all(mut, "[[a-z]]", "")) %>%
	mutate(label = case_when(
		mut %in% c("rtM184V", "rtK103N", "inT97A", "rtK65R", "rtY181C", "rtG190A") ~ mut))


e = ggplot(dr_muts_prev, aes(x=fit_pretreat*100, y=fit_treat*100, fill=protein, label=label, shape=protein)) +
	geom_abline(intercept = 0, color='#adadad', linetype="dashed") +
	geom_errorbar(aes(xmin=lwr_pretreat*100, xmax=upr_pretreat*100), color='#333333', size=0.25, width=0) +
	geom_errorbar(aes(ymin=lwr_treat*100, ymax=upr_treat*100), color='#333333', size=0.25, width=0) +
	geom_point(color='#333333', size=2, alpha=0.95) +
	geom_text_repel(color='#333333', box.padding = 1.2, nudge_x = 0.01,force=2, nudge_y = 0.01) +
	xlab('prevalence among viremic\npre-treatment PLHIV (%)') +
	ylab('prevalence among viremic\ntreatment-experienced PLHIV (%)') +
	scale_fill_manual(values=colors, guide="none") +
	scale_shape_manual(values=shapes, guide='none') +
	xlim(c(0, max(max(dr_muts_prev$upr_treat), max(dr_muts_prev$upr_pretreat))*100)) +
	ylim(c(0, max(max(dr_muts_prev$upr_treat), max(dr_muts_prev$upr_pretreat))*100)) +
	gtheme +
	coord_fixed()


p2 = plot_grid(b, c, labels=c("b", "c"), nrow=2)
bottom = plot_grid(p2, d, e, labels=c("", "d", "e"), nrow=1)
p = plot_grid(a, bottom, nrow=2, labels=c('a', ''), rel_heights=c(0.4, 0.6))

ggsave('figures/final/figure_4.pdf', p, width=14, height=10)


