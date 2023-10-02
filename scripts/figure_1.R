library(tidyverse)
source('scripts/utils.R')
library(cowplot)

dat_file = 'data/rakai_drug_resistance_categorized.tsv'
hiv_pred_file = 'models/HIV_prev_pred.tsv'
viremic_pred_file = 'models/HIV_viremic_prev_pred.tsv'
pretreat_pred_file = 'models/HIV_pretreat_viremic_prev_pred.tsv'
treat_pred_file = 'models/HIV_treat_viremic_prev_pred.tsv'

hiv_dr_cat = read_tsv(dat_file)
round_dates = hiv_dr_cat %>% select(round, round_mid_date) %>% unique() %>% mutate(round = as.character(round))

#### PANEL A ####
# prevalene of hiv over time
dodge_a = c('viremic treatment-experienced'=37.5, 
	"viremic pre-treatment" = 12.5, 
	"viremic" = -12.5, 
	"all PLHIV" = -37.5)
order_a = c('viremic treatment-experienced'=4, 
	"viremic pre-treatment" = 3, 
	"viremic" = -2, 
	"all PLHIV" = -1)

hiv_pred = read_tsv(hiv_pred_file)
viremic_pred = read_tsv(viremic_pred_file)
pretreat_pred = read_tsv(pretreat_pred_file)
treat_pred = read_tsv(treat_pred_file)
dat = 
	bind_rows(
		hiv_pred %>% mutate(type = 'all PLHIV', round = as.character(round)),
		viremic_pred %>% filter(round > 15) %>% mutate(type = 'viremic', round = as.character(round)),
		pretreat_pred %>% filter(round > 14) %>% mutate(type = 'viremic pre-treatment', round = as.character(round)),
		treat_pred %>% filter(round > 15) %>% mutate(type = 'viremic treatment-experienced', round = as.character(round))) %>%
	left_join(round_dates, by='round') %>%
	mutate(type = ordered(type, 
		levels=c('all PLHIV', 'viremic', 
			'viremic pre-treatment', 'viremic treatment-experienced'))) %>%
	mutate(round_mid_date_dodge = round_mid_date + dodge_a[type]) %>%
	arrange(by=order_a[type])


a = ggplot(dat, aes(x=round_mid_date_dodge, y=fit*100, fill=type, color=type, group=type)) +
	geom_line(color='#333333') + 
	geom_errorbar(aes(ymin=lwr*100, ymax=upr*100), color='#333333', linewidth=1, width=0) +
	geom_errorbar(aes(ymin=lwr*100, ymax=upr*100, color=type), width=0) +
	geom_point(color='#333333', shape=21, size=4) +
	scale_fill_manual(values=colors, name='') +
	scale_color_manual(values=colors) +
	guides(color = "none") +
	xlab('date') + 
	ylab('population prevalence in RCCS (%)') +
	ylim(c(NA, max(dat$upr)*1.2*100)) +
	gtheme +
	theme(legend.position=c(0.62, 0.9), legend.text=element_text(size=14),
		plot.margin = margin(0.3,0.3,0.3,0.3, "cm"))



#### PANEL B ####
# prevalence of pre-treatment resistance over time among all participants
dodge = c('nnrti'=25, "nrti" = 0, "pi" = -25)
order = c('nnrti' = 3, "nrti" = 2, "pi" = 1)
all_among_Par_file = 'models/all_amongPar_prev_pred.tsv'
max_y = max(read_tsv(all_among_Par_file)$upr) *100

pretreat_amongPar_file = 'models/pretreat_amongPar_prev_pred.tsv'
pretreat_amongPar = read_tsv(pretreat_amongPar_file) %>%
		mutate(round = as.character(round)) %>%
	left_join(round_dates, by="round") %>%
	mutate(round_mid_date_dodge = round_mid_date + dodge[class]) %>%
	arrange(by=order[class])

b = ggplot(pretreat_amongPar, 
			aes(x=round_mid_date_dodge, y=fit*100, color=class, shape=class, fill=class, group=class)) +
		geom_line(color='#333333') +
		geom_errorbar(aes(ymin=lwr*100, ymax=upr*100), 
			color='#333333', size=1, width=0) +
		geom_errorbar(aes(ymin=lwr*100, ymax=upr*100, color=class), width=0) +
		geom_point(color='#333333', size=4) +
		scale_shape_manual(values=shapes, name='') +
		scale_fill_manual(values=colors, name='') +
		scale_color_manual(values=colors, guide="none") +
		guides(override.aes = list(size = 5)) +
		xlab('date') + 
		ylab('population prevalence of\nviremic pre-treatment resistance (%)') +
		ylim(c(0, max_y)) +
		gtheme +
		theme(legend.justification = c(0, 1), legend.position = c(0, 1),
			plot.margin = margin(0.3,0.3,0.3,0.3, "cm"))


#### PANEL C ####
# prevalence of treatment experienced viremia among participants
treat_among_Par_file = 'models/treat_amongPar_prev_pred.tsv'
treat_amongPar = read_tsv(treat_among_Par_file) %>% 
	mutate(round = as.character(round)) %>% left_join(round_dates) %>%
	mutate(round_mid_date_dodge = round_mid_date + dodge[class]) %>%
	arrange(by=order[class])

c = ggplot(treat_amongPar, 
			aes(x=round_mid_date_dodge, y=fit*100, color=class, shape=class, fill=class, group=class)) +
		geom_line(color='#333333') +
		geom_errorbar(aes(ymin=lwr*100, ymax=upr*100), color='#333333', size=1, width=0) +
		geom_errorbar(aes(ymin=lwr*100, ymax=upr*100, color=class), width=0) +
		geom_point(color='#333333', size=4) +
		scale_shape_manual(values=shapes, guide='none') +
		scale_fill_manual(values=colors, name='') +
		scale_color_manual(values=colors) +
		guides(color = "none", fill="none") +
		xlim(c(min(treat_amongPar$round_mid_date)-180, 
			max(treat_amongPar$round_mid_date)+180)) +
		xlab('date') + 
		ylab('population prevalence of viremic\ntreatment-experienced resistance (%)') +
		ylim(c(0, max_y)) +
		gtheme +
		theme(axis.title.y = element_text(size=14))


#### PANEL D ####
# prevalence of pre-treatment resistance over time 
pretreat_amongPLHIV_file = 'models/pretreat_amongPLWHIV_prev_pred.tsv'
pretreat_amongPLHIV = read_tsv(pretreat_amongPLHIV_file) %>%
	mutate(round = as.character(round)) %>%
	left_join(round_dates, by="round") %>%
	mutate(round_mid_date_dodge = round_mid_date + dodge[class]) %>%
	arrange(by=order[class])

d = ggplot(pretreat_amongPLHIV, 
			aes(x=round_mid_date_dodge, y=fit*100, color=class, shape=class, fill=class, group=class)) +
		geom_line(color='#333333') +
		geom_errorbar(aes(ymin=lwr*100, ymax=upr*100), color='#333333', size=1, width=0) +
		geom_errorbar(aes(ymin=lwr*100, ymax=upr*100, color=class), width=0) +
		geom_point(color='#333333', size=4) +
		scale_shape_manual(values=shapes, name='', guide='none') +
		scale_fill_manual(values=colors, name='') +
		scale_color_manual(values=colors) +
		guides(color = "none", fill="none") +
		xlab('date') + 
		ylab('prevalence of resistance among\nviremic pre-treatment PLHIV (%)') +
		gtheme +
		theme(legend.justification = c(0, 1), legend.position = c(0, 1),legend.text=element_text(size=16),
			plot.margin = margin(0.3,0.3,0.3,0.3, "cm"))



#### PANEL E ####
# prevalence of resistance among treatment experienced PLHIV
treat_among_PLHIV_file = 'models/treat_amongPLWHIV_prev_pred.tsv'
treat_amongPLHIV = read_tsv(treat_among_PLHIV_file) %>% 
	mutate(round = as.character(round)) %>%
	left_join(round_dates) %>%
	mutate(round_mid_date_dodge = round_mid_date + dodge[class]) %>%
	arrange(by=order[class])

e = ggplot(treat_amongPLHIV, 
			aes(x=round_mid_date_dodge, y=fit*100, color=class, shape=class, fill=class, group=class)) +
		geom_line(color='#333333') +
		geom_errorbar(aes(ymin=lwr*100, ymax=upr*100), color='#333333', size=1, width=0) +
		geom_errorbar(aes(ymin=lwr*100, ymax=upr*100, color=class), width=0) +
		geom_point(color='#333333', size=4) +
		scale_shape_manual(values=shapes, guide='none') +
		scale_fill_manual(values=colors, name='') +
		scale_color_manual(values=colors) +
		guides(color = "none", fill="none") +
		xlim(c(min(treat_amongPLHIV$round_mid_date)-180, 
			max(treat_amongPLHIV$round_mid_date)+180)) +
		xlab('date') + 
		ylab('prevalence of resistance among viremic\ntreatment-experienced PLHIV (%)') +
		gtheme 


p1 = plot_grid(a, b, c, labels = c("a", "b", "c"), nrow=1, rel_widths=c(0.5, 0.5, 0.35))
p2 = plot_grid(NULL, d,e, NULL, labels = c('', "d", "e", ''), nrow=1, rel_widths=c(0.25, 0.5, 0.35, 0.25))
p = plot_grid(p1, p2, nrow=2, rel_widths=c(1, 0.5))

ggsave('figures/final/figure_1.pdf', p, width=16, height=10)

ggsave('figures/final/figure_1D.pdf', d, width=6.4, height=4.8)

