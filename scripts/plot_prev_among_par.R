suppressMessages(require(tidyverse))
suppressMessages(require(patchwork))
#library(cowplot)
suppressMessages(source('scripts/utils.R'))


dat_file = 'data/rakai_drug_resistance_categorized.tsv'
hiv_pred_file = 'models/HIV_prev_pred.tsv'
viremic_pred_file = 'models/HIV_viremic_prev_pred.tsv'
pretreat_pred_file = 'models/HIV_pretreat_prev_pred.tsv'
treat_pred_file = 'models/HIV_treat_prev_pred.tsv'

hiv_dr_cat = read_tsv(dat_file, show_col_types=FALSE)
round_dates = hiv_dr_cat %>% select(round, round_mid_date) %>% unique() %>% mutate(round = as.character(round))

#### PANEL A ####
# prevalene of hiv over time
dodge_a = c('treatment-experienced viremic'=37.5, 
	"pre-treatment viremic" = 12.5, 
	"viremic" = -12.5, 
	"all PLHIV" = -37.5)
order_a = c('treatment-experienced viremic'=4, 
	"pre-treatment viremic" = 3, 
	"viremic" = -2, 
	"all PLHIV" = -1)

hiv_pred = read_tsv(hiv_pred_file, show_col_types=FALSE)
viremic_pred = read_tsv(viremic_pred_file, show_col_types=FALSE)
pretreat_pred = read_tsv(pretreat_pred_file, show_col_types=FALSE)
treat_pred = read_tsv(treat_pred_file, show_col_types=FALSE)
dat = 
	bind_rows(
		hiv_pred %>% mutate(type = 'all PLHIV', round = as.character(round)),
		viremic_pred %>% filter(round > 15) %>% mutate(type = 'viremic', round = as.character(round)),
		pretreat_pred %>% filter(round > 14) %>% mutate(type = 'pre-treatment viremic', round = as.character(round)),
		treat_pred %>% filter(round > 15) %>% mutate(type = 'treatment-experienced viremic', round = as.character(round))) %>%
	left_join(round_dates, by='round') %>%
	mutate(type = ordered(type, 
		levels=c('all PLHIV', 'viremic', 
			'pre-treatment viremic', 'treatment-experienced viremic'))) %>%
	mutate(round_mid_date_dodge = round_mid_date + dodge_a[type]) %>%
	arrange(by=order_a[type])

a = ggplot(dat, aes(x=round_mid_date_dodge, y=fit*100, fill=type, color=type, group=type)) +
	geom_line(color='#333333') + 
	geom_errorbar(aes(ymin=lwr*100, ymax=upr*100), color='#333333', linewidth=1, width=0) +
	geom_errorbar(aes(ymin=lwr*100, ymax=upr*100, color=type), width=0) +
	geom_point(color='#333333', shape=21, size=4) +
	scale_fill_manual(values=colors, name=NULL) +
	scale_color_manual(values=colors) +
	guides(color = "none") +
	xlab('date') + 
	ylab('prevalence among\nRCCS participants (%)') +
	ylim(c(NA, max(dat$upr)*1.2*100)) +
	gtheme +
	theme(legend.position=c(1, 1),
		legend.justification=c(1,1),
		plot.margin = margin(0.3,0.3,0.3,0.3, "cm"))

#### PANEL B ####
# population prevalence of viremic resistance by round
dodge = c('nnrti'=25, "nrti" = 0, "pi" = -25)
order = c('nnrti' = 3, "nrti" = 2, "pi" = 1)
max_y = max(read_tsv('models/dr_amongPar_prev_pred_stratified.tsv', show_col_types=FALSE)$upr) *100

dr_among_par_file = 'models/dr_amongPar_prev_pred.tsv'
dr_among_par = read_tsv(dr_among_par_file, show_col_types=FALSE) %>%
		mutate(round = as.character(round)) %>%
	left_join(round_dates, by="round") %>%
	mutate(round_mid_date_dodge = round_mid_date + dodge[class]) %>%
	arrange(by=order[class])

b = ggplot(dr_among_par, 
			aes(x=round_mid_date_dodge, y=fit*100, color=class, shape=class, fill=class, group=class)) +
		geom_line(color='#333333') +
		geom_errorbar(aes(ymin=lwr*100, ymax=upr*100), 
			color='#333333', linewidth=1, width=0) +
		geom_errorbar(aes(ymin=lwr*100, ymax=upr*100, color=class), width=0) +
		geom_point(color='#333333', size=4) +
		scale_shape_manual(values=shapes, name=NULL) +
		scale_fill_manual(values=colors, name=NULL) +
		scale_color_manual(values=colors, guide="none") +
		guides(override.aes = list(size = 5)) +
		xlab('date') + 
		xlim(c(min(dr_among_par$round_mid_date)-180, 
			max(dr_among_par$round_mid_date)+180)) +
		ylab('population prevalence (%)') +
		ylim(c(0, max_y)) +
		ggtitle('viremic resistance\namong RCCS participants') +
		gtheme +
		theme(legend.justification = c(0, 1), legend.position = c(0, 1),
			plot.margin = margin(0.3,0.3,0.3,0.3, "cm"))


#### PANEL C ####
# prevalence of pre-treatment resistance over time among all participants
dodge = c('nnrti'=25, "nrti" = 0, "pi" = -25)
order = c('nnrti' = 3, "nrti" = 2, "pi" = 1)
dr_among_par_file = 'models/dr_amongPar_prev_pred.tsv'
max_y = 0.85

pretreat_among_par_file = 'models/pretreat_amongPar_prev_pred.tsv'
pretreat_among_par = read_tsv(pretreat_among_par_file, show_col_types=FALSE) %>%
		mutate(round = as.character(round)) %>%
	left_join(round_dates, by="round") %>%
	mutate(round_mid_date_dodge = round_mid_date + dodge[class]) %>%
	arrange(by=order[class])

c = ggplot(pretreat_among_par, 
			aes(x=round_mid_date_dodge, y=fit*100, color=class, shape=class, fill=class, group=class)) +
		geom_line(color='#333333') +
		geom_errorbar(aes(ymin=lwr*100, ymax=upr*100), 
			color='#333333', size=1, width=0) +
		geom_errorbar(aes(ymin=lwr*100, ymax=upr*100, color=class), width=0) +
		geom_point(color='#333333', size=4) +
		scale_shape_manual(values=shapes, name=NULL) +
		scale_fill_manual(values=colors, name=NULL) +
		scale_color_manual(values=colors, guide="none") +
		guides(override.aes = list(size = 5)) +
		xlab('date') + 
		ggtitle('pre-treatment viremic resistance\namong RCCS participants') +
		ylab('population prevalence (%)') +
		ylim(c(0, max_y)) +
		gtheme +
		theme(legend.justification = c(0, 1), legend.position = c(0, 1),
			plot.margin = margin(0.3,0.3,0.3,0.3, "cm"))


#### PANEL D ####
# prevalence of treatment experienced viremia among participants
treat_among_par_file = 'models/treat_amongPar_prev_pred.tsv'
treat_among_par = read_tsv(treat_among_par_file, show_col_types=FALSE) %>% 
	mutate(round = as.character(round)) %>% left_join(round_dates) %>%
	mutate(round_mid_date_dodge = round_mid_date + dodge[class]) %>%
	arrange(by=order[class])

d = ggplot(treat_among_par, 
			aes(x=round_mid_date_dodge, y=fit*100, color=class, shape=class, fill=class, group=class)) +
		geom_line(color='#333333') +
		geom_errorbar(aes(ymin=lwr*100, ymax=upr*100), color='#333333', size=1, width=0) +
		geom_errorbar(aes(ymin=lwr*100, ymax=upr*100, color=class), width=0) +
		geom_point(color='#333333', size=4) +
		scale_shape_manual(values=shapes, name=NULL) +
		scale_fill_manual(values=colors, name=NULL) +
		scale_color_manual(values=colors, guide="none") +
		xlim(c(min(treat_among_par$round_mid_date)-180, 
			max(treat_among_par$round_mid_date)+180)) +
		xlab('date') + 
		ggtitle('treatment-experienced\nviremic resistance\namong RCCS participants') +
		ylab('population prevalence (%)') +
		ylim(c(0, max_y)) +
		gtheme +
		theme(axis.title.y = element_text(size=14),
			legend.justification = c(0, 1), legend.position = c(0, 1))


p1 = a + b + plot_layout(widths = c(1,0.5))

p2 = c + d + plot_layout(widths=c(1, 0.5)) 

p = p1 / p2  + 
  plot_annotation(tag_levels = 'A')  & 
  theme(plot.tag = element_text(face="bold"))

if (!dir.exists('figures')){
dir.create('figures')
} 
ggsave('figures/prev_among_par.pdf', p, width=10, height=11)
