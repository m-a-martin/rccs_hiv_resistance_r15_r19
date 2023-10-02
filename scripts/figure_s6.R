library(tidyverse)
source('scripts/utils.R')
library(cowplot)

t97a_file = 'models/pretreat_int97a_amongPLWHIV_prev_pred.tsv'
t97a = read_tsv(t97a_file)

dat_file = 'data/rakai_drug_resistance_categorized.tsv'
mut_file = 'data/rakai_drug_resistance_mut_formatted.tsv'

hiv_mut = read_tsv(mut_file)
hiv_dr_cat = read_tsv(dat_file)


round_dates = hiv_dr_cat %>% select(round, round_mid_date) %>% unique()

t97a = t97a %>% left_join(round_dates, by="round")

max_y = max(t97a %>% filter(s == 'sex' | s == 'age_cat' | s == 'comm_type' | s == 'subtype') %>% select(upr))*100


#### PANEL A ####
# prevalence of T97A by sex
dodge = c('F'=12.5, "M" = -12.5)
order = c('F' = 2, 'M' = 1)
a_dat = t97a %>% filter(s == 'sex') %>%
	mutate(round_mid_date_dodge = round_mid_date + dodge[s_var]) %>%
	arrange(by=order[s_var])

a = ggplot(a_dat,
		aes(x=round_mid_date_dodge, y=fit*100, fill=s_var, color=s_var, group=s_var)) +
	geom_line(color='#333333') + 
	geom_errorbar(aes(ymin=lwr*100, ymax=upr*100), color='#333333', size=1, width=0) +
	geom_errorbar(aes(ymin=lwr*100, ymax=upr*100, color=s_var), width=0) +
	geom_point(color='#333333', shape=21, size=4) +
	scale_fill_manual(values=colors, labels=var_rename, name='sex') +
	scale_color_manual(values=colors) +
	guides(color = "none") +
	xlab('') + 
	ylab('prevalence of inT97A among viremic\npre-treatment PLWHIV (%)') +
	ylim(c(0, max_y)) + 
	gtheme +
	theme(legend.justification = c(0, 1), legend.position = c(0, 1), 
		legend.text=element_text(size=14),
		plot.margin = margin(0.3,0.3,0.3,0.3, "cm"))

#### PANEL B ####
# NNNRTI pre-treatment resistance by round and age
dodge = c("(14,24]" = -25, "(24,34]" = 0, '(34,50]'=25)
order = c('(14,24]' = 1, "(24,34]"  = 2, '(34,50]' = 3)
b_dat = t97a %>% filter(s == 'age_cat') %>%
	mutate(round_mid_date_dodge = round_mid_date + dodge[s_var]) %>%
	arrange(by=order[s_var])

b = ggplot(b_dat,
		aes(x=round_mid_date_dodge, y=fit*100, fill=s_var, color=s_var, group=s_var)) +
	geom_line(color='#333333') + 
	geom_errorbar(aes(ymin=lwr*100, ymax=upr*100), color='#333333', size=1, width=0) +
	geom_errorbar(aes(ymin=lwr*100, ymax=upr*100, color=s_var), width=0) +
	geom_point(color='#333333', shape=21, size=4) +
	scale_fill_manual(values=colors, labels=var_rename, name='age') +
	scale_color_manual(values=colors) +
	guides(color = "none") +
	xlab('') + 
	ylab('\n') +
	ylim(c(0, max_y)) + 
	gtheme +
	theme(legend.justification = c(0,1), legend.position = c(0, 1), 
		legend.text=element_text(size=14),
		legend.title = element_text(size=14),
		plot.margin = margin(0.3,0.3,0.3,0.3, "cm"))

#### PANEL C ####
# NNRTIpre-treatment resistance by round and community type
dodge = c("Agrarian" = -25, "Fishing" = 0, 'Trading'=25)
order = c('Agrarian' = 1, "Fishing"  = 2, 'Trading' = 3)
c_dat = t97a %>% filter(s == 'comm_type') %>%
	mutate(round_mid_date_dodge = round_mid_date + dodge[s_var]) %>%
	arrange(by=order[s_var])

c = ggplot(c_dat,
		aes(x=round_mid_date_dodge, y=fit*100, fill=s_var, color=s_var, group=s_var)) +
	geom_line(color='#333333') + 
	geom_errorbar(aes(ymin=lwr*100, ymax=upr*100), color='#333333', size=1, width=0) +
	geom_errorbar(aes(ymin=lwr*100, ymax=upr*100, color=s_var), width=0) +
	geom_point(color='#333333', shape=21, size=4) +
	scale_fill_manual(values=colors, labels=var_rename, name='community type') +
	scale_color_manual(values=colors) +
	guides(color = "none") +
	xlab('date') + 
	ylab('prevalence of inT97A among viremic\npre-treatment PLWHIV (%)') +
		ylim(c(0, max_y)) + 
	gtheme +
	theme(legend.justification = c(0,1), legend.position = c(0, 1), 
		legend.text=element_text(size=14),
		legend.title = element_text(size=14),
		plot.margin = margin(0.3,0.3,0.3,0.3, "cm"))


#### PANEL D ####
# NNRTIpre-treatment resistance by round and subtype
dodge = c("A1" = -25, "D" = 0, 'other (recombinant)'=25)
order = c('A1' = 1, "D"  = 2, 'ther (recombinant)' = 3)
d_dat = t97a %>% filter(s == 'subtype' & s_var != 'other (non-recombinant)' & s_var != 'B' & s_var != 'C') %>%
	mutate(round_mid_date_dodge = round_mid_date + dodge[s_var]) %>%
	arrange(by=order[s_var])

d = ggplot(d_dat,
		aes(x=round_mid_date_dodge, y=fit*100, fill=s_var, color=s_var, group=s_var)) +
	geom_line(color='#333333') + 
	geom_errorbar(aes(ymin=lwr*100, ymax=upr*100), color='#333333', size=1, width=0) +
	geom_errorbar(aes(ymin=lwr*100, ymax=upr*100, color=s_var), width=0) +
	geom_point(color='#333333', shape=21, size=4) +
	scale_fill_manual(values=colors, labels=var_rename, name='subtype') +
	scale_color_manual(values=colors) +
	guides(color = "none") +
	xlab('date') + 
	ylab('\n') +
		ylim(c(0, max_y)) + 
	gtheme +
	theme(legend.justification = c(0,1), legend.position = c(0, 1), 
		legend.text=element_text(size=14),
		legend.title = element_text(size=14),
		plot.margin = margin(0.3,0.3,0.3,0.3, "cm"))


p2 = plot_grid(a, b, c, d, labels=c('a', 'b', 'c', 'd'), nrow=2)
p = plot_grid(title, p2, nrow=2, rel_heights=c(0.1,0.9))
ggsave('figures/final/figure_s6.pdf', p2, width=12, height=9.6)


