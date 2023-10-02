library(tidyverse)
library(cowplot)
source('scripts/utils.R')

dat_file = 'data/rakai_drug_resistance_categorized.tsv'

hiv_dr_cat = read_tsv(dat_file)

round_dates = hiv_dr_cat %>% select(round, round_mid_date) %>% unique()

dat_file = 'models/pretreat_amongPLWHIV_prev_pred_stratified.tsv'
dat = read_tsv(dat_file) %>% filter(class == 'nnrti') %>% 
	left_join(round_dates, by="round") %>%
		mutate(fit = 100*fit, lwr=100*lwr, upr=100*upr)


max_y = max(dat %>% filter(s == 'sex' | s == 'age_cat' | s == 'comm_type') %>% select(upr))

#### PANEL A ####
# NNRTI pre-treatment resistance by round and sex
dodge = c('F'=12.5, "M" = -12.5)
order = c('F' = 2, 'M' = 1)
a_dat = dat %>% filter(s == 'sex') %>%
	mutate(round_mid_date_dodge = round_mid_date + dodge[s_var]) %>%
	arrange(by=order[s_var])

a = ggplot(a_dat,
		aes(x=round_mid_date_dodge, y=fit, fill=s_var, color=s_var, group=s_var)) +
	geom_line(color='#333333') + 
	geom_errorbar(aes(ymin=lwr, ymax=upr), color='#333333', linewidth=1, width=0) +
	geom_errorbar(aes(ymin=lwr, ymax=upr, color=s_var), width=0) +
	geom_point(color='#333333', shape=shapes['nnrti'], size=4) +
	scale_fill_manual(values=colors, labels=var_rename, name='sex') +
	scale_color_manual(values=colors) +
	guides(color = "none") +
	xlab('date') + 
	ylab('prevalence of resistance among\nviremic pre-treatment PLWHIV (%)') +
	ylim(c(0, max_y)) + 
	gtheme +
	theme(legend.justification = c(0, 1), legend.position = c(0, 1), 
		legend.text=element_text(size=14),
		plot.margin = margin(0.3,0.3,0.3,0.3, "cm"))

#### PANEL B ####
dodge = c("(14,24]" = -25, "(24,34]" = 0, '(34,50]'=25)
order = c('(14,24]' = 1, "(24,34]"  = 2, '(34,50]' = 3)
b_dat = dat %>% filter(s == 'age_cat') %>%
	mutate(round_mid_date_dodge = round_mid_date + dodge[s_var]) %>%
	arrange(by=order[s_var])

# NNNRTI pre-treatment resistance by round and age
b = ggplot(b_dat,
		aes(x=round_mid_date_dodge, y=fit, fill=s_var, color=s_var, group=s_var)) +
	geom_line(color='#333333') + 
	geom_errorbar(aes(ymin=lwr, ymax=upr), color='#333333', linewidth=1, width=0) +
	geom_errorbar(aes(ymin=lwr, ymax=upr, color=s_var), width=0) +
	geom_point(color='#333333', shape=shapes['nnrti'], size=4) +
	scale_fill_manual(values=colors, labels=var_rename, name='age') +
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

#### PANEL C ####
dodge = c("Agrarian" = -25, "Fishing" = 0, 'Trading'=25)
order = c('Agrarian' = 1, "Fishing"  = 2, 'Trading' = 3)
c_dat = dat %>% filter(s == 'comm_type') %>%
	mutate(round_mid_date_dodge = round_mid_date + dodge[s_var]) %>%
	arrange(by=order[s_var])

# NNRTIpre-treatment resistance by round and community type
c = ggplot(c_dat,
		aes(x=round_mid_date_dodge, y=fit, fill=s_var, color=s_var, group=s_var)) +
	geom_line(color='#333333') + 
	geom_errorbar(aes(ymin=lwr, ymax=upr), color='#333333', linewidth=1, width=0) +
	geom_errorbar(aes(ymin=lwr, ymax=upr, color=s_var), width=0) +
	geom_point(color='#333333', shape=shapes['nnrti'], size=4) +
	scale_fill_manual(values=colors, labels=var_rename, name='community type') +
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


p = plot_grid(a, b, c, labels=c('a', 'b', 'c'), nrow=1)
ggsave('figures/final/figure_s4.pdf', p, width=16, height=4.8)
