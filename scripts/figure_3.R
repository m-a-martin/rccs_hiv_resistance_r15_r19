library(tidyverse)
library(cowplot)
library(ggrepel)
library(patchwork)
source('scripts/utils.R')

dat_file = 'data/rakai_drug_resistance_formatted.tsv'
comm_type_prev_file = 'models/all_amongPar_prev_pred_stratified.tsv'
comm_num_prev_file = 'models/all_amongPar_comm_num_R18_prev_pred.tsv'
comm_v_prev_file = 'models/HIV_viremic_prev_pred_stratified.tsv'
class_prev_pred_file = 'models/all_amongPLWHIV_prev_pred_link.tsv'

comm_type_prev = read_tsv(comm_type_prev_file, show_col_types=FALSE) %>% 
	filter(round == 18, s == 'comm_type') %>%
	rename(comm_num=s_var) %>%
	mutate(
		type='pool',
		comm_num = case_when(
			comm_num == 'Fishing' ~ 'fishing',
			comm_num == 'Agrarian' ~ 'agrarian',
			comm_num == 'Trading' ~ 'trading')) %>%
	select(-s)

comm_num_prev = read_tsv(comm_num_prev_file, show_col_types=FALSE) %>% filter(round == 18) %>% 
	mutate(comm_num = as.character(comm_num), type='ind')

comm_prev = bind_rows(comm_type_prev, comm_num_prev)

num_type = read_tsv(dat_file, show_col_types=FALSE) %>% select(comm_num, comm_type) %>% unique() %>% 
	mutate(comm_num = as.character(comm_num))

comm_prev = comm_prev %>% left_join(num_type) %>% 
	mutate(
		comm_type = if_else(comm_num == 'fishing', 'Fishing', comm_type),
		comm_type = if_else(comm_num == 'trading', 'Trading', comm_type),
		comm_type = if_else(comm_num == 'agrarian', 'Agrarian', comm_type))

v_prev = read_tsv(comm_v_prev_file, show_col_types=FALSE) %>% filter(round == 18)

class_prev = read_tsv(class_prev_pred_file, show_col_types=FALSE)

max_y = ceiling(100*max(comm_prev$upr))
max_val = max_y
#### PANEL A ####
# nnrti dr prevalence by community
nnrti_comm_type_prev = comm_prev %>% filter(class == 'nnrti' & type == 'pool') %>%
	arrange(!(type == 'pool'), -fit) %>%
	mutate(comm_num = ordered(comm_num, levels=comm_num))

nnrti_comm_prev = comm_prev %>% filter(class == 'nnrti' & type != 'pool') %>%
	arrange(!(type == 'pool'), -fit) %>%
	mutate(comm_num = ordered(comm_num, levels=comm_num))

a1 = ggplot(nnrti_comm_type_prev, aes(x=comm_num, y=fit*100, fill=comm_type)) + 
	geom_errorbar(aes(ymin=lwr*100, ymax=upr*100), color='#333333', linewidth=0.25, width=0) +
	geom_point(shape=shapes['nnrti'], color='#333333', size=4) +
	xlab('community\ntype\n') +
	ylab('population prevalence of\nviremic resistance (%)') +	
	scale_fill_manual(values = colors,  name='', guide="none") +
	ylim(0,max_y) +
	gtheme + 
	theme(
		plot.margin = margin(0, 0, 0, 0, "pt"),
			plot.title=element_text(size=20),
			axis.text.x = element_text(angle=90, size=10, vjust = 0.5, hjust=1),
			legend.justification = c(1, 1), legend.position = c(1, 1),
					panel.grid.minor.x = element_blank(),
				panel.grid.major.x = element_blank())

a2 = ggplot(nnrti_comm_prev, aes(x=comm_num, y=fit*100, fill=comm_type)) + 
	geom_errorbar(aes(ymin=lwr*100, ymax=upr*100), color='#333333', linewidth=0.25, width=0) +
	geom_point(shape=shapes['nnrti'], color='#333333', size=4) +
	scale_x_discrete(breaks=NULL, labels=NULL) +
	scale_fill_manual(values = colors, name='') +
	xlab('community\n') +
	ylim(0,max_y) +
	ylab(NULL)+
	gtheme + 
	theme(
			axis.ticks.y = element_blank(),
			plot.margin = margin(0, 0, 0, 0, "cm"),
			plot.title=element_text(size=20),
			axis.text.x = element_text(angle=90, size=10, vjust = 0.5, hjust=1),
			axis.text.y = element_blank(),
			legend.justification = c(1, 1), legend.position = c(1, 1),
					panel.grid.minor.x = element_blank(),
				panel.grid.major.x = element_blank())

a = a1 + a2 + plot_layout(widths = c(1,  5)) + 
	plot_annotation(title='NNRTI', 
		theme = theme(plot.title = element_text(hjust = 0.5)))

#### PANEL B ####
nrti_comm_type_prev = comm_prev %>% filter(class == 'nrti' & type == 'pool') %>%
	arrange(!(type == 'pool'), -fit) %>%
	mutate(comm_num = ordered(comm_num, levels=comm_num))

nrti_comm_prev = comm_prev %>% filter(class == 'nrti' & type != 'pool') %>%
	arrange(!(type == 'pool'), -fit) %>%
	mutate(comm_num = ordered(comm_num, levels=comm_num))

b1 = ggplot(nrti_comm_type_prev, aes(x=comm_num, y=fit*100, fill=comm_type)) + 
	geom_errorbar(aes(ymin=lwr*100, ymax=upr*100), color='#333333', linewidth=0.25, width=0) +
	geom_point(shape=shapes['nrti'], color='#333333', size=4) +
	xlab('community\ntype\n') +
	ylab('') +	
	scale_fill_manual(values = colors,  name='', guide="none") +
	ylim(0,max_y) +
	gtheme + 
	theme(
		plot.margin = margin(0, 0, 0, 0, "pt"),
			plot.title=element_text(size=20),
			axis.text.x = element_text(angle=90, size=10, vjust = 0.5, hjust=1),
			legend.justification = c(1, 1), legend.position = c(1, 1),
					panel.grid.minor.x = element_blank(),
				panel.grid.major.x = element_blank())

b2 = ggplot(nrti_comm_prev, aes(x=comm_num, y=fit*100, fill=comm_type)) + 
	geom_errorbar(aes(ymin=lwr*100, ymax=upr*100), color='#333333', linewidth=0.25, width=0) +
	geom_point(shape=shapes['nrti'], color='#333333', size=4) +
	scale_x_discrete(breaks=NULL, labels=NULL) +
	scale_fill_manual(values = colors, name='', guide="none") +
	xlab('community\n') +
	ylim(0,max_y) +
	ylab(NULL)+
	gtheme + 
	theme(
			axis.ticks.y = element_blank(),
			plot.margin = margin(0, 0, 0, 0, "cm"),
			plot.title=element_text(size=20),
			axis.text.x = element_text(angle=90, size=10, vjust = 0.5, hjust=1),
			axis.text.y = element_blank(),
			legend.justification = c(1, 1), legend.position = c(1, 1),
					panel.grid.minor.x = element_blank(),
				panel.grid.major.x = element_blank())

b = b1 + b2 + plot_layout(widths = c(1,  5)) + 
	plot_annotation(title='NRTI', 
		theme = theme(plot.title = element_text(hjust = 0.5)))

#### PANEL C ####
# correlation between viremic HIV prevalence and resistant viremic HIV prevalence
nnrti_v_comm_prev = v_prev %>% 
	mutate(comm_num = case_when(
		comm_type == 'Agrarian' ~ 'agrarian',
		comm_type == 'Fishing' ~ 'fishing',
		comm_type == 'Trading' ~ 'trading')) %>%
	left_join(nnrti_comm_type_prev, suffix=c('_viremia', '_nnrti'), by=c('round', 'comm_num', 'comm_type')) 

# max val is in units of %
pred_vals = 
		tibble(fit_viremia=seq(0, max_val/100, .0001)) %>%
		mutate(
			fit_nnrti = 
				fit_viremia*exp((class_prev %>% filter(class == 'nnrti' & round == 18))$fit),
			lwr_nnrti = 
				fit_viremia*
					exp((class_prev %>% filter(class == 'nnrti' & round == 18))$fit + 
						qnorm(0.05/2) * (class_prev %>% filter(class == 'nnrti' & round == 18))$se.fit),
			upr_nnrti = 
				fit_viremia*
					exp((class_prev %>% filter(class == 'nnrti' & round == 18))$fit - 
						qnorm(0.05/2) * (class_prev %>% filter(class == 'nnrti' & round == 18))$se.fit))


c = ggplot(nnrti_v_comm_prev, aes(x=fit_viremia*100, y=fit_nnrti*100)) + 
		geom_ribbon(data=pred_vals, 
			aes(ymin=lwr_nnrti*100, ymax=upr_nnrti*100), color='#999999', fill='#d6d6d6') +
		geom_line(data=pred_vals, color='#333333') +
		#geom_abline(intercept = 0, color='#adadad', linetype="dashed") +
		geom_errorbar(aes(ymin=lwr_nnrti*100, ymax=upr_nnrti*100), color='#333333', linewidth=0.25, width=0) +
		geom_errorbar(aes(xmin=lwr_viremia*100, xmax=upr_viremia*100), color='#333333', linewidth=0.25, width=0) +
		geom_point(aes(fill=comm_type), shape=shapes['nnrti'], color='#333333', size=4) +
		xlab('population prevalence\nof viremia (%)') +
		ylab('population prevalence\nof resistant viremia (%)') +
		scale_fill_manual(values = colors, guide='none') +
		ylim(c(0,max_val)) +
		xlim(c(0,max_val)) +
		ggtitle('NNRTI') +
		coord_fixed() +
		gtheme+
		theme(plot.margin = margin(0, 0, 0, 0, "cm"))


#### PANEL D ####
# correlation between viremic HIV prevalence and resistant viremic HIV prevalence
nrti_v_comm_prev = v_prev %>% 
		mutate(comm_num = case_when(
			comm_type == 'Agrarian' ~ 'agrarian',
			comm_type == 'Fishing' ~ 'fishing',
			comm_type == 'Trading' ~ 'trading')) %>%
	left_join(nrti_comm_type_prev, suffix=c('_viremia', '_nrti'), by=c('round', 'comm_num', 'comm_type'))

# max val is in units of %
pred_vals = 
	tibble(fit_viremia=seq(0, max_val/100, .0001)) %>%
		mutate(
			fit_nrti = 
				fit_viremia*exp((class_prev %>% filter(class == 'nrti' & round == 18))$fit),
			lwr_nrti = 
				fit_viremia*
					exp((class_prev %>% filter(class == 'nrti' & round == 18))$fit + 
						qnorm(0.05/2) * (class_prev %>% filter(class == 'nrti' & round == 18))$se.fit),
			upr_nrti = 
				fit_viremia*
					exp((class_prev %>% filter(class == 'nrti' & round == 18))$fit - 
						qnorm(0.05/2) * (class_prev %>% filter(class == 'nrti' & round == 18))$se.fit))


d = ggplot(nrti_v_comm_prev, aes(x=fit_viremia*100, y=fit_nrti*100)) + 
		geom_ribbon(data=pred_vals, 
			aes(ymin=lwr_nrti*100, ymax=upr_nrti*100), color='#999999', fill='#d6d6d6') +
		geom_line(data=pred_vals, color='#333333') +
		#geom_abline(intercept = 0, color='#adadad', linetype="dashed") +
		geom_errorbar(aes(ymin=lwr_nrti*100, ymax=upr_nrti*100), 
			color='#333333', linewidth=0.25, width=0) +
		geom_errorbar(aes(xmin=lwr_viremia*100, xmax=upr_viremia*100), 
			color='#333333', linewidth=0.25, width=0) +
		geom_point(shape=shapes['nrti'], color='#333333', size=4, aes(fill=comm_type)) +
		xlab('popuation prevalence\nof viremia (%)') +
		ylab('\n') +
		scale_fill_manual(values = colors, guide='none') +
		ylim(c(0,max_val)) +
		xlim(c(0,max_val)) +
		coord_fixed() +
		ggtitle('NRTI') +
		gtheme +
		theme(plot.margin = margin(0, 0, 0, 0, "cm"))

top = plot_grid(a,b, nrow=1, labels=c('a', 'b'))
#bottom = plot_spacer() + c + d + plot_spacer()
bottom = plot_grid(NULL,c,d,NULL, nrow=1,labels=c('','c','d',''), rel_widths=c(0.5,1,1,0.5))
p = plot_grid(top, bottom,nrow=2)
ggsave('figures/final/figure_3.pdf', p, width=15, height=9.6,  bg = "transparent")


