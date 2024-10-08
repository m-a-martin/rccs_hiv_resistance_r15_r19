suppressMessages(library(tidyverse))
suppressMessages(library(cowplot))
suppressMessages(source('scripts/utils.R'))

dat_file = 'data/rakai_drug_resistance_categorized.tsv'
mut_dat_file = 'data/rakai_drug_resistance_mut_formatted.tsv'
hiv_pred_file = 'models/HIV_prev_pred.tsv'
viremic_pred_file = 'models/HIV_viremic_prev_pred.tsv'
pretreat_pred_file = 'models/HIV_pretreat_viremic_prev_pred.tsv'
treat_pred_file = 'models/HIV_treat_viremic_prev_pred.tsv'
mut_class_file = 'data/mut_class.tsv'
treat_dr_muts_prev_file = 'models/dr_muts_amongTreat_prev_pred.tsv'
treat_among_PLHIV_file = 'models/dr_amongTreat_prev_pred.tsv'


hiv_dr_cat = read_tsv(dat_file)
round_dates = hiv_dr_cat %>% select(round, round_mid_date) %>% unique() %>% mutate(round = as.character(round))

dodge = c('nnrti'=25, "nrti" = 0, "pi" = -25)
order = c('nnrti' = 3, "nrti" = 2, "pi" = 1)
#### PANEL E ####
# prevalence of resistance among treatment experienced PLHIV
treat_amongPLHIV = read_tsv(treat_among_PLHIV_file, show_col_types=FALSE) %>% 
	mutate(round = as.character(round)) %>%
	left_join(round_dates) %>%
	mutate(round_mid_date_dodge = round_mid_date + dodge[class]) %>%
	arrange(by=order[class])

a = ggplot(treat_amongPLHIV, 
			aes(x=round_mid_date_dodge, y=fit*100, color=class, shape=class, fill=class, group=class)) +
		geom_line(color='#333333') +
		geom_errorbar(aes(ymin=lwr*100, ymax=upr*100), color='#333333', linewidth=1, width=0) +
		geom_errorbar(aes(ymin=lwr*100, ymax=upr*100, color=class), width=0) +
		geom_point(color='#333333', size=4) +
		scale_shape_manual(values=shapes, name=NULL) +
		scale_fill_manual(values=colors, name=NULL) +
		scale_color_manual(values=colors) +
		guides(color = "none", shape = guide_legend(position = "inside"),
			fill = guide_legend(position = "inside")) +
		xlim(c(min(treat_amongPLHIV$round_mid_date)-180, 
			max(treat_amongPLHIV$round_mid_date)+180)) +
		xlab('date') + 
		ylab('prevalence among treatment-experienced\nviremic PLHIV (%)') +
		ylim(0, max(treat_amongPLHIV$upr)*1.225*100) + 
		gtheme  +
		theme(legend.justification.inside = c(1, 1),
			axis.title.y = element_text(size=14))


#### PANEL B ####
# prevalence of 15 most frequent mutations in treated data compared to R15
treat_dr_muts_prev = read_tsv(treat_dr_muts_prev_file, show_col_types=FALSE) 	%>% rename(mut = class)

mut_class = read_tsv(mut_class_file, show_col_types=FALSE)

plot_muts = treat_dr_muts_prev %>% filter(round == 18) %>% arrange(-fit) %>%
	slice(1:10) %>%
	select(mut) %>% mutate(idx = seq(1,n()))

treat_dr_muts_prev = treat_dr_muts_prev %>% inner_join(plot_muts, by='mut') %>%
	filter(round == 18) %>%
	left_join(mut_class, by='mut')

b = ggplot(treat_dr_muts_prev, aes(x=idx, y=fit*100, fill=class, color=class, shape=class)) +
	geom_errorbar(aes(ymin=lwr*100, ymax=upr*100), color='#333333', linewidth=1, width=0) +
	geom_errorbar(aes(ymin=lwr*100, ymax=upr*100), width=0) +
	geom_point(color='#333333', size=4) +
	scale_x_continuous(breaks = treat_dr_muts_prev$idx,labels = treat_dr_muts_prev$mut, expand = c(0.04,0.04)) +
	scale_fill_manual(values=colors, name=NULL) +
	scale_color_manual(values=colors, guide='none') +
	scale_shape_manual(values=shapes, name=NULL) +
	guides(shape = guide_legend(position = "inside"),
			fill = guide_legend(position = "inside")) +
	xlab('substitution') +
	ylab('prevalence among treatment-experienced\nviremic PLHIV (%)') +
	ggtitle(2017) +
	gtheme +
	theme(
		axis.title.y = element_text(size=14),
		axis.text.x = element_text(angle=90, size=10, vjust = 0.5, hjust=1),
		panel.grid.minor.x = element_blank(),
		panel.grid.major.x = element_blank(),
		legend.justification.inside = c(1, 1))


p = plot_grid(a,b,nrow=1,labels=c('A', 'B'))

ggsave('figures/prev_among_treat.pdf', p, width=8, height=4.8)
