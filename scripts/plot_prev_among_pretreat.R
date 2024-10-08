suppressMessages(library(tidyverse))
suppressMessages(library(cowplot))
suppressMessages(source('scripts/utils.R'))

dat_file = 'data/rakai_drug_resistance_categorized.tsv'
mut_class_file = 'data/mut_class.tsv'
hiv_dr_cat = read_tsv(dat_file, show_col_types=FALSE)
round_dates = hiv_dr_cat %>% select(round, round_mid_date) %>% unique() %>% mutate(round = as.character(round))
mut_class = read_tsv(mut_class_file, show_col_types=FALSE)


#### PANEL A ####
# prevalence of pre-treatment resistance over time 
dodge = c('nnrti'=25, "nrti" = 0, "pi" = -25)
order = c('nnrti' = 3, "nrti" = 2, "pi" = 1)
dr_amongPretreat_file = 'models/dr_amongPretreat_prev_pred.tsv'
dr_amongPretreat= read_tsv(dr_amongPretreat_file, show_col_types=FALSE) %>%
	mutate(round = as.character(round)) %>%
	left_join(round_dates, by="round") %>%
	mutate(round_mid_date_dodge = round_mid_date + dodge[class]) %>%
	arrange(by=order[class])


a = ggplot(dr_amongPretreat, 
			aes(x=round_mid_date_dodge, y=fit*100, color=class, shape=class, fill=class, group=class)) +
		geom_line(color='#333333') +
		geom_errorbar(aes(ymin=lwr*100, ymax=upr*100), color='#333333', linewidth=1, width=0) +
		geom_errorbar(aes(ymin=lwr*100, ymax=upr*100, color=class), width=0) +
		geom_point(color='#333333', size=4) +
		scale_shape_manual(values=shapes, name='') +
		scale_fill_manual(values=colors, name='') +
		scale_color_manual(values=colors) +
		guides(color = "none",
			shape = guide_legend(position = "inside"),
			fill = guide_legend(position = "inside"), override.aes = list(size = 5)) +
		xlab('date') + 
		ylab('\nprevalence among pre-treatment\nviremic PLHIV (%)') +
		gtheme +
		theme(legend.justification.inside = c(0, 1),
			legend.text=element_text(size=16),
			plot.margin = margin(0.3,0.3,0.3,0.3, "cm"))


#### PANEL C ####
# risk ratio of 15 most frequent mutations in pretreated data compared to R15
#pretreat_dr_muts_prev_file = 'models/dr_muts_amongPretreat_prev_pred.tsv'
#pretreat_dr_muts_prev = read_tsv(pretreat_dr_muts_prev_file, show_col_types=FALSE) %>%
#	rename(mut = class)
#plot_muts = pretreat_dr_muts_prev %>% filter(round == 18) %>% arrange(-fit) %>%
#	slice(1:10) %>% select(mut) %>% mutate(idx = seq(1,n()))

#pretreat_dr_muts_rr_file = 'models/dr_muts_amongPretreat_rr_pred.tsv'
#pretreat_dr_muts_rr = read_tsv(pretreat_dr_muts_rr_file, show_col_types=FALSE) %>%
#	rename(mut = class)
#pretreat_dr_muts_rr = pretreat_dr_muts_rr %>% inner_join(plot_muts, by="mut") %>%
#	left_join(mut_class, by='mut') %>%
#	filter(var == 'round18') %>%
#	arrange(-RR)


#y_breaks = 2^seq(floor(log2(min(pretreat_dr_muts_rr$LCI))), 
#	ceiling(log2(max(pretreat_dr_muts_rr$UCI))),1)
#y_breaks = y_breaks[2:length(y_breaks)]


#c = ggplot(pretreat_dr_muts_rr, aes(x=idx, y=RR, fill=class, color=class, shape=class)) +
#	geom_hline(yintercept=1, linetype="dashed", 
#                color = "#333333", linewidth=1) +
#	geom_errorbar(aes(ymin=LCI, ymax=UCI), color='#333333', linewidth=1, width=0) +
#	geom_errorbar(aes(ymin=LCI, ymax=UCI), width=0) +
#	geom_point(color='#333333', size=4) +
#	scale_x_continuous(breaks = pretreat_dr_muts_rr$idx,labels = pretreat_dr_muts_rr$mut, 
#		expand = c(0.04,0.04)) +
#	scale_y_continuous(trans='log2', breaks=y_breaks) +
#	scale_fill_manual(values=colors, name='') +
#	scale_color_manual(values=colors, guide=NULL) +
#	scale_shape_manual(values=shapes, name='') +
#	xlab('substitution') +
#	ylab('\nprev. ratio among viremic\npre-treatment PLHIV (v. R15)') +
#	ggtitle('2017') +
#	guides(
		#fill = guide_legend(position = "inside"), 
		#shape = guide_legend(position = "inside"), 
#		override.aes = list(size = 5)) +
#	gtheme +
#	theme(
#		axis.title.y = element_text(size=12),
#		axis.text.x = element_text(angle=90, size=10, vjust = 0.5, hjust=1),
#		panel.grid.minor.x = element_blank(),
#		panel.grid.major.x = element_blank(),
#		legend.justification.inside = c(0, 1))


#### PANEL B ####
# prevalence of 10 most frequent mutations
pretreat_dr_muts_prev_file = 'models/dr_muts_amongPretreat_prev_pred.tsv'
pretreat_dr_muts_prev = read_tsv(pretreat_dr_muts_prev_file, show_col_types=FALSE) %>%
	rename(mut = class)
plot_muts = pretreat_dr_muts_prev %>% filter(round == 18) %>% arrange(-fit) %>%
	slice(1:10) %>% select(mut) %>% mutate(idx = seq(1,n()))

pretreat_dr_muts_prev = pretreat_dr_muts_prev %>% inner_join(plot_muts) %>%
	filter(round == 18) %>%
	left_join(mut_class, by='mut')

b = ggplot(pretreat_dr_muts_prev, aes(x=idx, y=fit*100, fill=class, color=class, shape=class)) +
	geom_errorbar(aes(ymin=lwr*100, ymax=upr*100), color='#333333', linewidth=1, width=0) +
	geom_errorbar(aes(ymin=lwr*100, ymax=upr*100), width=0) +
	geom_point(color='#333333', size=4) +
	scale_x_continuous(breaks = pretreat_dr_muts_prev$idx,labels = pretreat_dr_muts_prev$mut, expand = c(0.04,0.04)) +
	scale_fill_manual(values=colors, name=NULL) +
	scale_color_manual(values=colors, guide='none') +
	scale_shape_manual(values=shapes, name=NULL) +
	ggtitle(2017) +
	xlab('substitution') +
	ylab('\nprevalence among pre-treatment\nviremic PLHIV (%)') +
	gtheme +
	theme(
		#axis.title.y = element_text(size=12),
		axis.text.x = element_text(angle=90, size=10, vjust = 0.5, hjust=1),
		panel.grid.minor.x = element_blank(),
		panel.grid.major.x = element_blank())

#p = plot_grid(
#	plot_grid(NULL, a, NULL, ncol=1, rel_heights=c(0.2, 0.6, 0.2), labels=c('', 'A', '')), 
#	plot_grid(d,c, nrow=2, labels=c('B', 'C')), nrow=1, rel_widths=c(0.5, 0.4))
#p = plot_grid(a,c,b,d, nrow=2, rel_widths=c(0.6, 0.4))
p = plot_grid(a,b,nrow=1,labels=c('A', 'B'), rel_widths=c(0.6, 0.4))

ggsave('figures/prev_among_pretreat.pdf', p, width=12, height=4.8)