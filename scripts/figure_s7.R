library(tidyverse)
source('scripts/utils.R')
library(cowplot)

t97a_file = 'models/pretreat_int97a_amongPLWHIV_prev_pred.tsv'
t97a = read_tsv(t97a_file)

dat_file = 'data/rakai_drug_resistance_categorized.tsv'
mut_file = 'data/rakai_drug_resistance_mut_formatted.tsv'

hiv_mut = read_tsv(mut_file)
hiv_dr_cat = read_tsv(dat_file)


# get samples with inT97A
# and merge with all data
hiv_dr_cat = hiv_dr_cat %>% left_join(hiv_mut %>% group_by(study_id, round) %>% summarise(t97a = 'inT97A' %in% mut) %>%
	filter(t97a)) %>%
	mutate(t97a = replace_na(t97a, FALSE))

# now filter
hiv_dr_cat = hiv_dr_cat %>%
	filter(
		!is.na(nnrti) & !is.na(nrti) & !is.na(pi) & !is.na(insti) &
		round > 14 &
		pre_treatment &
		viremic & 
		finalhiv == 'P' & 
		age_cat != '(50, 100]')


plotlist = list()
for (r_col in c('ref_subtype', 'genome_subtype', 'pol_subtype')){
	r_col_dat = hiv_dr_cat %>%
		rename(r = r_col) %>%
		group_by(r) %>%
		summarise(`ref samples` = sum(!t97a), `in97A samples` = sum(t97a)) %>%
		mutate(r = if_else(`ref samples` < 50, 'other', r))
	if (r_col == 'ref_subtype'){
		r_col_dat = r_col_dat %>%
		mutate(r = if_else(grepl('_', r), str_split(r, '_', simplify=T)[,2], r))
	}
	# todo do this in pipe?
	r_col_dat = r_col_dat %>%
		group_by(r) %>%
		summarise(
			`ref samples` = sum(`ref samples`), 
			`in97A samples` = sum(`in97A samples`), .groups='drop') %>%
		filter(r != "-") %>%
		mutate(
			`ref samples` = `ref samples`/sum(`ref samples`),
			`in97A samples` = `in97A samples`/sum(`in97A samples`))
	r_col_dat = r_col_dat %>% 
		mutate(r = ordered(r, levels=(r_col_dat %>% arrange(`ref samples`, `in97A samples`))$r)) %>%
		arrange(-`ref samples`, -`in97A samples`) %>%
		pivot_longer(-r) %>%
		mutate(
			x1 = if_else(name == 'ref samples', 1, 2),
			x2 = if_else(name == 'ref samples', 1.25, 1.75))
	plotlist[[length(plotlist) + 1]] = ggplot(r_col_dat, aes(x=x1, y=value, fill=r, group=r)) +
		geom_col(width=0.5, color='#333333') +
		geom_line(aes(x=x2), position = "stack", linewidth = 0.5, linetype='dashed') +
		scale_x_continuous(
			breaks=(r_col_dat %>% select(name, x1) %>% unique())$x1,
			labels=(r_col_dat %>% select(name, x1) %>% unique())$name) +
		scale_fill_manual(values=colors, name=var_rename[r_col]) +
		xlab('type') +
		ylab('proportion of samples') +
		gtheme +
		theme(
			panel.grid.major.x = element_blank(),
		  panel.grid.minor.x = element_blank())
	if (r_col != 'ref_subtype'){
		plotlist[[length(plotlist)]] = plotlist[[length(plotlist)]] + ylab('')
	}
}



p = plot_grid(plotlist=plotlist, labels=c('a', 'b', 'c'), nrow=1)


ggsave('figures/final/figure_s7.pdf', p, width=19, height=6.2)
