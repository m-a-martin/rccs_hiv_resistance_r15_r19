suppressMessages(require(tidyverse))
suppressMessages(require(cowplot))
suppressMessages(require(patchwork))
suppressMessages(require(ggpattern))
suppressMessages(require(scatterpie))
suppressMessages(source('scripts/utils.R'))


prev_dat_file = 'models/multi_amongPar_prev_pred.tsv'
dat_file = 'data/rakai_drug_resistance_categorized.tsv'


prev_dat = read_tsv(prev_dat_file, show_col_types=FALSE) %>%
	filter(round == 18)

prev_dat = prev_dat %>% arrange(-fit) %>%
	mutate(x = seq(1,n())) %>%
	separate(class, c('A', 'B', 'C'), remove=FALSE) %>%
	mutate(
		n = (1*!is.na(A)) + (1*!is.na(B)) + (1*!is.na(C)),
		A = replace_na(A, '-'),
		B = replace_na(B, '-'),
		C = replace_na(C, '-'),
		nnrti = 1*(A == 'nnrti' | B == 'nnrti' | C == 'nnrti')/n,
		nrti = 1*(A == 'nrti' | B == 'nrti' | C == 'nrti')/n,
		pi = 1*(A == 'pi' | B == 'pi' | C == 'pi')/n)

cat_dat = prev_dat %>% select(x, class) %>%
	separate(class, c("A", "B", "C")) %>%
	pivot_longer(-x) %>%
	filter(!is.na(value)) %>%
	mutate(value = ordered(value, levels=rev(unique(value))))

x_scale = 7.5
x_span = range(prev_dat$x/x_scale)
x_lims = c(x_span[1] - (x_span[2] - x_span[1])*0.05, 
	x_span[2] + (x_span[2] - x_span[1])*0.05)

p1 = 
	ggplot() + 
		geom_errorbar(data=prev_dat %>% filter(y > 0), aes(x = x/x_scale, ymin=lwr*100, ymax=upr*100), color='#333333', linewidth=1, width=0) +
		geom_scatterpie(data=prev_dat %>% mutate(x = x/x_scale), aes(x=x, y=fit*100, group=class), 
			pie_scale=1.25, cols=c('nnrti', 'nrti', 'pi'), color='#333333',  show.legend = F) + 
		geom_point(data = prev_dat %>% filter(!grepl('_', class)), aes(x=x/x_scale, y=fit*100, fill=class), 
			shape=21,size=6.5, color='#333333') +  
		coord_equal() +
		ylab('population prevalence (%)') +
		ggtitle('viremic resistance among\nRCCS participants (2017)') +
		scale_fill_manual(values = colors, name=NULL) +
		guides(fill = guide_legend(position = "inside")) +
		gtheme + 
		scale_x_continuous(name=NULL, labels=NULL, breaks=prev_dat$x/x_scale, limits = x_lims) + 
		theme(
			legend.justification.inside = c(1, 1),
			panel.grid.major.x = element_blank(), 
			panel.grid.minor.x = element_blank())


bg_points = tibble(
	x = rep(unique(cat_dat$x), each=length(unique(cat_dat$value))),
	value = rep(unique(cat_dat$value), length(unique(cat_dat$x))))
#bg_points = tibble(
#	x=rep(unique(cat_dat$x),each=length(unique(cat_dat$value))),
#	value = rep(unique(cat_dat$value), length(unique(cat_dat$value))*length(unique(cat_dat$x))))


p2 = ggplot(cat_dat, aes(x=x/x_scale, y=value, group=x)) +
	geom_point(data=bg_points, size=5, color='#eaeaea') + geom_line(color='#333333') +
	geom_point(size=4, color='#333333') + 
	geom_line(color='#333333', linewidth=1.25) +
	scale_x_continuous(breaks=NULL, name=NULL,  limits = x_lims) +
	ylab(NULL) +
	gtheme +
	theme(
		axis.ticks.y = element_blank(),
		panel.border = element_blank(),
		panel.grid.major.x = element_blank(),
		panel.grid.minor = element_blank()) 

hiv_dr_cat = read_tsv(dat_file, show_col_types=FALSE) %>% filter(round == 18)

pretreat_dat =  hiv_dr_cat %>% 
	filter(pre_treatment & !is.na(nnrti) & !is.na(nrti) & !is.na(pi)) %>%
	mutate(
		nnrti_r = (nnrti == 'intermediate') | (nnrti == 'high'),
		nrti_r = (nrti == 'intermediate') | (nrti == 'high'),
		pi_r = (pi == 'intermediate') | (pi == 'high'),
		r = case_when(
			nnrti_r & nrti_r & pi_r ~ 'NNRTI, NRTI, PI',
			nnrti_r & nrti_r & !pi_r ~ 'NNRTI, NRTI',
			nnrti_r & !nrti_r & pi_r ~ 'NNRTI, PI',
			nnrti_r & !nrti_r & !pi_r ~ 'NNRTI',
			!nnrti_r & nrti_r & pi_r ~ 'NRTI, PI',
			!nnrti_r & !nrti_r & pi_r ~ 'PI',
			!nnrti_r & nrti_r & !pi_r ~ 'NRTI')) %>% 
	filter(!is.na(r)) %>%
	group_by(r) %>%
	summarise(n=n()) %>%
	arrange(-n) %>%
	mutate(idx = seq(1,n()),
		class_1 = str_split(r, ',', simplify=T)[,1],
		class_2 = str_split(r, ', ', simplify=T)[,2],
		class_3 = str_split(r, ', ', simplify=T)[,3])

treat_dat =  hiv_dr_cat %>% 
	filter(!pre_treatment & !is.na(nnrti) & !is.na(nrti) & !is.na(pi)) %>%
	mutate(
		nnrti_r = (nnrti == 'intermediate') | (nnrti == 'high'),
		nrti_r = (nrti == 'intermediate') | (nrti == 'high'),
		pi_r = (pi == 'intermediate') | (pi == 'high'),
		r = case_when(
			nnrti_r & nrti_r & pi_r ~ 'NNRTI, NRTI, PI',
			nnrti_r & nrti_r & !pi_r ~ 'NNRTI, NRTI',
			nnrti_r & !nrti_r & pi_r ~ 'NNRTI, PI',
			nnrti_r & !nrti_r & !pi_r ~ 'NNRTI',
			!nnrti_r & nrti_r & pi_r ~ 'NRTI, PI',
			!nnrti_r & !nrti_r & pi_r ~ 'PI',
			!nnrti_r & nrti_r & !pi_r ~ 'NRTI')) %>% 
	filter(!is.na(r)) %>%
	group_by(r) %>%
	summarise(n=n()) %>%
	arrange(-n) %>%
	mutate(idx = seq(1,n()),
		class_1 = str_split(r, ',', simplify=T)[,1],
		class_2 = str_split(r, ', ', simplify=T)[,2],
		class_3 = str_split(r, ', ', simplify=T)[,2])

max_y = 10*ceiling(max(c(treat_dat$n, pretreat_dat$n))/10)

a = ggplot(pretreat_dat, aes(x=idx, y = n, group=r, pattern=r, fill=class_1, pattern_fill=class_2, label = n)) +
	geom_col_pattern(pattern_spacing=0.025, pattern_density=0.5, color='#333333') +
	geom_col_pattern(data = pretreat_dat %>% filter(r == 'NNRTI, NRTI, PI'),
		pattern_angle=-45, pattern_spacing=0.025, pattern_density=0.5, color='#333333',
		pattern_fill=colors['pi'], fill='white', alpha=0)+
	geom_text(color='#333333', vjust = -max_y*0.005) +
	scale_x_continuous(breaks=pretreat_dat$idx, labels=pretreat_dat$r) +
	scale_pattern_manual(values=patterns, guide='none') +
	scale_fill_manual(values=colors, guide='none') +
	scale_pattern_fill_manual(values=colors, guide='none') +
	gtheme +
	xlab('resistance class') +
	ylab('participants') +
	ylim(c(0, max_y)) +
	ggtitle('pre-treatment viremic\nresistance (2017)') +
	theme(
		panel.grid.major.x = element_blank(),
		  panel.grid.minor.x = element_blank(),
		  axis.text.x = element_text(angle=90, size=10, vjust = 0.5, hjust=1))

b = ggplot(treat_dat, aes(x=idx, y = n, group=r, pattern=r, fill=class_1, pattern_fill=class_2, label=n)) +
	geom_col_pattern(pattern_spacing=0.025, pattern_density=0.5, color='#333333') +
	geom_col_pattern(data = pretreat_dat %>% filter(r == 'NNRTI, NRTI, PI'),
		pattern_angle=-45, pattern_spacing=0.025, pattern_density=0.5, color='#333333',
		pattern_fill=colors['pi'], fill='white', alpha=0)+
	geom_text(color='#333333', vjust = -max_y*0.005) +
	scale_x_continuous(breaks=treat_dat$idx, labels=treat_dat$r) +
	scale_pattern_manual(values=patterns, guide='none') +
	scale_fill_manual(values=colors, guide='none') +
	scale_pattern_fill_manual(values=colors, guide='none') +
	gtheme +
	xlab('resistance class') +
	ylab('participants') +
	ylim(c(0, max_y)) +
	ggtitle('treatment-experienced viremic\nresistance (2017)') +
	theme(
		panel.grid.major.x = element_blank(),
		  panel.grid.minor.x = element_blank(), 
		  axis.text.x = element_text(angle=90, size=10, vjust = 0.5, hjust=1))


#p = p1/(plot_spacer() + p2 + plot_spacer() +
#	plot_layout(widths=c(0.16, 1,0.24))) +  plot_layout(heights = c(1, 0.25))
#ggsave('test.pdf', p, width=8.75, height=4.8)

p = plot_grid(
	p1/(plot_spacer() + p2 + plot_spacer() +
		plot_layout(widths=c(-0.035, 1, 0.05))) +  plot_layout(heights = c(1, 0.25)),
	a,b, nrow=1, rel_widths=c(0.5, 0.3, 0.315),  labels = c('A', 'B', 'C')) 
ggsave('figures/multiclass_prev_among_par.pdf', p, width=14, height=4.8)
