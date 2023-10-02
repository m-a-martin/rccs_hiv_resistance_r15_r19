library(tidyverse)
library(haven)
library(zoo)
source('scripts/utils.R')

hiv_dr_cat = read_tsv('data/rakai_drug_resistance_categorized.tsv')


n_per_month = hiv_dr_cat %>% 
	mutate(
		round = paste('Round ', round, sep=''),
		year = format(visit_dt,"%Y"),
		month = format(visit_dt, "%m"),
		label = if_else(month == "01", format(visit_dt, "%b\nY"), format(visit_dt, "%b")),
		yearmon = as.Date(paste(paste(year, month, sep='-'), '-01', sep=''), format='%Y-%m-%d')) %>% 
	group_by(round, yearmon, label) %>%
	summarise(n=n(), .groups="drop")

mylabels <- function(breaks){
    labels = vector("character", length(breaks))
    where_na = is.na(breaks)
    where_jan = format(breaks, "%m") == "01"
    labels[where_jan & !where_na] = format(breaks[where_jan & !where_na], "%b\n%Y")
    labels[!where_jan & !where_na] = format(breaks[!where_jan & !where_na], "%b")
    labels[where_na] = NA
    return(labels)
}


breaks = n_per_month$yearmon

g = ggplot(n_per_month, aes(x=yearmon, y=n)) +
	geom_bar(stat="identity", color='#333333', fill='#eaeaea') +
	facet_wrap(~round, ncol=1) +
	ylab('number of participates')+
	xlab('date') +
	scale_x_date(breaks="3 month", labels=mylabels) +
	gtheme

ggsave('figures/final/figure_s1.pdf', g, width=16, height=10)
