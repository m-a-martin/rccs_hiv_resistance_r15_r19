library(tidyverse)
source('scripts/utils.R')


all_dr_pred_file = 'models/all_amongPar_prev_pred.tsv'
all_dr_rr_file = 'models/all_amongPar_rr.tsv'

ts = class_round_prev_rr_table(all_dr_pred_file, all_dr_rr_file)


write_tsv(ts, 'tables/table_s23.tsv', na='')
