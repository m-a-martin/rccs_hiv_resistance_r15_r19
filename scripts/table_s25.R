library(tidyverse)
source('scripts/utils.R')


pred_file = 'models/pretreat_amongInc_pred.tsv'
rr_file = 'models/pretreat_amongInc_rr.tsv'

ts = class_round_prev_rr_table(pred_file, rr_file)

write_tsv(ts, 'tables/table_s25.tsv', na='')

read_tsv(pred_file)