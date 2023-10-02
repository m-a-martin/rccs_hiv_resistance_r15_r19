library(tidyverse)
source('scripts/utils.R')


pred_file = 'models/treat_amongPLWHIV_prev_pred.tsv'
rr_file = 'models/treat_amongPLWHIV_rr.tsv'

ts = class_round_prev_rr_table(pred_file, rr_file)

write_tsv(ts, 'tables/table_s25.tsv', na='')
