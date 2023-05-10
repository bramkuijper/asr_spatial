#!/usr/bin/env python3

import datetime as dt

exe_name = "./asr_spatial.exe"

# make the base name based on current date / time
basename = "sim_asr_"
current_date = dt.datetime.now()
basename += current_date.strftime("%Y%m%d_%H%M%S%f")

# 1st element mu_juv_f
# 2nd element mu_juv_m
mu_juv = [ [0.001,0.03], [ 0.001,0.001], [0.03,0.001]]

sigma = 0.2

nrep = 5

counter = 0

for replicate_i in range(0,nrep):
    for mu_juv_i in mu_juv:
        print(f"{exe_name} {mu_juv_i[0]} {mu_juv_i[1]} {basename}_{counter}")

        counter += 1
