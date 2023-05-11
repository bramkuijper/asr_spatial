#!/usr/bin/env python3

import datetime as dt

exe_name = "./asr_spatial.exe"

# make the base name based on current date / time
basename = "sim_asr_"
current_date = dt.datetime.now()
basename += current_date.strftime("%Y%m%d_%H%M%S%f")

# 1st element mu_juv_f
# 2nd element mu_juv_m
mu_juv = [ [0.001,0.001,0.001,0.001] ]

sigma = [0.0]

nrep = 5

counter = 0

envt_ch = [ 0.001,0.001]

init_prob = 0.5

def list_to_string(x):
    return [str(elmt) for elmt in x]

for replicate_i in range(0,nrep):
    for mu_juv_i in mu_juv:
        for sigma_i in sigma:
            print(f"{exe_name} " + " ".join(list_to_string(mu_juv_i)) + 
                  f" {sigma_i} {envt_ch[0]} {envt_ch[1]} {init_prob} {basename}_{counter}"  )

        counter += 1
