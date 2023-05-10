#!/usr/bin/env bash

find . -iname "sim_asr_*" -print0 | xargs -0 -P7 -I% ./plot_output.r %
