#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 23 19:10:08 2023

@author: raharinirina
"""
import os
import sys
import pandas as pd

stat = {}
try:
    os.remove("results/stat_prev.csv")
    stat["removed_prev"] = [True]
except:
    stat["removed_prev"] = ["Placeholder"]

try:
    os.remove(str(sys.argv[1])+str(sys.argv[2]))
    stat["removed_mean"] = [True]
except:
    stat["removed_mean"] = ["Placeholder"]
    
pd.DataFrame(stat).to_csv(sys.argv[3])

