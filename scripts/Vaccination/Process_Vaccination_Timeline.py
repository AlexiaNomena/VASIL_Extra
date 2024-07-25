#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 30 12:02:05 2023

@author: raharinirina
"""

import sys
import pandas as pd
import pdb
import numpy as np
import pickle

try:
    df_Vaccines = pd.read_csv(sys.argv[1])
    print(df_Vaccines["date"])
except:
    df_Vaccines = pd.read_csv(sys.argv[1], sep="\t")

date_start = str(sys.argv[2])
date_end = str(sys.argv[3])
switch_vacc = str(sys.argv[4])
vacc_considered = str(sys.argv[5]).split("/")

save_mut_to = str(sys.argv[-2])
save_to = str(sys.argv[-1])

### Drop dates between 
if date_start in df_Vaccines["date"].tolist():
    df_Vaccines.drop(index = df_Vaccines.index[:list(df_Vaccines["date"]).index(date_start)], inplace = True)

if date_end in df_Vaccines["date"].to_list():
    df_Vaccines.drop(index = df_Vaccines.index[list(df_Vaccines["date"]).index(date_end)+1:], inplace = True)


dates_vacc = df_Vaccines["date"].to_list()
### Extract relevant vaccine columns
Columns = df_Vaccines.columns

Timeline = {}
col_kept = []
for col in Columns:
    if col in vacc_considered:
    	col_kept.append(col)
    	col_new = col
    	vacc = df_Vaccines[col].to_numpy()
    	vacc[np.isnan(vacc)] = 0
    	vacc = np.diff(vacc) ### the vaccination data look like cumulative sums but we need the absolute daily numbers of vaccinated
    	vacc[vacc<0] = 0
    	if col_new == "people_fully_vaccinated":
    		Timeline["Wuhan-Hu-1*_as_"+col_new] = vacc ### always has to have the structure lineage*_as_...where lineage must be a variant within the spikegroups (otherwise it will not compute the cross-reactivity because of lack of mutation profile)
    	elif col_new == "total_boosters":
            Timeline["Wuhan-Hu-1*_as_"+col_new] = vacc ### here I used Wuhan because BA.5 does not exist in the UK timeline and it's tricky to implement another mutation profile
        #elif col_new == "another" #can add another one but must always be given as a parameter separated by / in parameter vacc_considered in vacc_config.yaml
                                  
print("Vaccines condidered:", col_kept)
Timeline["date"] = dates_vacc[:-1]
df_Timeline = pd.DataFrame(Timeline)
df_Timeline.to_csv(save_to+"/Vaccination_Timeline.csv")

Total_vacc = {}
Total_vacc["date"] = dates_vacc[:-1]
Total_vacc["vaccinated"] = np.sum(df_Timeline.to_numpy()[:, np.array(df_Timeline.columns) != "date"], axis = 1)

df_total = pd.DataFrame(Total_vacc)    
df_total.to_csv(save_to+"/Vaccination_Total.csv")