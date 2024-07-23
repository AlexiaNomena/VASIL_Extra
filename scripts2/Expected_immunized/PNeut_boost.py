#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  6 14:30:28 2024

@author: raharinirina
"""
import sys
import pickle
import pandas as pd
import numpy as np
import warnings
import re
import numpy.ma as ma
from functools import partial
import joblib as jb
import os
from scipy.integrate import solve_ivp
import pdb
# Disable
def blockPrint():
    sys.stdout = open(os.devnull, 'w')

# Restore
def enablePrint():
    sys.stdout = sys.__stdout__


def PNeut_boost(Cross_delta, Cross_Neut_File, dms_data, t_boost, clinical_VE, groups_list, antigen_list, result_location = ""):
    """Load delta relevant cross_neutralization files"""
    file1 = open(Cross_delta, "rb") # premade simulations
    Cross_with_delta_validation = pickle.load(file1)
    variants_x_names_show = Cross_with_delta_validation["variant_list"]
    Cross_with_delta_validation.pop("variant_list")
    file1.close()
    Ab_classes = list(Cross_with_delta_validation.keys())
    
    """Load relevant cross neutralization files"""
    file1 = open(Cross_Neut_File, "rb") # load cross reactivity if parameter is not a directory
    Cross_react_dic = pickle.load(file1)
    variants_in_cross = Cross_react_dic["variant_list"]
    Cross_react_dic.pop("variant_list")
    file1.close()
    
    """Escape fraction data """
    Escape_Fraction = pd.read_csv(dms_data)
    
    """Load vaccing efficacy data for fitting """
    Load_Delta = pd.read_excel(clinical_VE, engine='openpyxl')
    
    """Load variants list and antigen list"""
    Lin_list = groups_list
    cross_list = []
    run_k = True
    not_pres = []
    antigen_list = antigen_list
    
    def Antibody_ranges(thalf_vec, tmax_vec, t, Ab_classes):
        N = len(Ab_classes) # number of antibody classes
        is_log = False # if True, it returns the log of the antibody concentration
        dataname = "Ab_%d"%N
        solver = "lm" # root solver method for finding absorption rate ka (see scipy.optimize.root)
        c_t_vec = np.zeros((len(thalf_vec), len(tmax_vec), len(Ab_classes), len(t)))
        c_dframe_dic = {}
        for m in range(len(thalf_vec)):
            for n in range(len(tmax_vec)):
                t_half = thalf_vec[m]*np.ones(N) # antibody half-life for all antibody classes, respecively
                t_max = tmax_vec[n]*np.ones(N) # time to the peak of antibody concentration for each antibody class
                params_dic = {"t_max":t_max, "t_half":t_half}
                c_t, c_dframe, ka, ke, c_max_dic = Antibody(t = t, params_dic = params_dic, is_log = is_log, Ab_names = Ab_classes, ka_solver = solver)
                c_t_vec[m, n, :, :] = c_t
                c_dframe_dic["(%d, %d)"%(m, n)] = c_dframe
                
        return c_t_vec, c_dframe_dic, dataname
    
    def Antibody(t, params_dic, is_log = False, Ab_names = None, ka_solver = "lm"):
        """
        @brief: Compute Antibody Concentration as a function of time for N antibody classes
        
        Parameters
        ----------
        t : T time points ndarray (T) 
        params_dic : dictionary of parameters
                    params_dic["t_max"] = t_{max}, ndarray (N, )
                    params_dic["t_half"] = t_{1/2}, ndarray (N, ) 
        is_log : bool, optional
                 True if return the log of concentration. The default is False.
        ka_solver : root solver method for finding absorption rate k_a (see scipy.optimize.root), optional. The default is lm
    
        Returns
        -------
        Antibody concentration at time t.
    
        """
        t_max = params_dic["t_max"]
        t_half = params_dic["t_half"]
        
        # Antibody elimination rate
        ke = np.log(2)/t_half
        
        # Antibody absorption rate
        guess = np.ones(len(t_max))
        warnings.filterwarnings("ignore")
        ka = root(ka_solve, guess, args = (ke, t_max), method = ka_solver).x
        
        #if not np.all(np.isclose(ka_solve(ka, ke, t_max), np.zeros(len(t_max)))):
            #print("\n k_a was found correctly:", np.all(np.isclose(ka_solve(ka, ke, t_max), np.zeros(len(t_max)))), "\n", params_dic)
    
        # Compute Normalized Concentration
        c_max = (np.exp(- ke*t_max) - np.exp(- ka*t_max))
        #c_t = (np.exp(- ke[:, np.newaxis]*t) - np.exp(- ka[:, np.newaxis]*t))/c_max[:, np.newaxis]
        
        solve_cdiff = solve_ivp(Dc_boost, [min(t), max(t)], np.array([1, 0]), args = (ka[0], ke[0]), dense_output=True) ### ka and ke here are vectors giving the same number for all antibody classes (same PK for all ab)
        c_t_full = solve_cdiff.sol(t)
        c_t = c_t_full[1, :]
        
        solve_cdiff2 = solve_ivp(Dc_boost, [t_boost, max(t)], np.array([1, 0]), args = (ka[0], ke[0]), dense_output=True) ### ka and ke here are vectors giving the same number for all antibody classes (same PK for all ab)
        c_t_full2 = solve_cdiff2.sol(np.arange(t_boost, max(t)+1))
        c_t[list(t).index(t_boost):] += c_t_full2[1, :]
        
        if is_log:
            c_t = np.log(c_t)
    
        # Build pandas dataframe‚
        df = {}
        c_max_dic = {}
        df["Days"] = t
        for i in range(len(t_max)):
            if Ab_names is None:
                df["Ab class %d"%(i+1)] = c_t #[i, :] # same PK for all ab and the vector output here is already one-dimensional
                c_max_dic["Ab class %d"%(i+1)] = c_max[i]
            else:
                df[Ab_names[i]] = c_t#[i, :] # same PK for all ab and the vector output here is already one-dimensional
                c_max_dic[Ab_names[i]] = c_max[i]
            
            
        df = pd.DataFrame(df)
                 
        return c_t, df, ka, ke, c_max_dic
        
    
    def Dc_boost(t, X, ka, ke):
        dX0 = -ka*X[0]
        dX1 = ka*X[0] - ke*X[1]
        return np.array([dX0, dX1])
    
    def ka_solve(ka, ke, t_max):
        if np.all(ka)>0:
            res = np.divide(t_max*(ka - ke) - (np.log(ka) - np.log(ke)), (ka - ke), out = np.ones(len(ke)), where = (ka - ke)!=0)
        else:
            res = 1
        return res
            
    
    """Fit IC50 parameter"""
    from scipy.optimize import root
    
    """
    def vaccine_efficacy(x, ic50):
        return(x/(x + ic50))
    
    def efficacy_n_antibodies(x, ic50):
        res = 1
        for i in range(len(x)):
            ve = vaccine_efficacy(x[i], ic50[i])
            res *= (1 - ve)
        return(1 - res)
    """
    
    def sqrt_diff_FR(ic50, days_list, FR, ve_data, n, c_dframe):
        res = 0
        for d in range(len(ve_data)):
            data = ve_data[d]
            days = days_list[d]
            ve_estimate = np.zeros(len(days))
            for i in range(len(data)):
                antibody_level = c_dframe.loc[int(days[i]) - 1][1:n+1]
                #ve_estimate[i] = efficacy_n_antibodies(antibody_level, np.array(FR)*ic50)
                ve_estimate[i] = 1 - np.prod(1 - (antibody_level/(antibody_level + ic50*np.array(FR))))
            
            res += np.linalg.norm(data-ve_estimate[0:len(data)])
        return(res)
    
    
    """Compute fold change deviation for mean IC50"""
    IC50_group = Escape_Fraction.groupby('condition', as_index=False).first()[['condition', 'IC50', 'group']]
    mean_IC50_per_group = IC50_group.groupby('group')['IC50'].mean().reset_index()
    total_mean = mean_IC50_per_group['IC50'].mean()
    mean_IC50_per_group['fold_change'] = mean_IC50_per_group['IC50']/total_mean
    ntd_row = {'group': 'NTD','IC50': 1, 'fold_change': 1}
    mean_IC50_per_group = pd.concat([mean_IC50_per_group, pd.DataFrame(ntd_row, index=[10])])
    
    """Load fold change IC50 in data"""
    
    FC_ic50_dic = {Ab_classes[i]:(mean_IC50_per_group["fold_change"].values[mean_IC50_per_group["group"] == Ab_classes[i]])[0] for i in range(len(Ab_classes))}
    
    """Extract information from Clinical data"""
    def transform(x):
        x[x<0] = 0
        return x
    
    def extract_yerr(x, CI):
        lower_diff = np.minimum(transform(x), np.abs(x - CI[:, 0]))
        upper_diff = np.minimum(transform(x), np.abs(CI[:, 1] - x))
        return np.array([lower_diff, upper_diff]), CI[:, 1], CI[:, 0]
        
    
    """Fit Delta Vaccine Data """
    warnings.filterwarnings("ignore")
    days_fitting = []
    ve_fitting = []
    
    Delta_Sources = Load_Delta["Source"].values.astype(str)
    Delta_Vaccine = Load_Delta["Vaccine"].values.astype(str)
    
    All_Days_Delta = np.array([])
    All_Delta_Data = np.array([])
    All_Days_xerr_Delta = np.array([])
    All_Delta_yerr = np.array([])
    
    """ Set up Disease status to restrict the studies on any infected people"""
    keep_status_d = Load_Delta["Disease Status"].values.astype(str) == "Any infection"
    keep_method_d = (Load_Delta["Method"].values.astype(str) != "(1 - Adjusted OR)") & (Load_Delta["Method"].values.astype(str) != "(1 - OR)")
    
    u_Delta_Sources = np.unique(Delta_Sources[keep_status_d & keep_method_d])
    u_Delta_Vaccine = np.unique(Delta_Vaccine[keep_status_d & keep_method_d])
    
    u_Delta_Sources = u_Delta_Sources[~(u_Delta_Sources  == "nan")] 
    u_Delta_Vaccine = u_Delta_Vaccine[~(u_Delta_Vaccine  == "nan")]
    
    
    Delta_done = []
    for source in u_Delta_Sources:
        for vacc in u_Delta_Vaccine:
            where_source = (Delta_Sources == source)&(Delta_Vaccine == vacc)
            
            if (np.sum(where_source)!=0) and ("%s (%s)"%(vacc, source) not in Delta_done):
                
                #if (not re.search("Feikin", source)):
                days_fitting.append(Load_Delta["Days (Mid)"].values[where_source])
                ve_fitting.append(transform(Load_Delta["VE (value)"].values[where_source]))
                
                All_Days_Delta = np.concatenate((All_Days_Delta, Load_Delta["Days (Mid)"].values[where_source]))
                All_Delta_Data = np.concatenate((All_Delta_Data, Load_Delta["VE (value)"].values[where_source]))
                All_Days_xerr_Delta = np.concatenate((All_Days_xerr_Delta, Load_Delta["Days Err (+/-)"].values[where_source]))
                
                x = Load_Delta["VE (value)"].values[where_source]
                CI = np.array([Load_Delta["VE (Lower CI)"].values[where_source], Load_Delta["VE (Upper CI)"].values[where_source]]).T
                yerr, upper_CI, lower_CI = extract_yerr(x = x, CI = CI)
                
                if len(All_Delta_yerr.flatten()) !=0:
                    All_Delta_yerr = np.concatenate((All_Delta_yerr , yerr), axis = 1)
                else:
                    All_Delta_yerr = yerr
            
            Delta_done.append("%s (%s)"%(vacc, source))
    	
    def Fitting_IC50(thalf, tmax, t, Ab_classes, Cross_dic, quiet = False):  
        N = len(Ab_classes)
        t_max = tmax*np.ones(N) # time to the peak of antibody concentration for each antibody class
        t_half = thalf*np.ones(N) # antibody half-life for all antibody classes, respecively
        params_dic = {"t_max":t_max, "t_half":t_half}
        ### Compute  PK 
        is_log = False # if True, it returns the log of the antibody concentration
        solver = "lm" # root solver method for finding absorption rate ka (see scipy.optimize.root)
        c_t, c_dframe, ka, ke, c_max_dic = Antibody(t = t, params_dic = params_dic, is_log = is_log, Ab_names = Ab_classes, ka_solver = solver)
        
        #if not quiet:
        #    print("t_max = %.3f, t_half = %.3f"%(t_max[0], t_half[0]),"\n k_a:", ka, "\n k_e:", ke, "\n c_max", c_max_dic)
        
        ### select FR data for delta computed in the cross_reac_dic_show
        where_wt = list(variants_x_names_show).index("Wuhan-Hu-1")
        where_delta = list(variants_x_names_show).index("Delta: B.1.617.2")
        FR_delta = [Cross_dic[Ab_classes[i]][where_wt, where_delta] for i in range(len(Ab_classes))]
        ### Estimate one IC50 for all ABs
        guess = 0.3
        FC_ic50_list = [FC_ic50_dic[Ab_classes[i]] for i in range(len(Ab_classes))]
        IC50_data = root(sqrt_diff_FR, guess, args = (days_fitting, np.array(FC_ic50_list)*np.array(FR_delta), ve_fitting, len(Ab_classes), c_dframe), method = "lm").x
        #print("fitted", IC50_data)
        
        """Extract the IC50xx -- conserving the fold changes deviation from the mean of each epitope classes"""
        IC50xx = {Ab_classes[i]:IC50_data[0]*FC_ic50_dic[Ab_classes[i]] for i in range(len(Ab_classes))} 
        #print("IC50xx", IC50xx)
        
        return IC50xx, IC50_data[0], FC_ic50_list, FR_delta, c_dframe
    
    def Find_IC50_ranges(thalf_vec, tmax_vec, t, Ab_classes, Cross_dic):
        IC50xx_dic = {}
        for m in range(len(thalf_vec)):
            for n in range(len(tmax_vec)):
                thalf = thalf_vec[m]
                tmax = tmax_vec[n]
                IC50xx, IC50_data, FC_ic50_list, FR_delta, c_dframe = Fitting_IC50(thalf, tmax, t, Ab_classes, Cross_dic, quiet = True)           
                IC50xx_dic["(%d, %d)"%(m, n)] = IC50_data
    
        IC50xx = np.mean(list(IC50xx_dic.values()))
        mean_IC50xx_dic = {Ab_classes[i]:IC50xx*FC_ic50_dic[Ab_classes[i]] for i in range(len(Ab_classes))}
        return IC50xx_dic, mean_IC50xx_dic
    

    def PNeut_Envelope(s, t, variants, variant_x_names, Cross_react_dic, c_dframe_dic, IC50xx_dic, antigen_list = ["Wuhan-Hu-1"],mean_IC50xx = True):
        
        if mean_IC50xx:
            IC50xx = np.mean(list(IC50xx_dic.values()))
            to_print = "Computing P_Neut, used mean fitted IC50 %.5f"%IC50xx
        else:
            to_print = "Compute P_Neut, used fitted IC5 for each PK paramteters"
        
        to_print = to_print + " for %s vs. %s antigen"%("/".join(variants), antigen_list[s])
        
        num = "%d/%d"%(s+1, len(antigen_list))
        to_print = to_print + " (%s)"%num
        #print(to_print)
        
        splited_var = np.array(antigen_list[s].split("/"))
        splited_var = splited_var[~(splited_var == "")]
        splited_var = splited_var[~(splited_var == " ")]
        
        res = np.zeros((len(variants), len(list(c_dframe_dic.keys())), len(t)))
        ignore = np.zeros(res.shape).astype(bool)
        success = []
        for j in range(len(splited_var)):
            spl_sub = np.array(splited_var[j].split("="))
            spl_sub = spl_sub[~(spl_sub == "")]
            spl_sub = spl_sub[~(spl_sub == " ")]
    
            lin = spl_sub[0]
            if lin in variant_x_names:
                where_x = list(variant_x_names).index(lin)
            else:
                where_x = "Not Found"
                print("Error in antigen parameter: %s is not present in covsonar datat thus or Cross ALL simulated"%lin)
            
            if where_x != "Not Found":
                success.append(True)
                if len(spl_sub)==1:
                    prop_lin = 1/len(splited_var)
                else:
                    prop_0 = re.findall(r"[-+]?(?:\d*\.*\d+)", spl_sub[1])[0] ### anything else is error
                    prop_lin = float(prop_0)
                
                for i in range(len(variants)):
                    if variants[i] in variant_x_names:
                        where_y = list(variant_x_names).index(variants[i])
                    
                        for j in range(len(list(c_dframe_dic.keys()))):
                            if not mean_IC50xx:
                                IC50xx = IC50xx_dic[list(c_dframe_dic.keys())[j]]
                            
                            IC50xy = [Cross_react_dic[Ab_classes[i]][where_x, where_y]*IC50xx*FC_ic50_dic[Ab_classes[i]] for i in range(len(Ab_classes))]
                                
                            c_dframe = c_dframe_dic[list(c_dframe_dic.keys())[j]]
                            for l in range(len(t)): 
                                antibody_level = c_dframe.loc[l][1:]
                                #res[i, j, l] += prop_lin*efficacy_n_antibodies(antibody_level, IC50xy)
                                res[i, j, l] += prop_lin*(1 - np.prod(1 - (antibody_level/(antibody_level + np.array(IC50xy)))))
                    else:
                        #print("P_neut not computed for %s because it is not present in Cross ALL simulated"%variants[i])
                        ignore[i, j, :] = True
            else:
                success.append(False)
                
        if np.all(success):
            res = ma.masked_array(res, mask = ignore)
            return np.min(res, axis = (0, 1)), np.max(res, axis = (0, 1))
        else:
            return None, None
        
    """Compute Antibody concentration over time for a range of t_half and t_max"""
    thalf_vec = np.linspace(25, 69, 15) 
    tmax_vec = np.linspace(14, 28, 5)
    t_conc = np.arange(1, 656, 1) ### used in older codes
    c_t_vec, c_dframe_dic, dataname = Antibody_ranges(thalf_vec, tmax_vec, t_conc, Ab_classes)
    IC50xx_dic, mean_IC50xx_dic = Find_IC50_ranges(thalf_vec, tmax_vec, t_conc, Ab_classes,  Cross_with_delta_validation)
    
    """Save Dynamics: The Antibody PK is assumed to be the same for all Epitope Classes"""
    PK = {}
    PK["Day since activation"] = t_conc
    for key in c_dframe_dic.keys():
        key_num = np.array(re.findall(r"\d+", key)).astype(int)
        PK["t_half = %.3f \nt_max = %.3f"%(thalf_vec[key_num[0]], tmax_vec[key_num[1]])] = c_dframe_dic[key]["A"]
    
    PK_df = pd.DataFrame(PK)
   
    
    def pneut_util(Lin_name, variants_in_cross, antigen_list,
                   Cross_react_dic = None, var_name = None):
        
        PK_df.to_csv(result_location+"/PK_for_all_Epitopes_%s.csv"%Lin_name)
        if antigen_list == ["none"]:
            antigen_list = []
        
        if len(antigen_list) != 0:
            
            VE = {}
            VE["Day since infection"] = t_conc
            pfunc = partial(PNeut_Envelope, t=t_conc, 
                                variants=[Lin_name], 
                                variant_x_names = variants_in_cross, 
                                Cross_react_dic = Cross_react_dic,
                                c_dframe_dic = c_dframe_dic, 
                                IC50xx_dic = IC50xx_dic, 
                                antigen_list = antigen_list, 
                                mean_IC50xx = True, 
                                ) 
            status = False
            try:
                jb_res = list(jb.Parallel(n_jobs = -1, backend = "loky")(jb.delayed(pfunc)(d) for d in range(len(antigen_list))))
                status = True
                #print("run joblib.Parallel")
            except:
                try:
                    jb_res = list(jb.Parallel(n_jobs = -1, backend = "multiprocessing")(jb.delayed(pfunc)(d) for d in range(len(antigen_list))))
                    status=True
                    #print("run joblib.Parallel")
                except:
                    jb_res = list(jb.Parallel(n_jobs = -1, prefer = "threads")(jb.delayed(pfunc)(d) for d in range(len(antigen_list))))
                    status=True
                    #print("run joblib.Parallel")
            
            if status:
                 for i in range(len(antigen_list)): ## "Wuhan-Hu-1" is always in each cross reactivity files produced by our pipeline
                     antigen = antigen_list[i]
                     EnvD_Min,EnvD_Max = jb_res[i]
                     if EnvD_Min is not None:
                         VE["Proba Neut Min\n vs. %s antigen"%antigen] = EnvD_Min
                         VE["Proba Neut Max\n vs. %s antigen"%antigen] = EnvD_Max
            
            else:
                #print("joblib.Parallel failed running, using brute force looping")
                for i in range(len(antigen_list)): ## "Wuhan-Hu-1" is always in each cross reactivity files produced by our pipeline
                    antigen = antigen_list[i]
                    EnvD_Min,EnvD_Max = PNeut_Envelope(1, t_conc, [Lin_name], variants_in_cross, Cross_react_dic, c_dframe_dic, IC50xx_dic, antigen_list = antigen_list, mean_IC50xx = True)
                    if EnvD_Min is not None:
                        VE["Proba Neut Min\n vs. %s antigen"%antigen] = EnvD_Min
                        VE["Proba Neut Max\n vs. %s antigen"%antigen] = EnvD_Max
                
            """ Save P_Neut ranges"""
            VE_df = pd.DataFrame(VE)
            VE_df.to_csv(result_location+"/P_neut_"+Lin_name+".csv")
        
    
    #blockPrint()
    for lin in Lin_list:
        pneut_util(lin, variants_in_cross, antigen_list, Cross_react_dic)
    #enablePrint()