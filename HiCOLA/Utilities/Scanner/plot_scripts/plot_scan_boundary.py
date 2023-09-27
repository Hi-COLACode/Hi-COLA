#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 03:11:29 2023

@author: ashimsg
"""


EdS = 0.6423229
k1seed = -5.16652

import HiCOLA.Frontend.expression_builder as eb
import HiCOLA.Utilities.Other.alpha_and_sound_speed as al
import HiCOLA.Frontend.numerical_solver as ns
import sympy as sym
import numpy as np
import matplotlib.pyplot as plt
import glob
from HiCOLA.Utilities.Other.support import compute_fphi, ESS_direct_to_seed, ESS_seed_to_direct
from HiCOLA.Frontend.read_parameters import read_scan_result
#import HiCOLA.Utilities.Scanner.scan_to_plot as stp

to_exec = eb.declare_symbols()
exec(to_exec)
Xreal = 0.5*(E**2.)*phiprime**2.

k1, k2, g31, g32= sym.symbols('k_1 k_2 g_{31} g_{32}')
symbol_list = [k1, k2, g31, g32]
odeint_parameter_symbols = [E, phi, phiprime, omegar, omegam]
closure_declaration = ['odeint_parameters',2]

def find_nearest_seed_result(k1seed, EdS, seedlist):
    seedlist = np.array(seedlist)
    k1sdiff_arr = [np.abs(k1seed - k1s) for k1s in seedlist[:,0]]
    #print(k1sdiff_arr)
    minindices = np.where(k1sdiff_arr == np.min(k1sdiff_arr))[0]#[index for index, diff in enumerate(k1sdiff_arr) if diff == np.min(k1sdiff_arr)]
    #print(minindices)
    k1s_matches = np.array([seedlist[minind] for minind in minindices])
    #print(k1s_matches)
    EdSdiff_arr = [np.abs(EdSv - EdS) for EdSv in k1s_matches[:,4]]
    minindex = np.array(EdSdiff_arr).argmin()
    best_match = k1s_matches[minindex]
    return best_match

directory = '/home/ashimsg/Documents/QMUL_Desktop/Horndeski_COLA/Hi-COLA/Output/Scanner/from_apocrita/2023-06-14/00-51_45/'

colours = ['green', 'grey', 'black', 'magenta', 'pink', 'red', 'blue', 'yellow']

files = glob.glob(directory+'*.txt')

#collect only filenames that pertain to scan results
results = {}
for colour in colours:
    for file in files:
        if colour in file:
            results.update({colour:file})
            
green_results = {}
grey_results = {}
black_results = {}
magenta_results = {}
pink_results = {}
red_results = {}
blue_results = {}
yellow_results = {}

for key, result in results.items():
    
    if key == 'green':
        result_dict = read_scan_result(results[key])
        green_results.update(result_dict)
    # if key == 'grey':
    #     result_dict = read_scan_result(results[key])
    #     grey_results.update(result_dict)
    # if key == 'black':
    #     result_dict = read_scan_result(results[key])
    #     black_results.update(result_dict)
    # if key == 'magenta':
    #     result_dict = read_scan_result(results[key])
    #     magenta_results.update(result_dict)
    if key == 'pink':
        result_dict = read_scan_result(results[key])
        pink_results.update(result_dict)
    # if key == 'red':
    #     result_dict = read_scan_result(results[key])
    #     red_results.update(result_dict)
    # if key == 'blue':
    #     result_dict = read_scan_result(results[key])
    #     blue_results.update(result_dict)
    # if key == 'yellow':
    #     result_dict = read_scan_result(results[key])
    #     yellow_results.update(result_dict)
            
        
green_EdS = [1/U0 for U0 in green_results['U0_arr'] ]
green_f = [compute_fphi(omega_l, omega_m, omega_r) for omega_l, omega_m, omega_r in zip(green_results['Omega_l0_arr'], green_results['Omega_m0_arr'], green_results['Omega_r0_arr']  )]
green_k1seed = [ESS_direct_to_seed(k1, g31, omega_l0, f_phi, EdS)[0] for k1, g31, omega_l0, f_phi, EdS in zip(green_results['k1_arr'], green_results['g31_arr'], green_results['Omega_l0_arr'], green_f, green_EdS )]
green_g31seed = [ESS_direct_to_seed(k1, g31, omega_l0, f_phi, EdS)[1] for k1, g31, omega_l0, f_phi, EdS in zip(green_results['k1_arr'], green_results['g31_arr'], green_results['Omega_l0_arr'], green_f, green_EdS )]
green_seeds = [ESS_direct_to_seed(k1, g31, omega_l0, f_phi, EdS) for k1, g31, omega_l0, f_phi, EdS in zip(green_results['k1_arr'], green_results['g31_arr'], green_results['Omega_l0_arr'], green_f, green_EdS )]

pink_EdS = [1/U0 for U0 in pink_results['U0_arr'] ]
pink_f = [compute_fphi(omega_l, omega_m, omega_r) for omega_l, omega_m, omega_r in zip(pink_results['Omega_l0_arr'], pink_results['Omega_m0_arr'], pink_results['Omega_r0_arr']  )]
pink_k1seed = [ESS_direct_to_seed(k1, g31, omega_l0, f_phi, EdS)[0] for k1, g31, omega_l0, f_phi, EdS in zip(pink_results['k1_arr'], pink_results['g31_arr'], pink_results['Omega_l0_arr'], pink_f, pink_EdS )]
pink_g31seed = [ESS_direct_to_seed(k1, g31, omega_l0, f_phi, EdS)[1] for k1, g31, omega_l0, f_phi, EdS in zip(pink_results['k1_arr'], pink_results['g31_arr'], pink_results['Omega_l0_arr'], pink_f, pink_EdS )]
pink_seeds = [ESS_direct_to_seed(k1, g31, omega_l0, f_phi, EdS) for k1, g31, omega_l0, f_phi, EdS in zip(pink_results['k1_arr'], pink_results['g31_arr'], pink_results['Omega_l0_arr'], pink_f, pink_EdS )]


nearest_seed_entry = find_nearest_seed_result(k1seed, EdS, pink_seeds)

def find_nearest_scan_result(nearest_seed_entry, result_dictionary):
    nearest_k1seed = nearest_seed_entry[0]
    nearest_g31seed = nearest_seed_entry[1]
    nearest_omegal0 = nearest_seed_entry[2]
    nearest_fphi = nearest_seed_entry[3]
    nearest_EdS = nearest_seed_entry[4]
    nearest_U0 = 1./nearest_EdS
    
    nearest_k1, nearest_g31 = ESS_seed_to_direct(nearest_k1seed, nearest_g31seed, nearest_omegal0, nearest_fphi, nearest_EdS)
    
    k1diff_arr = [np.abs(nearest_k1 - k1) for k1 in result_dictionary['k1_arr']]
    mink1indices = np.where(k1diff_arr == np.min(k1diff_arr))[0]#[index for index, diff in enumerate(k1diff_arr) if diff == np.min(k1diff_arr)]
    
    g31_k1indexed = np.array([result_dictionary['g31_arr'][minind] for minind in mink1indices])
    g31diff_arr = [np.abs(nearest_g31 - g31) for g31 in g31_k1indexed]
    ming31indices = np.where(g31diff_arr == np.min(g31diff_arr))[0]#[index for index, diff in enumerate(g31diff_arr) if diff == np.min(g31diff_arr)]
    
    U0_g31indexed = np.array([result_dictionary['U0_arr'][minind] for minind in ming31indices])
    U0diff_arr = [np.abs(nearest_U0 - U0) for U0 in U0_g31indexed]
    minU0index = np.array(U0diff_arr).argmin()
    
    nearest_scan_res = []
    for i in result_dictionary.values():
        nearest_scan_res.append(i[minU0index])
    U0, phi0, phiprime0, Omega_m0, Omega_r0, Omega_l0, fsolveier, k1, k2, g31, g32 = nearest_scan_res
    scan_res_reordered = U0, phi0, phiprime0, Omega_r0, Omega_m0, Omega_l0, fsolveier, [k1, k2, g31, g32]
    return scan_res_reordered

scan_entry = find_nearest_scan_result(nearest_seed_entry, pink_results)
print(scan_entry)