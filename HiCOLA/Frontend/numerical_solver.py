###################
# Loading modules #
###################

import numpy as np
import scipy.integrate as integrate
from scipy.integrate import odeint
from scipy.optimize import fsolve
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from matplotlib.widgets import Slider, Button, RadioButtons
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext
from HiCOLA.Frontend.expression_builder import *
import sympy as sym
import sys
import itertools as it
import time
import os
from HiCOLA.Utilities.Other.suppressor import *
#import timing #need to include the timing.py in working directory for this to work in apocrita


##################
# Plotting style #
##################

#plt.rcParams['lines.linewidth'] = 4
plt.rcParams['axes.linewidth'] = 2
plt.rcParams['font.size'] = 16
#plt.rc('xtick.major', pad=20)
#plt.rc('ytick.major', pad=15)
plt.rcParams['lines.linewidth'] = 4
plt.rcParams['xtick.labelsize'] = 16 #28#32
plt.rcParams['ytick.labelsize'] = 16 #28#32
plt.rcParams['axes.labelpad'] = 15
plt.rcParams["errorbar.capsize"] = 8
plt.rcParams["lines.markersize"] = 10 # 8 # 5
plt.rcParams["lines.markeredgewidth"] = 1
#plt.rcParams["axes.markeredgecolor"] = 'black'


##### Miscellaneous functions ####

def nearest_index(arr, val):
    '''
   Finds index of entry in arr with the nearest value *above* val

    '''
    arr = np.asarray(arr)
    subtracted_arr = arr - val
    for index,val in enumerate(subtracted_arr):
        newval=1000.
        if val < 0:
            subtracted_arr[index] = newval	
    near_index = subtracted_arr.argmin()
    return near_index



##########################
# Cosmological functions #
##########################

def comp_H_LCDM(z, Omega_r0, Omega_m0, H0):
    Omega_L0 = 1.-Omega_m0-Omega_r0
    H = H0*np.sqrt(Omega_m0*(1.+z)**3. + Omega_r0*(1.+z)**4. + Omega_L0)
    return H

def comp_E_LCDM(z, Omega_r0, Omega_m0):
    Omega_L0 = 1.-Omega_m0-Omega_r0
    E = np.sqrt(Omega_m0*(1.+z)**3. + Omega_r0*(1.+z)**4. + Omega_L0)
    return E

def comp_E_LCDM_DE(z, Omega_r0, Omega_m0):
    Omega_L0 = 1.-Omega_m0-Omega_r0
    E = np.sqrt(Omega_m0*(1.+z)**3. + Omega_r0*(1.+z)**4. + Omega_L0)
    return E

def comp_Omega_r_LCDM(z, Omega_r0, Omega_m0):
    Omega_L0 = 1.-Omega_m0-Omega_r0
    E = comp_E_LCDM(z, Omega_r0, Omega_m0)
    Omega_r = Omega_r0*(1.+z)**4./E/E
    return Omega_r

def comp_Omega_m_LCDM(z, Omega_r0, Omega_m0):
    Omega_L0 = 1.-Omega_m0-Omega_r0
    E = comp_E_LCDM(z, Omega_r0, Omega_m0)
    Omega_m = Omega_m0*(1.+z)**3./E/E
    return Omega_m

def comp_Omega_L_LCDM(z, Omega_r0, Omega_m0):
    Omega_L0 = 1.-Omega_m0-Omega_r0
    E = comp_E_LCDM(z, Omega_r0, Omega_m0)
    Omega_L = Omega_L0/E/E
    return Omega_L

def comp_Omega_DE_LCDM(x, Omega_r0, Omega_m0):
    Omega_DE0 = 1. - Omega_m0 - Omega_r0
    term1 = Omega_r0*np.exp(-4.*x)+Omega_m0*np.exp(-3.*x)+Omega_DE0
    Omega_DE = Omega_DE0/term1
    return Omega_DE

def comp_E_prime_E_LCDM(x, Omega_r0, Omega_m0):
    Omega_DE0 = 1. - Omega_m0 - Omega_r0
    term1 = Omega_r0*np.exp(-4.*x)+Omega_m0*np.exp(-3.*x)+Omega_DE0
    term2 = 4.*Omega_r0*np.exp(-4.*x)+3.*Omega_m0*np.exp(-3.*x)
    E_prime_E = -0.5*term2/term1
    return E_prime_E

def comp_Omega_DE_prime_LCDM(E_prime_E, Omega_DE):
    Omega_DE_prime = -2.*E_prime_E*Omega_DE
    return Omega_DE_prime

def comp_alpha_M_propto_Omega_DE_LCDM(c_M, Omega_DE):
    alpha_M_DE = c_M*Omega_DE
    return alpha_M_DE

def comp_alpha_M_prime_propto_Omega_DE_LCDM(c_M, Omega_DE_prime):
    alpha_M_prime_DE = c_M*Omega_DE_prime
    return alpha_M_prime_DE

def alpha_M_int_propto_Omega_DE_LCDM(x, Omega_r0, Omega_m0, c_M):
    Omega_DE = comp_Omega_DE_LCDM(x, Omega_r0, Omega_m0)
    alpha_M_int = c_M * Omega_DE
    return alpha_M_int

def comp_Meffsq_x2_x1_propto_Omega_DE_LCDM(x1, x2, Omega_r0, Omega_m0, c_M):
    Meffsq_x2_x1 = np.exp(integrate.quad(alpha_M_int_propto_Omega_DE_LCDM, x1, x2, args=(Omega_r0, Omega_m0, c_M))[0])
    return Meffsq_x2_x1

def comp_Omega_r_prime(Omega_r, E, E_prime):
    E_prime_E = E_prime/E
    Omega_r_prime = -Omega_r*(4.+2.*E_prime_E)
    return Omega_r_prime

def comp_Omega_l_prime(Omega_l0, E, Eprime):
    return -2.*Omega_l0*Eprime/E/E/E

def comp_Omega_m_prime(Omega_m, E, E_prime):
    E_prime_E = E_prime/E
    Omega_m_prime = -Omega_m*(3.+2.*E_prime_E)
    return Omega_m_prime

#fried1 is ultimately the closure equation for the density parameters
def fried1(phi_prime, k1, g1, Omega_r, Omega_m, E, alpha_M, Ms_Mp, Meffsq_Mpsq):#this is eq. 50 in Bill's non-tracker overleaf with the 1 moved over to right-side so that RHS =0
    zer0 = (Omega_r+Omega_m)/Meffsq_Mpsq
    zer1a = 3.*Meffsq_Mpsq/Ms_Mp/Ms_Mp
    zer1b = 0.5*k1*phi_prime**2.
    zer1c = 3.*g1*Ms_Mp*E**2.*phi_prime**3.
    zer1 = (zer1b+zer1c)/zer1a
    zer2 = -alpha_M
    zer = zer0 + zer1 + zer2 - 1.
    return zer


def fried_RHS_wrapper(cl_variable, cl_declaration, fried_RHS_lambda, E, phi, phi_prime, Omega_r, Omega_m, Omega_l, parameters):
    argument_no = 6 + len(parameters)
    if cl_declaration[1] > argument_no -1:
        raise Exception('Invalid declaration - there is no valid argument index for the declaration')
    if cl_declaration[0] == 'odeint_parameters':
        if cl_declaration[1] == 0:
            return fried_RHS_lambda(cl_variable, phi, phi_prime,Omega_r,Omega_m,Omega_l, *parameters) #Closure used to compute E0
        if cl_declaration[1] == 1:
            return fried_RHS_lambda(E, cl_variable,phi_prime, Omega_r,Omega_m,Omega_l,*parameters) #Closure used to compute phi0
        if cl_declaration[1] == 2:
            return fried_RHS_lambda(E, phi, cl_variable, Omega_r,Omega_m,Omega_l,*parameters) #Closure used to compute phi_prime0
        if cl_declaration[1] == 3:
            return fried_RHS_lambda(E,phi_prime,cl_variable,Omega_m,Omega_l,*parameters) #Closure used to compute Omega_r0
        if cl_declaration[1] == 4:
            return fried_RHS_lambda(E,phi_prime,Omega_r,cl_variable,Omega_l,*parameters) #Closure used to compute Omega_m0
        # if cl_declaration[1] == 5:
        #     return fried_RHS_lambda(E,phi,phi_prime,Omega_r,Omega_m,cl_variable,*parameters) #Closure used to compute Omega_l0
    if cl_declaration[0] == 'parameters':
        parameters[cl_declaration[1]] = cl_variable
        return fried_RHS_lambda(E, phi, phi_prime, Omega_r, Omega_m, Omega_l, *parameters)
    else:
        raise Exception('Invalid string in declaration list. Must be either \'odeint_parameters\' or \'parameters\'')



def comp_param_close(fried_closure_lambda, cl_declaration, E0, phi0, phi_prime0, Omega_r0, Omega_m0, Omega_l0,parameters): #calling the closure-fixed parameter "k1" is arbitrary, the choice of which parameter to fix is determined by lambdification or fried_RHS

    cl_guess = 1.0 #this may need to be changed depending on what is being solved for through closure, if fsolve has trouble
    if cl_declaration[0] == 'odeint_parameters':
        if cl_declaration[1] == 0:
            cl_guess = E0
        if cl_declaration[1] == 1:
            cl_guess = phi0
        if cl_declaration[1] == 2:
            cl_guess = phi_prime0
        if cl_declaration[1] == 3:
            cl_guess = Omega_r0
        if cl_declaration[1] == 4:
            cl_guess = Omega_m0
        # if cl_declaration[1] ==5:
        #     cl_guess = Omega_l0
    if cl_declaration[0] == 'parameters':
        cl_guess = parameters[cl_declaration[1]]


    cl_variable = fsolve(fried_RHS_wrapper, cl_guess, args=(cl_declaration, fried_closure_lambda, E0, phi0, phi_prime0, Omega_r0,Omega_m0, Omega_l0,parameters), xtol=1e-6,full_output=True) #make sure arguments match lambdification line in run_builder.py
    # print('f solve integer is '+str(fsolveier))
    # print('f solve message is '+fsolvemsg)
    # print('f solve found '+str(cl_variable))
    cl_variable0 = cl_variable[0][0] #to unpack the fsolve solution/iteration

    return cl_variable0, cl_variable

def comp_E_prime_E(k1, g1, Omega_r, E, phi_prime, alpha_M, alpha_M_prime, Ms_Mp, Meffsq_Mpsq):
    denom1 = g1*Ms_Mp**3.*E**2.*phi_prime**3./Meffsq_Mpsq - 2. - alpha_M
    denom2 = g1*E**2.*Ms_Mp**3.*phi_prime**2./Meffsq_Mpsq
    denom3 = k1 + 6.*g1*Ms_Mp*E**2.*phi_prime
    denom4a = -phi_prime*(k1+6.*g1*Ms_Mp*E**2.*phi_prime)
    denom4b = 3.*Meffsq_Mpsq*alpha_M/phi_prime/Ms_Mp/Ms_Mp
    denom4c = -3.*g1*Ms_Mp*E**2.*phi_prime**2.
    denom4 = denom4a+denom4b+denom4c
    if denom3 == 0:
        print('warning, zeros: ', k1, 6.*g1*Ms_Mp*E**2.*phi_prime, denom3)
    denom = 1. + denom2*denom4/denom1/denom3
    numer1 = g1*Ms_Mp**3.*E**2.*phi_prime**3./Meffsq_Mpsq - 2. - alpha_M
    numer2 = 3. + 0.5*k1*Ms_Mp**2.*phi_prime**2./Meffsq_Mpsq
    numer3 = 2.*alpha_M + alpha_M**2. + alpha_M_prime + Omega_r/Meffsq_Mpsq
    numer4a = -g1*Ms_Mp**3.*E**2.*phi_prime**2./(Meffsq_Mpsq**1.5) #Meffsq_Mpsq?
    numer4b = k1 + 6.*g1*Ms_Mp*E**2.*phi_prime
    numer4c = 6.*Meffsq_Mpsq*alpha_M/phi_prime/Ms_Mp/Ms_Mp - 3.*phi_prime*k1
    numer4d = -9.*g1*Ms_Mp*E**2.*phi_prime**2.
    numer4 = numer4a*(numer4c+numer4d)/numer4b
    numer = (numer2+numer3+numer4)/numer1
    E_prime_E = numer/denom
    return E_prime_E


def comp_phi_primeprime(k1, g1, Omega_r, E, E_prime, phi_prime, alpha_M, alpha_M_prime, Ms_Mp, Meffsq_Mpsq):
    E_prime_E = E_prime/E
    term1 = -E_prime_E*phi_prime
    numer2 = 3.*Meffsq_Mpsq*(E_prime_E+2.)*alpha_M
    denom2 = Ms_Mp**2.*phi_prime*(k1+6.*g1*Ms_Mp*E**2.*phi_prime)
    term2 = numer2/denom2
    numer3 = -3.*k1*phi_prime
    denom3 = k1 + 6.*g1*Ms_Mp*E**2.*phi_prime
    term3 = numer3/denom3
    numer4 = -3.*g1*Ms_Mp*E**2.*phi_prime**2.*(E_prime_E+3.)
    denom4 = k1 + 6.*g1*Ms_Mp*E**2.*phi_prime
    term4 = numer4/denom4
    phi_primeprime = term1 + term2 + term3 + term4
    return phi_primeprime

def comp_primes(x, Y, E0, Omega_r0, Omega_m0, Omega_l0, E_prime_E_lambda, E_prime_E_safelambda, phi_primeprime_lambda, phi_primeprime_safelambda, A_lambda, cl_declaration, parameters,threshold=1e-3,GR_flag=False): #x, Y swapped for solve_ivp ###ADD LAMBDA FUNCTION AS ARGUMENT###
    '''
    x: the
    '''
    a = np.exp(x) #unused, this (or x) was previously used to compute LCDM quantities, which are no longer used.
    #print('a is '+str(a))
    #print([Y,a],file=Y_vs_a_txt)
    #print('Y is '+str(Y))
    phiY, phi_primeY, EUY, Omega_rY, Omega_mY, Omega_lY = Y
    A_value = A_lambda(EUY,phiY,phi_primeY,*parameters)
    if A_value - abs(A_value) == 0:
        A_sign = 1.
    elif A_value - abs(A_value) != 0:
        A_sign = -1.


    if (abs(A_value) >= threshold and GR_flag==False) or (threshold==0. and GR_flag==False):
        E_prime_E_evaluated = E_prime_E_lambda(EUY,phiY, phi_primeY,Omega_rY,Omega_lY,*parameters)
        E_prime_evaluated = E_prime_E_evaluated*EUY
        phi_primeprime_evaluated = phi_primeprime_lambda(EUY,E_prime_evaluated, phiY, phi_primeY,*parameters)
    if (abs(A_value) < threshold and GR_flag==False):
        E_prime_E_evaluated = E_prime_E_safelambda(EUY,phiY, phi_primeY,Omega_rY,Omega_lY, threshold,A_sign,*parameters)
        E_prime_evaluated = E_prime_E_evaluated*EUY
        phi_primeprime_evaluated = phi_primeprime_safelambda(EUY,E_prime_evaluated,phiY, phi_primeY,threshold,A_sign,*parameters)
    if GR_flag==True:
        E_prime_E_evaluated = comp_E_prime_E_LCDM(x,Omega_r0,Omega_m0)
        E_prime_evaluated = E_prime_E_evaluated*EUY
        phi_primeprime_evaluated = 0.
    if cl_declaration[0] == 'odeint_parameters': #usually indicates dS approach, so we must convert U back to E, since this is what the Omega_prime functions use
        EY = EUY/E0
        EYprime = E_prime_evaluated/E0
    if cl_declaration[0] == 'parameters': #usually indicates 'today' approach, no need to change the Hubble variable, it is already E
        EY = EUY
        EYprime = E_prime_evaluated
    Omega_r_prime = comp_Omega_r_prime(Omega_rY, EY, EYprime)
    Omega_m_prime = comp_Omega_m_prime(Omega_mY, EY, EYprime)
    Omega_l_prime = comp_Omega_l_prime(Omega_l0,EY, EYprime)
    return [phi_primeY, phi_primeprime_evaluated, E_prime_evaluated, Omega_r_prime, Omega_m_prime, Omega_l_prime]

def chi_over_delta(a_arr, E_arr, calB_arr, calC_arr, Omega_m0): #the E_arr is actual E, not U! Convert U_arr to E_arr!
    chioverdelta = np.array(calB_arr)*np.array(calC_arr)*Omega_m0/np.array(E_arr)/np.array(E_arr)/np.array(a_arr)/np.array(a_arr)/np.array(a_arr)
    return chioverdelta

def beta_DGP(rcH0, E, Eprime):
    beta = 1 + 2*E*rcH0*(1 + Eprime/3./E)
    return beta

def coupling_DGP(rcH0, E, Eprime):
    beta = beta_DGP(rcH0, E, Eprime)
    return 1./3./beta

def chi_over_delta_DGP(rcH0, a, E, Eprime, Omega_m0):
    beta = beta_DGP(rcH0, E, Eprime)
    chioverdelta = 8.*rcH0*rcH0*Omega_m0/9./beta/beta/a/a/a
    return chioverdelta

def coupling_kmou(a,betaK):
    arr = np.ones(len(a))
    return arr*2*betaK*betaK

def screening_kmou(a,K0,H0,betaK,M_pG4,M_KG4):
    term1 = 6.*abs(K0)*((2/3)**3.)
    term2 = 3*betaK*M_pG4/2./np.pi/H0/M_KG4
    term2m = term2*term2
    screening = (term1*term2m)**0.25

    arr = np.ones(len(a))
    return arr*screening


def run_solver(read_out_dict):
# (z_num, z_ini, Hubble0, phi_prime0, E_prime_E_lambda, E_prime_E_safelambda, phi_primeprime_lambda, phi_primeprime_safelambda,
#  omega_phi_lambda, A_lambda, fried_RHS_lambda, calB_lambda, calC_lambda, coupling_factor, cl_declaration, parameters,
#  Omega_r0, Omega_m0, Omega_l0, parameter_symbols, odeint_parameter_symbols, suppression_flag, threshold,GR_flag):


    [Omega_r0, Omega_m0, Omega_l0] = read_out_dict['cosmological_parameters']
    [Hubble0, phi0, phi_prime0] = read_out_dict['initial_conditions']
    [Npoints, z_max, suppression_flag, threshold, GR_flag] = read_out_dict['simulation_parameters']
    parameters = read_out_dict['Horndeski_parameters']

    
    E_prime_E_lambda = read_out_dict['E_prime_E_lambda']
    E_prime_E_safelambda = read_out_dict['E_prime_E_safelambda']
    phi_primeprime_lambda = read_out_dict['phi_primeprime_lambda']
    phi_primeprime_safelambda = read_out_dict['phi_primeprime_safelambda']
    omega_phi_lambda = read_out_dict['omega_phi_lambda']
    A_lambda = read_out_dict['A_lambda']
    fried_RHS_lambda = read_out_dict['fried_RHS_lambda']
    calB_lambda = read_out_dict['calB_lambda']
    calC_lambda = read_out_dict['calC_lambda']
    coupling_factor = read_out_dict['coupling_factor']

    #debug
    EpE_debug = read_out_dict['EpE_debug']
    omega_phi_debug = read_out_dict['omega_phi_debug']
    K_lambda = read_out_dict['Horndeski_K']

    parameter_symbols = read_out_dict['symbol_list']
    odeint_parameter_symbols = read_out_dict['odeint_parameter_symbols']
    cl_declaration = read_out_dict['closure_declaration']

    z_final = 0.
    x_ini = np.log(1./(1.+z_max))
    x_final = np.log(1./(1.+z_final))
    x_arr = np.linspace(x_ini, x_final, Npoints)
    a_arr = [np.exp(x) for x in x_arr]
    x_arr_inv = x_arr
    a_arr_inv = a_arr


    # print('phi_prime imminent')
    if GR_flag is True:
        phi_prime0 = 0.
    #print('phi prime0 is '+str(phi_prime0))

    #omega_i(z=z_max) for forward solver ICs
    Omega_l_ini = Omega_l0/Hubble0/Hubble0
    Omega_m_ini = Omega_m0*((1+z_max)**3.)/Hubble0/Hubble0
    Omega_r_ini = Omega_r0*((1+z_max)**4)/Hubble0/Hubble0

    cl_var, cl_var_full = comp_param_close(fried_RHS_lambda, cl_declaration, Hubble0, phi0, phi_prime0, Omega_r_ini, Omega_m_ini, Omega_l_ini, parameters)

    if cl_declaration[0] == 'odeint_parameters':
        if cl_declaration[1] == 0:
            Hubble0_closed = cl_var
            Y0 = [phi0, phi_prime0,Hubble0_closed,Omega_r_ini,Omega_m_ini, Omega_l_ini]
        if cl_declaration[1] == 1:
            phi0_closed = cl_var
            Y0 = [phi0_closed, phi_prime0, Hubble0, Omega_r_ini,Omega_m_ini, Omega_l_ini]
        if cl_declaration[1] == 2:
            phi_prime0_closed = cl_var
            if 1.-Omega_r0 - Omega_m0 == Omega_l0:
                    phi_prime0_closed = 0.
            Y0 = [phi0, phi_prime0_closed,Hubble0,Omega_r_ini,Omega_m_ini, Omega_l_ini]
        if cl_declaration[1] == 3:
            Omega_r_ini_closed = cl_var
            Y0 =[phi0, phi_prime0,Hubble0,Omega_r_ini_closed,Omega_m_ini, Omega_l_ini]
        if cl_declaration[1] == 4:
            Omega_m_ini_closed = cl_var
            Y0 = [phi0, phi_prime0,Hubble0,Omega_r_ini,Omega_m_ini_closed, Omega_l_ini]
        #print('Closure parameter is '+ str(odeint_parameter_symbols[cl_declaration[1]])+' = ' +str(cl_var))
    if cl_declaration[0] == 'parameters':
        parameters[cl_declaration[1]] = cl_var
        Y0 = [phi0, phi_prime0,Hubble0,Omega_r_ini,Omega_m_ini, Omega_l_ini]
        #print('Closure parameter is '+ str(parameter_symbols[cl_declaration[1]])+' = ' +str(cl_var))
    # print('Y0 is '+str(Y0))

    if suppression_flag is True:
        with stdout_redirected():
            ans = odeint(comp_primes, Y0, x_arr_inv, args=(Hubble0, Omega_r_ini,Omega_m_ini, Omega_l_ini, E_prime_E_lambda, E_prime_E_safelambda, phi_primeprime_lambda, phi_primeprime_safelambda, A_lambda, cl_declaration, parameters,threshold,GR_flag), tfirst=True)#, hmax=hmaxv) #k1=-6, g1 = 2
    else:
        ans = odeint(comp_primes, Y0, x_arr_inv, args=(Hubble0, Omega_r_ini,Omega_m_ini, Omega_l_ini, E_prime_E_lambda, E_prime_E_safelambda, phi_primeprime_lambda, phi_primeprime_safelambda, A_lambda, cl_declaration, parameters,threshold,GR_flag), tfirst=True)#, hmax=hmaxv)
    #ans = solve_ivp(comp_primes, t_eval=x_arr_inv, t_span=(x_arr_inv[0],x_arr_inv[-1]), y0=Y0, args=(Omega_r0, Omega_m0, c_M), method='Radau',max_step=hmaxv) #k1=-6, g1 = 2
    #print('odeint stage computed')
    phi_arr = ans[:,0]
    phi_prime_arr = ans[:,1]
    Hubble_arr = ans[:,2]

    Omega_r_arr = ans[:,3]
    Omega_m_arr = ans[:,4]
    Omega_l_arr = ans[:,5]

    # phi_prime_arr = ans.y[0]
    # Hubble_arr = ans.y[1]
    # print('shape of Hubble_arr is '+str(np.shape(Hubble_arr)))
    # print('Hubble_arr is '+str(Hubble_arr[0]))
    # Omega_r_arr = ans.y[2]
    # Omega_m_arr = ans.y[3]


    E_prime_E_LCDM_arr = [comp_E_prime_E_LCDM(xv, Omega_r0, Omega_m0) for xv in x_arr_inv]
    Omega_DE_LCDM_arr = [comp_Omega_DE_LCDM(xv, Omega_r0, Omega_m0) for xv in x_arr_inv]
    Omega_DE_prime_LCDM_arr = [comp_Omega_DE_prime_LCDM(E_prime_E_LCDMv, Omega_DE_LCDMv) for E_prime_E_LCDMv, Omega_DE_LCDMv in zip(E_prime_E_LCDM_arr, Omega_DE_LCDM_arr)]

    E_prime_E_arr = []
    phi_primeprime_arr = []
    for Omega_rv, Omega_lv, Ev, phiv, phi_primev in zip(Omega_r_arr, Omega_l_arr, Hubble_arr, phi_arr, phi_prime_arr):
        if GR_flag is True:
            E_prime_E_arr = E_prime_E_LCDM_arr
        else:
            E_prime_E_arr.append(E_prime_E_lambda(Ev,phiv, phi_primev,Omega_rv, Omega_lv, *parameters))
    Hubble_prime_arr = [E_prime_Ev*Ev for E_prime_Ev, Ev in zip(E_prime_E_arr, Hubble_arr)]
    for Omega_rv, Ev, E_primev, phiv, phi_primev in zip(Omega_r_arr, Hubble_arr, Hubble_prime_arr, phi_arr, phi_prime_arr):
        phi_primeprime_arr.append(phi_primeprime_lambda(Ev,E_primev,phiv, phi_primev,*parameters))

    A_arr = []
    for Ev, phiv, phi_primev,  in zip(Hubble_arr,phi_arr,phi_prime_arr):
        A_arr.append(A_lambda(Ev, phiv, phi_primev, *parameters))

    Omega_phi_arr = []
    for Ev, phiv, phiprimev, omegalv, omegamv, omegarv in zip(Hubble_arr,phi_arr,phi_prime_arr, Omega_l_arr, Omega_m_arr, Omega_r_arr):
        Omega_phi_arr.append(omega_phi_lambda(Ev,phiv,phiprimev,omegalv, omegamv, omegarv,*parameters))
    #######

    Omega_DE_arr = []
    Omega_phi_diff_arr = []
    Omega_r_prime_arr = []
    Omega_m_prime_arr= []
    Omega_l_prime_arr = []
    for Ev, E_primev, Omega_rv, Omega_mv, Omega_lv in zip(Hubble_arr, Hubble_prime_arr, Omega_r_arr,Omega_m_arr, Omega_l_arr):
        Omega_DE_arr.append(1. - Omega_rv - Omega_mv)
        Omega_phi_diff_arr.append(1. - Omega_rv - Omega_mv - Omega_lv)
        Omega_r_prime_arr.append(comp_Omega_r_prime(Omega_rv, Ev, E_primev))
        Omega_m_prime_arr.append(comp_Omega_m_prime(Omega_mv, Ev, E_primev))
        Omega_l_prime_arr.append(comp_Omega_l_prime(Omega_l0,Ev, E_primev))


    array_output = []
    for i in [a_arr_inv, Hubble_arr, E_prime_E_arr, Hubble_prime_arr, phi_prime_arr,  phi_primeprime_arr, Omega_r_arr, Omega_m_arr, Omega_DE_arr, Omega_l_arr, Omega_phi_arr, Omega_phi_diff_arr, Omega_r_prime_arr, Omega_m_prime_arr, Omega_l_prime_arr, A_arr]:
        i = np.array(i)
        array_output.append(i)
    [a_arr_inv, Hubble_arr, E_prime_E_arr, Hubble_prime_arr, phi_prime_arr,  phi_primeprime_arr, Omega_r_arr, Omega_m_arr, Omega_DE_arr, Omega_l_arr, Omega_phi_arr, Omega_phi_diff_arr, Omega_r_prime_arr, Omega_m_prime_arr, Omega_l_prime_arr, A_arr] = array_output

    calB_arr = []
    calC_arr = []
    coupling_factor_arr = []
    for UEv, UEprimev, phiv, phiprimev, phiprimeprimev in zip(Hubble_arr, Hubble_prime_arr, phi_arr, phi_prime_arr, phi_primeprime_arr):
        calB_arr.append(calB_lambda(UEv,UEprimev,phiv,phiprimev,phiprimeprimev, *parameters))
        calC_arr.append(calC_lambda(UEv,UEprimev,phiv,phiprimev,phiprimeprimev, *parameters))
        coupling_factor_arr.append(coupling_factor(UEv,UEprimev,phiv,phiprimev,phiprimeprimev,*parameters))

    E_arr = Hubble_arr/Hubble0
    chioverdelta_arr = chi_over_delta(a_arr_inv, E_arr, calB_arr, calC_arr, Omega_m0)

    #debug
    debug_dict = {}
    Adebug_arr = []
    A1debug_arr = []
    A2debug_arr = []
    B1debug_arr = []
    B2debug_arr = []
    B21debug_arr = []
    B22debug_arr = []
    B23debug_arr = []
    EpE_term1_arr = []
    EpE_term21_arr = []
    EpE_term22_arr = []
    EpE_term2_arr = []
    EpE_term3_arr = []
    EpE_term4_arr = []
    for key, term in EpE_debug.items():
        for Omega_rv, Omega_lv, Ev, phiv, phi_primev in zip(Omega_r_arr, Omega_l_arr, Hubble_arr, phi_arr, phi_prime_arr):
            if key == 'A':
                Adebug_arr.append(term(Ev,phiv, phi_primev,Omega_rv, Omega_lv, *parameters))
                debug_dict.update({key:Adebug_arr})
            if key == 'A1':
                A1debug_arr.append(term(Ev,phiv, phi_primev,Omega_rv, Omega_lv, *parameters))
                debug_dict.update({key:A1debug_arr})
            if key == 'A2':
                A2debug_arr.append(term(Ev,phiv, phi_primev,Omega_rv, Omega_lv, *parameters))
                debug_dict.update({key:A2debug_arr})
            if key == 'B1':
                B1debug_arr.append(term(Ev,phiv, phi_primev,Omega_rv, Omega_lv, *parameters))
                debug_dict.update({key:B1debug_arr})
            if key == 'B2':
                B2debug_arr.append(term(Ev,phiv, phi_primev,Omega_rv, Omega_lv, *parameters))
                debug_dict.update({key:B2debug_arr})
            if key == 'B21':
                B21debug_arr.append(term(Ev,phiv, phi_primev,Omega_rv, Omega_lv, *parameters))
                debug_dict.update({key:B21debug_arr})
            if key == 'B22':
                B22debug_arr.append(term(Ev,phiv, phi_primev,Omega_rv, Omega_lv, *parameters))
                debug_dict.update({key:B22debug_arr})
            if key == 'B23':
                B23debug_arr.append(term(Ev,phiv, phi_primev,Omega_rv, Omega_lv, *parameters))
                debug_dict.update({key:B23debug_arr})
            if key == 'term1':
                EpE_term1_arr.append(term(Ev,phiv, phi_primev,Omega_rv, Omega_lv, *parameters))
                debug_dict.update({key:EpE_term1_arr})
            if key == 'term21':
                EpE_term21_arr.append(term(Ev,phiv, phi_primev,Omega_rv, Omega_lv, *parameters))
                debug_dict.update({key:EpE_term21_arr})
            if key =='term22':
                EpE_term22_arr.append(term(Ev,phiv, phi_primev,Omega_rv, Omega_lv, *parameters))
                debug_dict.update({key:EpE_term22_arr})
            if key == 'term2':
                EpE_term2_arr.append(term(Ev,phiv, phi_primev,Omega_rv, Omega_lv, *parameters))
                debug_dict.update({key:EpE_term2_arr})
            if key == 'term3':
                EpE_term3_arr.append(term(Ev,phiv, phi_primev,Omega_rv, Omega_lv, *parameters))
                debug_dict.update({key:EpE_term3_arr})
            if key == 'term4':
                EpE_term4_arr.append(term(Ev,phiv, phi_primev,Omega_rv, Omega_lv, *parameters))
                debug_dict.update({key:EpE_term4_arr})
    
    debug_dict2 = {}
    op_term1_arr = []
    op_term21_arr = []
    op_term22_arr = []
    op_term2_arr = []
    for key, term in omega_phi_debug.items():
        for Ev, phiv, phiprimev, omegalv, omegamv, omegarv in zip(Hubble_arr,phi_arr,phi_prime_arr, Omega_l_arr, Omega_m_arr, Omega_r_arr):
            if key == 'term1':
                    op_term1_arr.append(term(Ev, phiv, phiprimev, omegalv, omegamv, omegarv,*parameters))
                    debug_dict2.update({key:op_term1_arr})
            if key == 'term21':
                    op_term21_arr.append(term(Ev, phiv, phiprimev, omegalv, omegamv, omegarv,*parameters))
                    debug_dict2.update({key:op_term21_arr})
            if key == 'term22':
                    op_term22_arr.append(term(Ev, phiv, phiprimev, omegalv, omegamv, omegarv,*parameters))
                    debug_dict2.update({key:op_term22_arr})
            if key == 'term2':
                    op_term2_arr.append(term(Ev, phiv, phiprimev, omegalv, omegamv, omegarv,*parameters))
                    debug_dict2.update({key:op_term2_arr})

    K_arr = []
    for Ev, phiv, phi_primev in zip(Hubble_arr, phi_arr, phi_prime_arr):
        K_arr.append(K_lambda(Ev, phiv, phi_primev, *parameters))        

    solution_arrays = {'a':a_arr_inv, 'Hubble':Hubble_arr, 'Hubble_prime':Hubble_prime_arr,'E_prime_E':E_prime_E_arr, 'scalar':phi_arr,'scalar_prime':phi_prime_arr,'scalar_primeprime':phi_primeprime_arr}
    cosmological_density_arrays = {'omega_m':Omega_m_arr,'omega_r':Omega_r_arr,'omega_l':Omega_l_arr,'omega_phi':Omega_phi_arr, 'omega_DE':Omega_DE_arr}
    cosmo_density_prime_arrays = {'omega_m_prime':Omega_m_prime_arr,'omega_r_prime':Omega_r_prime_arr,'omega_l_prime':Omega_l_prime_arr}
    force_quantities = {'A':A_arr, 'calB':calB_arr, 'calC':calC_arr, 'coupling_factor':coupling_factor_arr, 'chi_over_delta':chioverdelta_arr}
    result = {}

    for i in [solution_arrays, cosmological_density_arrays, cosmo_density_prime_arrays,force_quantities]:
        result.update(i)
    result.update({'closure_value':cl_var, 'closure_declaration':cl_declaration, 'closure_full':cl_var_full})
    result.update({'EpE_debug':debug_dict, 'omega_phi_debug':debug_dict2, 'Horndeski_K':K_arr})

    return result #a_arr_inv, Hubble_arr, E_prime_E_arr, Hubble_prime_arr, phi_prime_arr,  phi_primeprime_arr, Omega_r_arr, Omega_m_arr, Omega_DE_arr, Omega_l_arr, Omega_phi_arr, Omega_phi_diff_arr, Omega_r_prime_arr, Omega_m_prime_arr, Omega_l_prime_arr, A_arr, calB_arr, calC_arr, coupling_factor_arr, chioverdelta_arr #phi_prime_check_arr




###############################

def comp_almost_track(E_dS_fac, Omega_r0, Omega_m0, f_phi_value, almost, fried_RHS_lambda):
    E0 = comp_E_LCDM(0., Omega_r0, Omega_m0)
    EdS = E0*E_dS_fac #0.93 #*1e-30
    U0=1./EdS

    Omega_DE0 = 1. - Omega_r0 - Omega_m0
    Omega_l0 = (1.-f_phi_value)*Omega_DE0
    alpha_param = 1. - Omega_l0/EdS/EdS
    k1_dS, g31_dS = -6.*alpha_param, 2.*alpha_param
    parameters = [k1_dS, g31_dS]

    y0 = comp_param_close(fried_RHS_lambda, ['odeint_parameters',1], U0, 0.9, Omega_r0, Omega_m0, Omega_l0, parameters)
    track = (E0/EdS)**2.*y0-1.
    almost_track = track - almost
    return almost_track

def comp_E_dS_max(E_dS_max_guess, Omr, Omm, f_phi, almost, fried_RHS_lambda):
    # if f_phi == 0:
    #     if almost < 1e-3:
    #         f_phi = 1e-3
    #     else:
    #         f_phi = almost
    E_dS_max = fsolve(comp_almost_track, E_dS_max_guess, args=(Omr, Omm, f_phi, almost,fried_RHS_lambda))[0]
    return E_dS_max


### growth factor
def compute_growthfac_primes(Y, x, x_arr, Omega_m_arr, E_prime_E_arr, coupling_arr):
    """Second order differential equation for growth factor D"""
    D, Dprime = Y
    Omega_m_spl = interp1d(x_arr, Omega_m_arr, fill_value="extrapolate")
    E_prime_E_spl = interp1d(x_arr, E_prime_E_arr, fill_value="extrapolate")
    coupling_spl = interp1d(x_arr, coupling_arr, fill_value="extrapolate")
    Omega_m = Omega_m_spl(x)
    E_prime_E = E_prime_E_spl(x)
    coupling = coupling_spl(x)
    GeffOverG = 1.+coupling
    Dprimeprime = -1.*(2.+E_prime_E)*Dprime + 1.5*GeffOverG*Omega_m*D
    return [Dprime, Dprimeprime]

def solve_1stgrowth_factor(a_arr_inv, Omega_m_arr, E_prime_E_arr, coupling_arr):
    """Solve for growth mode"""
    x_arr_inv = [np.log(av) for av in a_arr_inv]
    x_arr = x_arr_inv[::-1]
    x_ini = np.min(x_arr)
    a_ini = np.exp(x_ini)
    z_ini = 1./a_ini - 1.
    #print('z_ini=', z_ini)
    D_ini_growth = 1. #np.exp(x_ini) #1.
    Dprime_ini_growth = 1. #np.exp(x_ini) #1.
    Y_ini_growth = [D_ini_growth, Dprime_ini_growth]
    #print(Y_ini_growth)
    ans = odeint(compute_growthfac_primes, Y_ini_growth, x_arr, args=(x_arr_inv, Omega_m_arr, E_prime_E_arr, coupling_arr))
    D_growth_arr = ans[:,0]
    Dprime_growth_arr = ans[:,1]
    D_growth_arr_inv = D_growth_arr[::-1]
    Dprime_growth_arr_inv = Dprime_growth_arr[::-1]
    return D_growth_arr_inv, Dprime_growth_arr_inv
