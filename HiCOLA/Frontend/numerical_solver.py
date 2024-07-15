###################
# Loading modules #
###################

import numpy as np
import scipy.integrate as integrate
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
from scipy.optimize import curve_fit
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

def make_timeout(start_time, max_time):
    def timeout(t, y, Hubble0, Omega_r_ini, Omega_m_ini, Omega_l_ini, E_prime_E_lambda, E_prime_E_safelambda, phi_primeprime_lambda, phi_primeprime_safelambda, A_lambda, cl_declaration, parameters,threshold,GR_flag):
        elapsed = time.time() - start_time
        return max_time - elapsed
    
    timeout.terminal = True
    return timeout

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

#def comp_alpha_M_propto_Omega_DE_LCDM(c_M, Omega_DE):
#    alpha_M_DE = c_M*Omega_DE
#    return alpha_M_DE

#def comp_alpha_M_prime_propto_Omega_DE_LCDM(c_M, Omega_DE_prime):
#    alpha_M_prime_DE = c_M*Omega_DE_prime
#    return alpha_M_prime_DE

#def alpha_M_int_propto_Omega_DE_LCDM(x, Omega_r0, Omega_m0, c_M):
#    Omega_DE = comp_Omega_DE_LCDM(x, Omega_r0, Omega_m0)
#    alpha_M_int = c_M * Omega_DE
#    return alpha_M_int

#def comp_Meffsq_x2_x1_propto_Omega_DE_LCDM(x1, x2, Omega_r0, Omega_m0, c_M):
#    Meffsq_x2_x1 = np.exp(integrate.quad(alpha_M_int_propto_Omega_DE_LCDM, x1, x2, args=(Omega_r0, Omega_m0, c_M))[0])
#    return Meffsq_x2_x1

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

def comp_LCDM(read_out_dict):
    [Npoints, z_max, suppression_flag, threshold, GR_flag] = read_out_dict['simulation_parameters']
    [Omega_r0, Omega_m0, Omega_l0] = read_out_dict['cosmological_parameters']

    z_final = 0.
    x_ini = np.log(1./(1.+z_max))
    x_final = np.log(1./(1.+z_final))
    x_arr = np.linspace(x_ini, x_final, Npoints)
    a_arr = [np.exp(x) for x in x_arr]
    a_arr = np.array(a_arr)
    z_arr = 1/a_arr - 1

    E = comp_E_LCDM(z_arr, Omega_r0, Omega_m0)
    E_prime_E = comp_E_prime_E_LCDM(x_arr, Omega_r0, Omega_m0)
    Omega_m = comp_Omega_m_LCDM(z_arr, Omega_r0, Omega_m0)
    Omega_r = comp_Omega_r_LCDM(z_arr, Omega_r0, Omega_m0)
    Omega_l = comp_Omega_L_LCDM(z_arr, Omega_r0, Omega_m0)
    return [E, E_prime_E, Omega_m, Omega_r, Omega_l]

#fried1 is ultimately the closure equation for the density parameters
def fried1(phi_prime, k1, g1, Omega_r, Omega_m, E, alpha_M, Ms_Mp, Meffsq_Mpsq):
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
            return fried_RHS_lambda(cl_variable,phi,phi_prime,Omega_r,Omega_m,Omega_l, *parameters) #Closure used to compute E0
        if cl_declaration[1] == 1:
            return fried_RHS_lambda(E,cl_variable,phi_prime,Omega_r,Omega_m,Omega_l,*parameters) #Closure used to compute phi0
        if cl_declaration[1] == 2:
            return fried_RHS_lambda(E, phi, cl_variable, Omega_r,Omega_m,Omega_l,*parameters) #Closure used to compute phi_prime0
        if cl_declaration[1] == 3:
            return fried_RHS_lambda(E,phi,phi_prime,cl_variable,Omega_m,Omega_l,*parameters) #Closure used to compute Omega_r0
        if cl_declaration[1] == 4:
            return fried_RHS_lambda(E,phi,phi_prime,Omega_r,cl_variable,Omega_l,*parameters) #Closure used to compute Omega_m0
    if cl_declaration[0] == 'parameters':
        parameters[cl_declaration[1]] = cl_variable
        return fried_RHS_lambda(E,phi,phi_prime,Omega_r,Omega_m,Omega_l,*parameters)
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

    cl_variable0 = cl_variable[0][0] #to unpack the fsolve solution/iteration

    return cl_variable0, cl_variable

def comp_E_closure(fried_closure_lambda, cl_declaration, E0, phi, phi_prime, Omega_r, Omega_m, Omega_l, a_arr, parameters): #calling the closure-fixed parameter "k1" is arbitrary, the choice of which parameter to fix is determined by lambdification or fried_RHS
    cl_dec = cl_declaration.copy()
    cl_dec[1] = 0
    E_cl = np.zeros(len(a_arr))

    for a in a_arr:
        full = fsolve(fried_RHS_wrapper, E0, args=(cl_dec, fried_closure_lambda, E0, phi[a_arr==a], phi_prime[a_arr==a], Omega_r[a_arr==a],Omega_m[a_arr==a], Omega_l[a_arr==a],parameters), xtol=1e-6,full_output=True) #make sure arguments match lambdification line in run_builder.py
        #if full[2] != 1: print(full[3])
        E_cl[a_arr==a] = full[0][0]
        E0 = E_cl[a_arr==a]

    return E_cl

def comp_primes(x, Y, E0, Omega_r0, Omega_m0, Omega_l0, E_prime_E_lambda, E_prime_E_safelambda, phi_primeprime_lambda, phi_primeprime_safelambda, A_lambda, cl_declaration, parameters,threshold=1e-3,GR_flag=False): #x, Y swapped for solve_ivp ###ADD LAMBDA FUNCTION AS ARGUMENT###


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


def run_solver(read_out_dict):


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

    parameter_symbols = read_out_dict['symbol_list']
    odeint_parameter_symbols = read_out_dict['odeint_parameter_symbols']
    cl_declaration = read_out_dict['closure_declaration']

    z_final = 0.
    x_ini = np.log(1./(1.+z_max))
    x_final = np.log(1./(1.+z_final))
    x_arr = np.linspace(x_ini, x_final, Npoints)
    a_arr = [np.exp(x) for x in x_arr]
    if read_out_dict['forwards_flag'] is True:
        x_arr_inv = x_arr
        a_arr_inv = a_arr
        x_ini = np.log(1./(1.+z_final))
        x_final = np.log(1./(1.+z_max))

        #omega_i(z=z_max) for forward solver ICs
        Omega_l_ini = Omega_l0/Hubble0/Hubble0
        Omega_m_ini = Omega_m0*((1+z_max)**3.)/Hubble0/Hubble0
        Omega_r_ini = Omega_r0*((1+z_max)**4)/Hubble0/Hubble0
    else:
        x_arr_inv = x_arr[::-1]
        a_arr_inv = a_arr[::-1]

        Omega_l_ini = Omega_l0
        Omega_m_ini = Omega_m0
        Omega_r_ini = Omega_r0

    if GR_flag is True:
        phi_prime0 = 0.
        
    cl_var, cl_var_full = comp_param_close(fried_RHS_lambda, cl_declaration, Hubble0, phi0, phi_prime0, Omega_r_ini, Omega_m_ini, Omega_l_ini, parameters)

    #print(cl_var_full)

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

    max_time = 10.0
    start_time = time.time()
    timeout = make_timeout(start_time, max_time)

    if suppression_flag is True:
        with stdout_redirected():
            ans = solve_ivp(comp_primes,[x_final, x_ini], Y0, t_eval=x_arr_inv, method='Radau', args=(Hubble0, Omega_r_ini, Omega_m_ini, Omega_l_ini, E_prime_E_lambda, E_prime_E_safelambda, phi_primeprime_lambda, phi_primeprime_safelambda, A_lambda, cl_declaration, parameters,threshold,GR_flag), rtol = 1e-10, events=timeout)#, hmax=hmaxv) #k1=-6, g1 = 2
    else:
        ans = solve_ivp(comp_primes,[x_final,x_ini], Y0, t_eval=x_arr_inv, method='Radau', args=(Hubble0, Omega_r_ini, Omega_m_ini, Omega_l_ini, E_prime_E_lambda, E_prime_E_safelambda, phi_primeprime_lambda, phi_primeprime_safelambda, A_lambda, cl_declaration, parameters,threshold,GR_flag), rtol = 1e-10, events=timeout)#, hmax=hmaxv)

    ans = ans["y"].T
    phi_arr = ans[:,0]
    phi_prime_arr = ans[:,1]
    Hubble_arr = ans[:,2]
    Omega_r_arr = ans[:,3]
    Omega_m_arr = ans[:,4]
    Omega_l_arr = ans[:,5]

    #simple check for numerical discontinuity which returns False if one is found
    if len(Hubble_arr) != 1000:
        return False

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

    E_arr = Hubble_arr/Hubble0 #if E_ini is normalised to 1
    chioverdelta_arr = chi_over_delta(a_arr_inv, E_arr, calB_arr, calC_arr, Omega_m0)

    solution_arrays = {'a':a_arr_inv, 'Hubble':Hubble_arr, 'Hubble_prime':Hubble_prime_arr,'E_prime_E':E_prime_E_arr, 'scalar':phi_arr,'scalar_prime':phi_prime_arr,'scalar_primeprime':phi_primeprime_arr}
    cosmological_density_arrays = {'omega_m':Omega_m_arr,'omega_r':Omega_r_arr,'omega_l':Omega_l_arr,'omega_phi':Omega_phi_arr, 'omega_DE':Omega_DE_arr}
    cosmo_density_prime_arrays = {'omega_m_prime':Omega_m_prime_arr,'omega_r_prime':Omega_r_prime_arr,'omega_l_prime':Omega_l_prime_arr}
    force_quantities = {'A':A_arr, 'calB':calB_arr, 'calC':calC_arr, 'coupling_factor':coupling_factor_arr, 'chi_over_delta':chioverdelta_arr}
    result = {}

    for i in [solution_arrays, cosmological_density_arrays, cosmo_density_prime_arrays,force_quantities]:
        result.update(i)
    result.update({'closure_value':cl_var, 'closure_declaration':cl_declaration})

    return result

def comp_alphas(read_out_dict, background_quantities):
    M_star_sqrd_lamb = read_out_dict['M_star_sqrd_lambda']
    alpha_M_lamb = read_out_dict['alpha_M_lambda']
    alpha_B_lamb = read_out_dict['alpha_B_lambda']
    alpha_K_lamb = read_out_dict['alpha_K_lambda']
    parameters = read_out_dict['Horndeski_parameters']

    E = background_quantities['Hubble']
    phi = background_quantities['scalar']
    phi_prime = background_quantities['scalar_prime']
    
    M_star_sqrd_evaluated = M_star_sqrd_lamb(phi, *parameters)
    if not isinstance(M_star_sqrd_evaluated, np.ndarray):
        M_star_sqrd_evaluated = np.ones(len(E))*M_star_sqrd_evaluated
    alpha_M_evaluated = alpha_M_lamb(E, phi, phi_prime, *parameters)
    if not isinstance(alpha_M_evaluated, np.ndarray):
        alpha_M_evaluated = np.ones(len(E))*alpha_M_evaluated
    alpha_B_evaluated = alpha_B_lamb(E, phi, phi_prime, *parameters)
    if not isinstance(alpha_B_evaluated, np.ndarray):
        alpha_B_evaluated = np.ones(len(E))*alpha_B_evaluated
    alpha_K_evaluated = alpha_K_lamb(E, phi, phi_prime, *parameters)
    if not isinstance(alpha_K_evaluated, np.ndarray):
        alpha_K_evaluated = np.ones(len(E))*alpha_K_evaluated
    alphas = {'M_star_sq':M_star_sqrd_evaluated, 'alpha_M':alpha_M_evaluated, 'alpha_B':alpha_B_evaluated,'alpha_K':alpha_K_evaluated}
    return alphas

#based on already calculated background evolution only
def comp_w_phi(read_out_dict, background_quantities):
    """
    Computes equation of state for scalar field using the Friedman equations. 
    Code equivalent of 3.5 in https://iopscience.iop.org/article/10.1088/1475-7516/2014/07/050.
    """
    E = background_quantities['Hubble']
    E_prime = background_quantities['Hubble_prime']
    Omega_m = background_quantities['omega_m']
    Omega_r = background_quantities['omega_r']
    Omega_l = background_quantities['omega_l']
    M_star_sq = background_quantities['M_star_sq']

    H_0 = 100*read_out_dict['little_h']

    P_phi = -2*H_0**2*E*E_prime - 3*H_0**2*E**2*(1+(Omega_r*1/3-Omega_l)/M_star_sq)
    rho_phi = 3*H_0**2*E**2*(1-(Omega_m+Omega_r+Omega_l)/M_star_sq)
    w_phi = P_phi/rho_phi
    return w_phi, P_phi, rho_phi

#based on P and rho funcs of background evolution
def comp_w_phi2(read_out_dict, background_quantities):
    """
    Computes equation of state for scalar field using expressions for rho and P in terms of Horndeski functions.
    """
    P_phi_lamb = read_out_dict['P_DE_lambda']
    rho_phi_lamb = read_out_dict['rho_DE_lambda']
    parameters = read_out_dict['Horndeski_parameters']

    E = background_quantities['Hubble']
    E_prime = background_quantities['Hubble_prime']
    phi = background_quantities['scalar']
    phi_prime = background_quantities['scalar_prime']
    phi_primeprime = background_quantities['scalar_primeprime']

    P_phi_evaluated = P_phi_lamb(E, E_prime, phi, phi_prime, phi_primeprime, *parameters)
    rho_phi_evaluated = rho_phi_lamb(E, phi, phi_prime, *parameters)
    w_phi = P_phi_evaluated/rho_phi_evaluated
    return w_phi, P_phi_evaluated, rho_phi_evaluated

def comp_w_eff(background_quantities):
    """
    Computes effective equation of state.
    """
    EprimeE = background_quantities['E_prime_E']
    w_eff = -1-EprimeE*2/3
    return w_eff

def comp_stability(read_out_dict, background_quantities):
    """
    Computes Q_s and c_s^2*D and checks stability conditions.
    """
    parameters = read_out_dict['Horndeski_parameters']
    Q_s_lamb = read_out_dict['Q_s_lambda']
    c_s_sq_lamb = read_out_dict['c_s_sq_lambda']

    a = background_quantities['a']
    E = background_quantities['Hubble']
    Eprime = background_quantities['Hubble_prime']
    phi = background_quantities['scalar']
    phi_prime = background_quantities['scalar_prime']
    Omega_m = background_quantities['omega_m']
    Omega_r = background_quantities['omega_r']
    alpha_B_evaluated = background_quantities['alpha_B']

    lna = np.log(a)
    alphaBprime = np.gradient(alpha_B_evaluated, lna)

    Q_s_evaluated = Q_s_lamb(E, phi, phi_prime, *parameters)
    c_s_sq_evaluated = c_s_sq_lamb(E, Eprime, phi, phi_prime, Omega_m, Omega_r, alphaBprime, *parameters)

    if (Q_s_evaluated>0).all() and (c_s_sq_evaluated>0).all():
        unstable = False
        #print('Stability conditions satisified')
    elif (Q_s_evaluated>0).all():
        unstable = 2
        #print('Warning: Stability condition not satisfied: c_s_sq not always > 0')
    elif (c_s_sq_evaluated>0).all():
        unstable = 3
        #print('Warning: Stability condition not satisfied: Q_s not always > 0')
    else:
        unstable = 1
        #print('Warning: Stability conditions not satisfied: Q_s and c_s_sq not always > 0')
    
    return Q_s_evaluated, c_s_sq_evaluated, unstable
    
def alpha_X1(a, alpha_X0, Omega_m0, Omega_r0):
    """
    Function used for the first parameterisation of property functions.
    """
    z = 1/a - 1

    Omega_l = comp_Omega_L_LCDM(z, Omega_r0, Omega_m0)

    alph_X = alpha_X0*Omega_l/(1 - Omega_m0 - Omega_r0)
    return alph_X

def alpha_X2(a, alpha_M0, alpha_B0, alpha_K0, q):
    """
    Function used for the second parameterisation of property functions.
    """
    alph_M = alpha_M0*a**q
    alph_B = alpha_B0*a**q
    alph_K = alpha_K0*a**q
    return np.concatenate((alph_M, alph_B, alph_K))

def alpha_X3(a,  alpha_X0, q_X):
    """
    Function used for the third parameterisation of property functions.
    """
    alph_X = alpha_X0*a**q_X
    return alph_X

def r_chi2(y, y_model):
    """
    This function calculates the value of reduced chi squared between a given model and a set of data.
    """
    if (y_model==0).all():
        return 0
    else:
        chisq = np.sum(((y-y_model)**2.0)/(y_model**2.0))
        return chisq/len(y)

def parameterise1(a_arr, z_max, Omega_m0, Omega_r0, alpha_M_arr, alpha_B_arr, alpha_K_arr):
    """
    Finds best fit parameters for the first parameterisation and evaluates goodness of fit using Bayesian information criterion (BIC).
    """
    z_arr = 1/a_arr - 1
    a_fit = a_arr[z_arr<z_max]
    alpha_M_fit = alpha_M_arr[z_arr<z_max]
    alpha_B_fit = alpha_B_arr[z_arr<z_max]
    alpha_K_fit = alpha_K_arr[z_arr<z_max]

    lower = [-50,Omega_m0-1e-10,Omega_r0-1e-10]
    upper = [50,Omega_m0+1e-10,Omega_r0+1e-10]
    (alpha_M01, junk1, junk2), junk = curve_fit(alpha_X1, a_fit, alpha_M_fit, bounds=(lower,upper))
    (alpha_B01, junk1, junk2), junk = curve_fit(alpha_X1, a_fit, alpha_B_fit, bounds=(lower,upper))
    (alpha_K01, junk1, junk2), junk = curve_fit(alpha_X1, a_fit, alpha_K_fit, bounds=(lower,upper))

    alpha_M_param1_arr = alpha_X1(a_arr, alpha_M01, Omega_m0, Omega_r0)
    alpha_B_param1_arr = alpha_X1(a_arr, alpha_B01, Omega_m0, Omega_r0)
    alpha_K_param1_arr = alpha_X1(a_arr, alpha_K01, Omega_m0, Omega_r0)

    BIC_M = r_chi2(alpha_M_param1_arr[z_arr<z_max], alpha_M_fit) + np.log(len(a_fit))/len(a_fit)
    BIC_B = r_chi2(alpha_B_param1_arr[z_arr<z_max], alpha_B_fit) + np.log(len(a_fit))/len(a_fit)
    BIC_K = r_chi2(alpha_K_param1_arr[z_arr<z_max], alpha_K_fit) + np.log(len(a_fit))/len(a_fit)

    model_vals = (alpha_M_param1_arr, alpha_B_param1_arr, alpha_K_param1_arr)
    model_params = (alpha_M01, alpha_B01, alpha_K01)
    fit_goodness = (BIC_M, BIC_B, BIC_K)
    return model_vals, model_params, fit_goodness

def parameterise2(a_arr, z_max, alpha_M_arr, alpha_B_arr, alpha_K_arr):
    """
    Finds best fit parameters for the second parameterisation and evaluates goodness of fit using Bayesian information criterion (BIC).
    """
    z_arr = 1/a_arr - 1
    a_fit = a_arr[z_arr<z_max]
    alpha_M_fit = alpha_M_arr[z_arr<z_max]
    alpha_B_fit = alpha_B_arr[z_arr<z_max]
    alpha_K_fit = alpha_K_arr[z_arr<z_max]

    alpha_X2_fit = np.concatenate((alpha_M_fit, alpha_B_fit, alpha_K_fit))
    (alpha_M02, alpha_B02, alpha_K02, q), junk = curve_fit(alpha_X2, a_fit, alpha_X2_fit, bounds=([-50,-50,-50,0],[50,50,50,6]))

    alpha_M_param2_arr = alpha_X3(a_arr, alpha_M02, q)
    alpha_B_param2_arr = alpha_X3(a_arr, alpha_B02, q)
    alpha_K_param2_arr = alpha_X3(a_arr, alpha_K02, q)

    BIC_M = r_chi2(alpha_M_param2_arr[z_arr<z_max], alpha_M_fit) + np.log(len(a_fit))/len(a_fit)
    BIC_B = r_chi2(alpha_B_param2_arr[z_arr<z_max], alpha_B_fit) + np.log(len(a_fit))/len(a_fit)
    BIC_K = r_chi2(alpha_K_param2_arr[z_arr<z_max], alpha_K_fit) + np.log(len(a_fit))/len(a_fit)

    model_vals = (alpha_M_param2_arr, alpha_B_param2_arr, alpha_K_param2_arr)
    model_params = (alpha_M02, alpha_B02, alpha_K02, q)
    fit_goodness = (BIC_M, BIC_B, BIC_K)
    return model_vals, model_params, fit_goodness

def parameterise3(a_arr, z_max, alpha_M_arr, alpha_B_arr, alpha_K_arr):
    """
    Finds best fit parameters for the third parameterisation and evaluates goodness of fit using Bayesian information criterion (BIC).
    """
    z_arr = 1/a_arr - 1
    a_fit = a_arr[z_arr<z_max]
    alpha_M_fit = alpha_M_arr[z_arr<z_max]
    alpha_B_fit = alpha_B_arr[z_arr<z_max]
    alpha_K_fit = alpha_K_arr[z_arr<z_max]

    (alpha_M03, q_M), junk = curve_fit(alpha_X3, a_fit, alpha_M_fit, bounds=([-50,0],[50,6]))
    (alpha_B03, q_B), junk = curve_fit(alpha_X3, a_fit, alpha_B_fit, bounds=([-50,0],[50,6]))
    (alpha_K03, q_K), junk = curve_fit(alpha_X3, a_fit, alpha_K_fit, bounds=([-50,0],[50,6]))

    alpha_M_param3_arr = alpha_X3(a_arr, alpha_M03, q_M)
    alpha_B_param3_arr = alpha_X3(a_arr, alpha_B03, q_B)
    alpha_K_param3_arr = alpha_X3(a_arr, alpha_K03, q_K)

    BIC_M = r_chi2(alpha_M_param3_arr[z_arr<z_max], alpha_M_fit) + np.log(len(a_fit))/len(a_fit)
    BIC_B = r_chi2(alpha_B_param3_arr[z_arr<z_max], alpha_B_fit) + np.log(len(a_fit))/len(a_fit)
    BIC_K = r_chi2(alpha_K_param3_arr[z_arr<z_max], alpha_K_fit) + np.log(len(a_fit))/len(a_fit)

    model_vals = (alpha_M_param3_arr, alpha_B_param3_arr, alpha_K_param3_arr)
    model_params = (alpha_M03, alpha_B03, alpha_K03, q_M, q_B, q_K)
    fit_goodness = (BIC_M, BIC_B, BIC_K)
    return model_vals, model_params, fit_goodness