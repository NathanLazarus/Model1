#todo
#epsilon distinguishable set
#parallelize

import numpy as np
import matplotlib.pyplot as plt
from casadi import *
from scipy.optimize import fsolve
from scipy import integrate
from scipy.special import h_roots
import csv
from multiprocessing import Pool, freeze_support, cpu_count
import os
import sys
# from numpy.polynomial.hermite_e import *
import functools
import itertools
import operator
import warnings
from polyfit_helper_funcs import *
# from statsmodels.tsa.arima_process import arma_generate_sample as AR_gen



T = 30
periods = 15
k0s = np.linspace(0.1,5.1,num = 50)
z0s = np.linspace(0.8,1.2,num = 50)
degree = 6
n_iter = 20
P = 1.25
nrow_sol_path = 6
sol_path = np.array([]).reshape(nrow_sol_path,0)

ncol_reg = 5
reg_dat = np.array([]).reshape(0,ncol_reg)



alpha = 0.32
g = 1.014
delta = 0.08
beta = 0.98
lambda_zeta = 0.92
eta = 2
sigma_rho = 0.0072

lambda_tikhonov = 0.01


L_max = np.inf
L_min = 1e-6
C_min = 1e-6
K_min = -np.inf #borrowing constraint




# solve for steady state
cstar = SX.sym('cstar')
lstar = SX.sym('lstar')
kstar = SX.sym('kstar')

obj = 1


def steadystateconstraint(cstar,lstar,kstar):
    c1 = cstar - (kstar**alpha*lstar**(1-alpha) + (1-delta)*kstar - g*kstar)
    c2 = lstar - (((1-alpha)/P)*kstar**alpha*(1/cstar))**(1/(eta+alpha))
    c3 = kstar - ((g/beta - (1 - delta))*(P/alpha))**(1/(alpha-1))*lstar
    return vertcat(c1,c2,c3)


starconstraint = steadystateconstraint(cstar,lstar,kstar)
star_x_0 = DM.ones(3)

star_nlp = {'x':vertcat(cstar,lstar,kstar), 'f':obj, 'g':starconstraint}
star_solver = nlpsol('star_solver', 'ipopt', star_nlp,{'ipopt.print_level':0})
with suppress_stdout_stderr():
    star_solution = star_solver(x0=star_x_0,lbg=-1e-14,ubg=1e-14)
ssc, ssl, ssk = vertsplit(star_solution['x'])
print(ssc, ssl, ssk)


   






consumption = SX.sym('consumption', T, 1)
labor = SX.sym('labor', T, 1)
capital = SX.sym('capital', T+1, 1)


# Define the start point
x_0 = vertcat(DM.ones(T), DM.ones(T),DM.ones(T+1))


def condition1(consumption, labor, capital, zeta):
    return ((1/g)*(zeta[:T]*capital[:T]**alpha*labor[:T]**(1-alpha) + (1-delta) * capital[:T] - consumption[:T]) - capital[1:T+1])
    
def condition2(consumption, labor, capital, zeta, lfunc, cfunc):
    return (((beta/g)*((1/P)*(zeta[:T-1]**lambda_zeta) * alpha*capital[1:T]**(alpha-1)* labor[1:T]**(1-alpha) + 1 - delta)*(1/consumption[1:T]))
    -  1/consumption[:T-1])

    # (logquad(capital[:T],zeta[:T],8,lfunc,cfunc)
         # - 1/consumption[:T]) #length = T

def condition3(consumption, labor, capital, zeta):
    return ((1/P)*(1-alpha)*zeta[:T]*capital[:T]**alpha*labor[:T]**(-alpha)
        - consumption[:T]*labor[:T]**eta) #length = T

def FOCs(consumption, labor, capital, zeta, lfunc = 1, cfunc = 1):
    return vertcat(
        condition1(consumption, labor, capital, zeta),
        condition2(consumption, labor, capital, zeta, lfunc, cfunc),
        condition3(consumption, labor, capital, zeta)
    )

def dotheThing(k0, z0):

    zeta = z0 ** (lambda_zeta**DM(range(T)))



    lower_bound_C = vertcat(DM.zeros(T) + C_min)    # lower bound on the consumption -> not binding anyway
    upper_bound_C = vertcat(DM.zeros(T) + np.inf)

    lower_bound_L = vertcat(DM.zeros(T) + L_min)
    upper_bound_L = vertcat(DM.zeros(T) + L_max) # upper bound on labor also doesn't bind

    lower_bound_K = vertcat(k0, DM.zeros(T) + K_min)
    upper_bound_K = vertcat(k0, DM.zeros(T) + np.inf)

    lb_x = vertcat(lower_bound_C, lower_bound_L, lower_bound_K)
    ub_x = vertcat(upper_bound_C, upper_bound_L, upper_bound_K)
    
    objective = (capital[T] - ssk) ** 2
    nonlin_con = FOCs(consumption, labor, capital, zeta)
    nlp = {'x':vertcat(consumption, labor, capital), 'f':objective, 'g':nonlin_con}
    solver = nlpsol('solver', 'ipopt', nlp,{'ipopt.print_level':5})
    solution = solver(x0=x_0,lbx=lb_x,ubx=ub_x,lbg=-1e-10,ubg=1e-10)
    sol = solution['x']
    c_sim, l_sim, k_sim = vertsplit(sol,[0,consumption.shape[0],consumption.shape[0] + labor.shape[0], consumption.shape[0] + labor.shape[0] + capital.shape[0]])

    c = np.array([]).reshape(0,1)
    l = np.array([]).reshape(0,1)
    k = np.array([]).reshape(0,1)
    z = np.array([]).reshape(0,1)
    c = np.concatenate([c,c_sim[:periods]])
    l = np.concatenate([l,l_sim[:periods]])
    k = np.concatenate([k,k_sim[:periods]])
    z = np.concatenate([z,zeta[:periods]])
    return np.hstack([k,z,c,l])

# np.random.seed(1)
# zeta = DM(exp(AR_gen([1,-lambda_zeta],[1],T,burnin=0,scale = sigma_rho)))

grid = list(itertools.product(k0s,z0s))
    
if __name__ == "__main__":
    with suppress_stdout_stderr():
        with Pool(cpu_count()) as p:
            data = np.vstack(
                p.starmap(
                    dotheThing,
                    grid,
                )
            )

    np.savetxt(
        "NLCEQdata.csv",
        data,
        delimiter=",",
        comments="",
        header="k,z,c,l",
    )

    [k, z, c, l] = np.hsplit(data, 4)
    # print(np.vstack(data))
    # print(k,z)

    with suppress_stdout_stderr():
        c_coefs = herme2d_fit([k,z],c,degree).T
    if c_coefs.shape[0] > 1:
        c_coefs = c_coefs.T
    print(c_coefs)
    def c_function(k,z):
        c_poly = vander_nd_flat_SYM((hermevander_casadiSYM,hermevander_casadiSYM),[k,z],[degree,degree]) @ c_coefs.T
        output_length = c_poly.shape[0]
        return fmax(c_poly,SX.zeros(output_length)+C_min)

    l_coefs = herme2d_fit([k,z],l,degree).T
    if l_coefs.shape[0] > 1:
        l_coefs = l_coefs.T
    
    def l_function(k,z):
        l_poly = vander_nd_flat_SYM((hermevander_casadiSYM,hermevander_casadiSYM),[k,z],[degree,degree]) @ l_coefs.T
        output_length = l_poly.shape[0]
        return fmax(l_poly,SX.zeros(output_length)+L_min)

    print(FOCs(c_function(DM.ones(40)*ssk,DM.ones(40)),l_function(DM.ones(40)*ssk,DM.ones(40)),DM.ones(40)*ssk,DM.ones(40),l_function,c_function))