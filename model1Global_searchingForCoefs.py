#todo
#switch from LS to LAD
#epsilon distinguishable set
#parallelize
import sys
print(sys.version)
import platform
print(platform.architecture())
import numpy as np
import matplotlib.pyplot as plt
from casadi import *
from scipy.optimize import fsolve
from scipy import integrate
from scipy.special import h_roots
import csv
from multiprocessing import Pool
import os
from numpy.polynomial.hermite_e import *
import functools
import operator
import warnings
from statsmodels.tsa.arima_process import arma_generate_sample as AR_gen

# print(GlobalOptions())
class RankWarning(UserWarning):
    """Issued by chebfit when the design matrix is rank deficient."""
    pass



# periods = 20
k0 = DM(3)
# k0s = np.linspace(0.1,5.1,num = 20)
# z0s = np.linspace(0.8,1.2,num = 41)
degree = 4
# n_iter = 15
# print(k0s)
P = 1.3206
# nrow_sol_path = 6
# sol_path = np.array([]).reshape(nrow_sol_path,0)

# ncol_reg = 5
# reg_dat = np.array([]).reshape(0,ncol_reg)


T = 100

alpha = 0.32
g = 1.014
delta = 0.08
beta = 0.98
lambda_zeta = 0.92
eta = 2
sigma_rho = 0.0072

# lambda_tikhonov = 0.01


L_max = 1
L_min = 1e-6
C_min = 1e-6




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
star_solver = nlpsol('star_solver', 'snopt', star_nlp,{'ipopt.print_level':0})
star_solution = star_solver(x0=star_x_0,lbg=-1e-14,ubg=1e-14)
ssc, ssl, ssk = vertsplit(star_solution['x'])
print(ssc, ssl, ssk)





# def l_function(k,z,l_coefs):
    # return 0.7933207-0.0271536*k+0.1006549*z
# def c_function(k,z,c_coefs):
    # return -0.1570917+0.1948116*k+0.6911522*z

def transformede(rho,capital,zeta,lfunc, lcoefs, cfunc, ccoefs):
    rho_length = rho.shape[0]
    zeta_length = zeta.shape[0]
    zeta_reshaped = reshape(zeta,zeta_length,1)
    capital_reshaped = reshape(capital,zeta_length,1)
    return ((beta/g)*((1/P)*
        (zeta_reshaped**lambda_zeta @ rho.reshape(1,rho_length)**sigma_rho) @ 
        alpha*capital_reshaped**(alpha-1)*
        (lfunc(capital_reshaped,(zeta_reshaped**lambda_zeta @ rho.reshape(1,rho_length)**sigma_rho),lcoefs))**(1-alpha) +
        1 - delta)*
        (1/(cfunc(capital_reshaped,(zeta_reshaped**lambda_zeta @ rho.reshape(1,rho_length)**sigma_rho),ccoefs))))


def logquad(capital,zeta,n,lfunc, lcoefs, cfunc, ccoefs):
    samplepoints,sampleweights = h_roots(n)
    return (transformede(np.exp(np.sqrt(2)*samplepoints),capital,zeta,lfunc, lcoefs, cfunc, ccoefs)/np.sqrt(np.pi))@sampleweights

def _nth_slice(i, ndim):
    sl = [np.newaxis] * ndim
    sl[i] = slice(None)
    return tuple(sl)


def _vander_nd(vander_fs, points, degrees):
    
    n_dims = len(vander_fs)
    if n_dims != len(points):
        raise ValueError(
            f"Expected {n_dims} dimensions of sample points, got {len(points)}")
    if n_dims != len(degrees):
        raise ValueError(
            f"Expected {n_dims} dimensions of degrees, got {len(degrees)}")
    if n_dims == 0:
        raise ValueError("Unable to guess a dtype or shape when no points are given")

    # convert to the same shape and type
    # points = tuple(np.array(tuple(points), copy=False) + 0.0)

    # produce the vandermonde matrix for each dimension, placing the last
    # axis of each in an independent trailing axis of the output
    vander_arrays = (
        vander_fs[i](points[i], degrees[i])[(...,) + _nth_slice(i, n_dims)]
        for i in range(n_dims)
    )

    # we checked this wasn't empty already, so no `initial` needed
    return functools.reduce(operator.mul, vander_arrays)


def _vander_nd_flat(vander_fs, points, degrees):
    """
    Like `_vander_nd`, but flattens the last ``len(degrees)`` axes into a single axis
    Used to implement the public ``<type>vander<n>d`` functions.
    """
    v = _vander_nd(vander_fs, points, degrees)
    v = v.reshape(v.shape[:-len(degrees)] + (-1,))
    if len(v.shape)>2.5:
        v = v.reshape(v.shape[0],v.shape[2])
    upperleft = np.zeros(((degrees[0]+1),(degrees[0]+1)),dtype=bool)
    for i in range(degrees[0]+1):
        upperleft[i,0:((degrees[0]+1)-i)]=True
    return v[0:v.shape[0],np.arange((degrees[0]+1)**2).reshape((degrees[0]+1),(degrees[0]+1))[upperleft]]
    


def hermevander_casadiSYM(x, deg):
    
    ideg = operator.index(deg)
    dims = (ideg + 1,) + x.shape
    v = SX.zeros(dims[:2])
    v[0,0:dims[1]] = 1
    if ideg > 0:
        v[1,0:dims[1]] = x.T
        for i in range(2, ideg + 1):
            v[i,0:dims[1]] = (v[i-1,0:dims[1]]*x.T - v[i-2,0:dims[1]]*(i - 1))
    return v.T

def _vander_nd_flat_SYM(vander_fs, points, degrees):
    
    n_dims = len(vander_fs)
    v = SX.zeros(points[0].shape[0],(degrees[0]+1)**2)
    for i in range(points[0].shape[0]):
        v[i,0:(degrees[0]+1)**2] = reshape((vander_fs[0](points[0][i],degrees[0]).T@vander_fs[1](points[1][i],degrees[1])).T,1,(degrees[0]+1)**2)
    upperleft = np.zeros(((degrees[0]+1),(degrees[0]+1)),dtype=bool)
    for i in range(degrees[0]+1):
        upperleft[i,0:((degrees[0]+1)-i)]=True
    return v[0:v.shape[0],np.arange((degrees[0]+1)**2).reshape((degrees[0]+1),(degrees[0]+1))[upperleft]]
   


   

# capital = SX.sym('capital', T+1, 1)


# lower_bound_C = vertcat(DM.zeros(T) + C_min)    # lower bound on the consumption -> not binding anyway
# lower_bound_L = vertcat(DM.zeros(T) + L_min)

# upper_bound_C = vertcat(DM.zeros(T) + np.inf)
# upper_bound_L = vertcat(DM.zeros(T) + L_max) # upper bound on labor also doesn't bind



  
# for iteration in range(n_iter):

#     c = np.array([]).reshape(0,1)
#     l = np.array([]).reshape(0,1)
#     k = np.array([]).reshape(0,1)
#     z = np.array([]).reshape(0,1)

#     for k0 in k0s:
np.random.seed(10)
rho = np.random.randn(T+1)*sigma_rho
toZeta = [1]
for i in range(T+1):
    toZeta.append(toZeta[i]**lambda_zeta*np.exp(rho[i]))
zeta = DM(np.array(toZeta[1:T+2]))
# zeta = exp(AR_gen([1,-0.9],[1],T,burnin=0,scale = sigma_rho))



def condition1(consumption, labour, capital):
    temp = ((1/g)*(zeta[:T]*capital[:T]**alpha*labour[:T]**(1-alpha) + (1-delta) * capital[:T] - consumption[:T]) - capital[1:T+1])
    return temp
    
def condition2(consumption, labour, capital, lfunc, lcoefs, cfunc, ccoefs):
    temp = (logquad(capital[:T], zeta[:T], 8, lfunc, lcoefs, cfunc, ccoefs)
         - 1/consumption[:T])
    return temp #length = T

def condition3(consumption, labour, capital):
    temp = ((1/P)*(1-alpha)*zeta[:T]*capital[:T]**alpha*labour[:T]**(-alpha)
        - consumption[:T]*labour[:T]**eta)
    return temp #length = T

def FOCs(capital, lcoefs, ccoefs):
    consumption = c_function(reshape(capital[:T],T,1), reshape(zeta[:T],T,1), ccoefs)
    labour = l_function(reshape(capital[:T],T,1), reshape(zeta[:T],T,1), lcoefs)
    return vertcat(condition1(consumption, labour, capital),condition2(consumption, labour, capital, l_function, lcoefs, c_function, ccoefs),
        condition3(consumption, labour, capital))

def c_function(k,z,ccoefs):
    c_poly = _vander_nd_flat_SYM((hermevander_casadiSYM,hermevander_casadiSYM),[k,z],[degree,degree]) @ ccoefs
    output_length = c_poly.shape[0]
    return fmax(c_poly,SX.zeros(output_length)+C_min)

def l_function(k,z,lcoefs):
    l_poly = _vander_nd_flat_SYM((hermevander_casadiSYM,hermevander_casadiSYM),[k,z],[degree,degree]) @ lcoefs
    output_length = l_poly.shape[0]
    return fmax(l_poly,SX.zeros(output_length)+L_min)


capital_as_a_consequence_of_l_and_c = SX.sym('capital_as_a_consequence_of_l_and_c',T+1,1)
lower_bound_K = vertcat(k0, -DM.ones(T)*np.inf)
upper_bound_K = vertcat(k0, DM.ones(T)*np.inf)
l_coefs, c_coefs = [DM([0.577233, 0.473765, -0.124375, 1.22109, 0.119698, 0.849019, 0.649891, 0.26678, 0.295526, 0.80309, 0.394419, -0.865858, -0.318923, -1.05689, 0.496659]),
    DM([0.656539, 0.206159, -0.706945, -0.943086, -0.72148, 0.151223, -1.77132, -3.8663, -1.06527, -1.00243, -2.4653, -3.82312, -2.46326, 5.6092, -1.05324])]

objective = 1
nonlin_con = condition1(c_function(reshape(capital_as_a_consequence_of_l_and_c[:T],T,1),reshape(zeta[:T],T,1),c_coefs),
                        l_function(reshape(capital_as_a_consequence_of_l_and_c[:T],T,1),reshape(zeta[:T],T,1),l_coefs),
                        capital_as_a_consequence_of_l_and_c)
nlp = {'x':capital_as_a_consequence_of_l_and_c, 'f':objective, 'g':nonlin_con}
solver = nlpsol('solver', 'ipopt', nlp,{'ipopt.print_level':5})
solution = solver(x0=DM.ones(T+1), lbx=lower_bound_K, ubx=upper_bound_K, lbg=-1e-12, ubg=1e-12)

print(solution['x'])
print(FOCs(solution['x'], l_coefs, c_coefs))

l_coefs = SX.sym('l_coefs', (degree + 1) * (degree + 2) // 2, 1)
c_coefs = SX.sym('c_coefs', (degree + 1) * (degree + 2) // 2, 1)
capital = SX.sym('capital',T+1,1)
# Define the start point
x_0 = vertcat(DM([0.222395, 0.224298, -0.58592, -0.109818, -0.10865, 0.0214373, 0.0216639, 0.0201385, -0.0100878, -0.00776822, -0.00785973, 0.0208302, -0.00393461, -0.00453676, 0.00133056]),
    DM([0.0504284, 0.0493043, -0.162497, -0.0230653, -0.0214815, 0.043345, 0.0435092, -0.0181494, -0.0207192, 0.0255858, 0.0250109, -0.051569, 0.00956733, 0.00855651, -0.000673861]),
    DM.ones(T+1))

# lower_bound_K = vertcat(k0, 1e-9 * DM.ones(T-1), ssk/2)
# upper_bound_K = vertcat(k0, DM.ones(T-1)*np.inf, ssk/2)
lower_bound_K = vertcat(k0, -DM.ones(T-1)*np.inf, ssk/2)
upper_bound_K = vertcat(k0, DM.ones(T-1)*np.inf, ssk/2)
lb_x = vertcat(DM.ones(c_coefs.shape[0])*-np.inf, DM.ones(l_coefs.shape[0])*-np.inf, lower_bound_K)
ub_x = vertcat(DM.ones(c_coefs.shape[0])*np.inf, DM.ones(l_coefs.shape[0])*np.inf, upper_bound_K)
objective = sum1(vertcat(l_coefs, c_coefs)**2)
nonlin_con = FOCs(capital, l_coefs, c_coefs)
nlp = {'x':vertcat(l_coefs, c_coefs, capital), 'f':objective, 'g':nonlin_con}
solver = nlpsol('solver', 'ipopt', nlp, {'ipopt.print_level':5,'ipopt.max_iter':1000})
solution = solver(x0 = x_0,lbx = lb_x,ubx = ub_x,lbg = -1e-12,ubg = 1e-12)

sol = solution['x']
print(sol)
# print(FOCs(sol[:T+1],sol[T+1:2*T+1],sol[2*T+1:]))
# c_sim, l_sim, k_sim = sol[:T],sol[T:2*T],sol[2*T:]
# # r = (1/P)*(alpha)*k[:T]**(alpha-1)*l**(1-alpha)-delta

# c = np.concatenate([c,c_sim[:periods]])
# l = np.concatenate([l,l_sim[:periods]])
# k = np.concatenate([k,k_sim[:periods]])
# z = np.concatenate([z,zeta[:periods]])

# print(l_coefs)
# print(c_coefs)
# print(l_function(DM.ones(5)*ssk,DM.ones(5),l_coefs))
# print(c_function(DM.ones(2)*ssk,DM.ones(2),c_coefs))
