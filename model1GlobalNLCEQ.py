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
from numpy.polynomial.hermite_e import *
import functools
import operator
import warnings
from statsmodels.tsa.arima_process import arma_generate_sample as AR_gen

class RankWarning(UserWarning):
    """Issued by chebfit when the design matrix is rank deficient."""
    pass



periods = 5
k0s = np.linspace(0.1,5.1,num = 20)
z0s = np.linspace(0.8,1.2,num = 10)
degree = 3
n_iter = 20
P = 1.25
nrow_sol_path = 6
sol_path = np.array([]).reshape(nrow_sol_path,0)

ncol_reg = 5
reg_dat = np.array([]).reshape(0,ncol_reg)


T = 20

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
star_solution = star_solver(x0=star_x_0,lbg=-1e-14,ubg=1e-14)
ssc, ssl, ssk = vertsplit(star_solution['x'])
print(ssl,ssc,ssk)





# def l_function(k,z):
#         return 0.7933207-0.0271536*k+0.1006549*z
# def c_function(k,z):
#     return -0.1570917+0.1948116*k+0.6911522*z

# def transformede(rho,capital,zeta,lfunc,cfunc):
#     rho_length = rho.shape[0]
#     zeta_length = zeta.shape[0]
#     return ((beta/g)*((1/P)*
#         (reshape(zeta,zeta_length,1)**lambda_zeta @ rho.reshape(1,rho_length)**sigma_rho) @ 
#         alpha*reshape(capital,zeta_length,1)**(alpha-1)*
#         (lfunc(reshape(capital,zeta_length,1),(reshape(zeta,zeta_length,1)**lambda_zeta @ rho.reshape(1,rho_length)**sigma_rho)))**(1-alpha) +
#         1 - delta)*
#         (1/(cfunc(reshape(capital,zeta_length,1),(reshape(zeta,zeta_length,1)**lambda_zeta @ rho.reshape(1,rho_length)**sigma_rho)))))


# def logquad(capital,zeta,n,lfunc,cfunc):
#     samplepoints,sampleweights = h_roots(n)
#     return (transformede(np.exp(np.sqrt(2)*samplepoints),capital,zeta,lfunc,cfunc)/np.sqrt(np.pi))@sampleweights

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
    

def herme2d_fit(x, y, deg, rcond=None, full=False, w=None):
    """
    Helper function used to implement the ``<type>fit`` functions.
    Parameters
    ----------
    vander_f : function(array_like, int) -> ndarray
        The 1d vander function, such as ``polyvander``
    c1, c2 :
        See the ``<type>fit`` functions for more detail
    """
    # x = np.asarray(x) + 0.0
    y = np.asarray(y) + 0.0
    deg = np.asarray(deg)

    # check arguments.
    if deg.ndim > 1 or deg.dtype.kind not in 'iu' or deg.size == 0:
        raise TypeError("deg must be an int or non-empty 1-D array of int")
    if deg.min() < 0:
        raise ValueError("expected deg >= 0")
    # if x.ndim != 1:
        # raise TypeError("expected 1D vector for x")
    # if x.size == 0:
    #     raise TypeError("expected non-empty vector for x")
    if y.ndim < 1 or y.ndim > 2:
        raise TypeError("expected 1D or 2D array for y")
    # if len(x) != len(y):
        # raise TypeError("expected x and y to have same length")

    # if deg.ndim == 0:
    #     lmax = deg
    #     order = lmax + 1
    #     van = vander_f(x, lmax)

    # else:
    #     deg = np.sort(deg)
    #     lmax = deg[-1]
    #     order = len(deg)
    #     van = vander_f(x, lmax)[:, deg]
    van = _vander_nd_flat((hermevander,hermevander),x,[deg,deg])
    order = van.shape[1] #(deg+1)*(deg+2)/2


    # set up the least squares matrices in transposed form
    lhs = van.T
    rhs = y.T
    if w is not None:
        w = np.asarray(w) + 0.0
        if w.ndim != 1:
            raise TypeError("expected 1D vector for w")
        if len(x) != len(w):
            raise TypeError("expected x and w to have same length")
        # apply weights. Don't use inplace operations as they
        # can cause problems with NA.
        lhs = lhs * w
        rhs = rhs * w

    # set rcond
    # if rcond is None:
    #     rcond = len(y)*np.finfo(y.dtype).eps

    # Determine the norms of the design matrix columns.
    if issubclass(lhs.dtype.type, np.complexfloating):
        scl = np.sqrt((np.square(lhs.real) + np.square(lhs.imag)).sum(1))
    else:
        scl = np.sqrt(np.square(lhs).sum(1))
    scl[scl == 0] = 1
    lhsT = lhs.T
    if len(lhsT.shape)>2.5:
        lhsT = lhsT.reshape(lhsT.shape[0],lhsT.shape[2])
        scl = np.sqrt(np.sum(np.square(lhs),axis=2)).T
        lhsT_over_scl = lhsT/scl
    else:
        lhsT_over_scl = lhsT/scl

    # Solve the least squares problem.
    # c, resids, rank, s = np.linalg.lstsq(lhsT_over_scl, rhs.T, rcond)
    # c = (c.T/scl).T
    beta = SX.sym('beta',lhsT_over_scl.shape[1],1)
    objective_fit = sum1(fabs(lhsT_over_scl@beta - rhs.T)**1) #+ lambda_tikhonov * sum1(beta**2)
    nlp_fit = {'x':beta, 'f':objective_fit}
    solvr = nlpsol('solvr', 'ipopt', nlp_fit, {'ipopt.print_level':0,'ipopt.tol':1e-12,'ipopt.acceptable_tol':1e-12})
    sol = solvr(x0=DM.ones(beta.size1()))['x']
    # my_resids = lhsT_over_scl@sol - rhs.T
    # alt_resids = lhsT_over_scl@c - rhs.T
    # my_sol = (sol/scl)
    # c = my_sol
    return sol/scl # I don't really understand this
    # why divide by scl twice?


    # # Expand c to include non-fitted coefficients which are set to zero
    # if deg.ndim > 0:
    #     if c.ndim == 2:
    #         cc = np.zeros((lmax+1, c.shape[1]), dtype=c.dtype)
    #     else:
    #         cc = np.zeros(lmax+1, dtype=c.dtype)
    #     cc[deg] = c
    #     c = cc

    # # warn on rank reduction
    # if rank != order and not full:
    #     msg = "The fit may be poorly conditioned"
    #     warnings.warn(msg, RankWarning, stacklevel=2)

    # if full:
    #     return c, [resids, rank, s, rcond]
    # else:
    #     return c


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
   


   
consumption = SX.sym('consumption', T, 1)
labour = SX.sym('labour', T, 1)
capital = SX.sym('capital', T+1, 1)


lower_bound_C = vertcat(DM.zeros(T) + C_min)    # lower bound on the consumption -> not binding anyway
lower_bound_L = vertcat(DM.zeros(T) + L_min)

upper_bound_C = vertcat(DM.zeros(T) + np.inf)
upper_bound_L = vertcat(DM.zeros(T) + L_max) # upper bound on labor also doesn't bind






# Define the start point
x_0 = vertcat(DM.ones(T), DM.ones(T),DM.ones(T+1))
l_function = 1
c_function = 1
    
for iteration in range(n_iter):

    c = np.array([]).reshape(0,1)
    l = np.array([]).reshape(0,1)
    k = np.array([]).reshape(0,1)
    z = np.array([]).reshape(0,1)

    for k0 in k0s:
        # np.random.seed(1)
        # zeta = DM(exp(AR_gen([1,-lambda_zeta],[1],T,burnin=0,scale = sigma_rho)))

        for z0 in z0s:

            zeta = []
            zeta.append(z0)
            for i in range(T-1):
                zeta.append(zeta[i]**lambda_zeta)
            zeta = DM(zeta)
    
    
    
            def condition1(consumption, labour, capital):
                return ((1/g)*(zeta[:T]*capital[:T]**alpha*labour[:T]**(1-alpha) + (1-delta) * capital[:T] - consumption[:T]) - capital[1:T+1])
                
            def condition2(consumption, labour, capital, lfunc, cfunc):
                return (((beta/g)*((1/P)*(zeta[:T-1]**lambda_zeta) * alpha*capital[1:T]**(alpha-1)* labour[1:T]**(1-alpha) + 1 - delta)*(1/consumption[1:T]))
                -  1/consumption[:T-1])

                # (logquad(capital[:T],zeta[:T],8,lfunc,cfunc)
                     # - 1/consumption[:T]) #length = T
            
            def condition3(consumption, labour, capital):
                return ((1/P)*(1-alpha)*zeta[:T]*capital[:T]**alpha*labour[:T]**(-alpha)
                    - consumption[:T]*labour[:T]**eta) #length = T
            
            def FOCs(consumption, labour, capital, lfunc, cfunc):
                return vertcat(condition1(consumption, labour, capital),condition2(consumption, labour, capital, lfunc, cfunc),
                    condition3(consumption, labour, capital))
    
    
    
            lower_bound_K = vertcat(k0, DM.zeros(T) - np.inf)
            upper_bound_K = vertcat(k0, DM.zeros(T) + np.inf)
    
            lb_x = vertcat(lower_bound_C, lower_bound_L, lower_bound_K)
            ub_x = vertcat(upper_bound_C, upper_bound_L, upper_bound_K)
            
            objective = (capital[T] - ssk) ** 2
            nonlin_con = FOCs(consumption, labour, capital, l_function, c_function)
            nlp = {'x':vertcat(consumption, labour, capital), 'f':objective, 'g':nonlin_con}
            solver = nlpsol('solver', 'ipopt', nlp,{'ipopt.print_level':5})
            solution = solver(x0=x_0,lbx=lb_x,ubx=ub_x,lbg=-1e-10,ubg=1e-10)
            sol = solution['x']
            c_sim, l_sim, k_sim = vertsplit(sol,[0,T,2*T,3*T+1])
    
            c = np.concatenate([c,c_sim[:periods]])
            l = np.concatenate([l,l_sim[:periods]])
            k = np.concatenate([k,k_sim[:periods]])
            z = np.concatenate([z,zeta[:periods]])

    c_coefs = herme2d_fit([k,z],c,degree).T
    if c_coefs.shape[0] > 1:
        c_coefs = c_coefs.T
    
    def c_function(k,z):
        c_poly = _vander_nd_flat_SYM((hermevander_casadiSYM,hermevander_casadiSYM),[k,z],[degree,degree]) @ c_coefs.T
        output_length = c_poly.shape[0]
        return fmax(c_poly,SX.zeros(output_length)+C_min)

    l_coefs = herme2d_fit([k,z],l,degree).T
    if l_coefs.shape[0] > 1:
        l_coefs = l_coefs.T
    vals = _vander_nd_flat((hermevander,hermevander),[k[:periods],zeta[:periods]],[degree,degree])
    
    def l_function(k,z):
        l_poly = _vander_nd_flat_SYM((hermevander_casadiSYM,hermevander_casadiSYM),[k,z],[degree,degree]) @ l_coefs.T
        output_length = l_poly.shape[0]
        return fmax(l_poly,SX.zeros(output_length)+L_min)

    print(k,z)
    print(c)
    print(l)
    np.savetxt(
        "NLCEQdata.csv",
        np.hstack([k,z,c,l]),
        delimiter=",",
        comments="",
        header="k,z,c,l",
    )
    print(l_coefs)
    print(c_coefs)
    print(l_function(DM.ones(5)*ssk,DM.ones(5)))
    print(c_function(DM.ones(2)*ssk,DM.ones(2)))