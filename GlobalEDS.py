from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
# import pandas as pd
import numpy as np
from numba import jit
from multiprocessing import Pool, freeze_support, cpu_count
from statsmodels.tsa.arima_process import arma_generate_sample as AR_gen
from casadi import *
from scipy.special import h_roots

# model_outputs = pd.read_csv(
#     "C:/Users/Nathan/Downloads/PerturbationMethods/Model1/NLCEQdata.csv"
# )

T = 40
burnin = 5
sigma_rho = 0.072
lambda_zeta = 0.96
delta = 0.08
alpha = 0.32
g = 1.014

beta = 0.98
eta = 2
P = 1
L_max = np.inf
L_min = 1e-6
C_min = 1e-6

# # @jit(cache = True)
# def l_function(k,z):
#         return 1

# # @jit(cache = True)
# def c_function(k,z):
#     return 0.5*k
def k_function(k,z):
    return k ** 0.9 + 0.1 * (z - 1)

k = np.zeros(T+1)
c = np.zeros(T)
l = np.zeros(T)

np.random.seed(1)
zeta = np.array(np.exp(AR_gen(ar = [1,-lambda_zeta], ma = [1],nsample = T,burnin=0,scale = sigma_rho)))
k[0] = 1
for i in range(T):
    # c[i] = c_function(k[i], zeta[i])
    # l[i] = l_function(k[i], zeta[i])
    # k[i+1] = (k[i] * (1 - delta) + zeta[i] * k[i]**alpha * l[i] ** (1-alpha) - c[i])/g
    k[i+1] = k_function(k[i],zeta[i])


def condition1(consumption, labor, capital, zeta, capitalplus):
    return ((1/g)*(zeta*capital**alpha*labor**(1-alpha) + (1-delta) * capital - consumption) - capitalplus)

def condition3(consumption, labor, capital, zeta):
    return ((1/P)*(1-alpha)*zeta*capital**alpha*labor**(-alpha)
        - consumption*labor**eta) #length = T

consumption = SX.sym('consumption',T,1)
labor = SX.sym('labor',T,1)
objective = 1
x_0 = DM.ones(vertcat(consumption,labor).shape[0])

lower_bound_C = vertcat(DM.zeros(T) + C_min)    # lower bound on the consumption -> not binding anyway
upper_bound_C = vertcat(DM.zeros(T) + np.inf)
lower_bound_L = vertcat(DM.zeros(T) + L_min)
upper_bound_L = vertcat(DM.zeros(T) + L_max) # upper bound on labor also doesn't bind

nonlin_con = vertcat(condition1(consumption, labor, k[:-1], zeta, k[1:]),condition3(consumption, labor, k[:-1], zeta))
nlp = {'x':vertcat(consumption, labor), 'f':objective, 'g':nonlin_con}
solver = nlpsol('solver', 'ipopt', nlp,{'ipopt.print_level':0})
solution = solver(x0=x_0,lbx=vertcat(lower_bound_C,lower_bound_L),ubx=vertcat(upper_bound_C,upper_bound_L),lbg=-1e-10,ubg=1e-10)
print(solution['x'])
c, l = vertsplit(solution['x'],[0,consumption.shape[0],consumption.shape[0]+labor.shape[0]])
print(k,c,l)
state_data = np.vstack([
    k[:T],
    zeta,
    c.T,
    l.T
    ]).T
state_data = state_data[burnin:,:]

# def model_outputs_to_df_of_states(model_outputs):
#     state_df = pd.concat(
#         [
#             model_outputs["k"],
#             model_outputs["z"],
#         ],
#         axis=1,
#     ).dropna()
#     return state_df


# @jit(cache = True)
def sqdist(x, y):
    return ((x - y) ** 2).sum(-1)

def eds(data_as_principal_components, epsilon):
    candidate_points = np.hstack([
        np.arange(data_as_principal_components.shape[0]).reshape(data_as_principal_components.shape[0],1),
        data_as_principal_components])
    print(candidate_points.shape)
    inEDS = np.empty(candidate_points.shape)
    counter = 0
    while candidate_points.size > 0:
        addingtoEDS = candidate_points[0]
        inEDS[counter] = addingtoEDS
        candidate_points = candidate_points[
            sqdist(candidate_points[:,1:], addingtoEDS[1:]) > epsilon ** 2
        ]
        counter += 1
    else:
        inEDS = inEDS[0:counter]

    return inEDS


def examine_different_epsilon_values(
    data, desiredEDSpoints=14, startepsilon=0.9, endepsilon=5, epsilonstep=0.1
):
    standardize = StandardScaler()
    normalized = standardize.fit_transform(data)
    pca = PCA(n_components = data.shape[1])
    as_principal_components = pca.fit_transform(normalized)

    for epsilon in np.arange(startepsilon, endepsilon, epsilonstep):
        EDSapplied = eds(as_principal_components, epsilon)
        if EDSapplied.shape[0] <= desiredEDSpoints:
            break
    return standardize.inverse_transform(
            pca.inverse_transform(EDSapplied)
        )

def eds_fixed_epsilon(
    data, epsilon = 0.1
):
    standardize = StandardScaler()
    normalized = standardize.fit_transform(data)
    pca = PCA(n_components = data.shape[1])
    as_principal_components = pca.fit_transform(normalized)
    indices, entries = np.hsplit(eds(as_principal_components, epsilon),np.array([1]))

    # return standardize.inverse_transform(
    #         pca.inverse_transform(
    #             entries
    #         )
    #     )

    return indices

# if __name__ == "__main__":
#     with Pool(cpu_count()) as p:
#         pointstouse = p.map(
#             eds_fixed_epsilon,
#             (
#                 state_data
#             ),
#         )
#     print(pointstouse)

#     np.savetxt("eds_points.csv", np.vstack(pointstouse))

n_quadrature_nodes = 5
points, weights = h_roots(n_quadrature_nodes)
def integrationnodes(znow, lambda_zeta, nodes):
    return np.outer(znow ** lambda_zeta, exp(nodes)).flatten(order = 'F')

print(np.outer(np.array([1,2,3]),np.array([4,5,6])))
print(np.outer(np.array([1,2,3]),np.array([4,5,6])).flatten(order = 'C'))
print(eds_fixed_epsilon(state_data[:,0:2],1))
indices = eds_fixed_epsilon(state_data[:,0:2],1)
indices = np.squeeze(indices.astype('int'))
n_eds_points = indices.shape[0]
know, znow, cnow, lnow = np.hsplit(state_data[indices],4)


zplus = integrationnodes(znow, lambda_zeta, points)
kplus = np.repeat(state_data[indices+1,0],n_quadrature_nodes)
kplusplus = k_function(kplus,zplus)


cplusSX = SX.sym('cplusSX',n_eds_points*n_quadrature_nodes,1)
lplusSX = SX.sym('lplusSX',n_eds_points*n_quadrature_nodes,1)

objective = 1
x_0 = DM.ones(vertcat(cplusSX,lplusSX).shape[0])

lower_bound_C = vertcat(DM.zeros(cplusSX.shape[0]) + C_min)    # lower bound on the consumption -> not binding anyway
upper_bound_C = vertcat(DM.zeros(cplusSX.shape[0]) + np.inf)
lower_bound_L = vertcat(DM.zeros(lplusSX.shape[0]) + L_min)
upper_bound_L = vertcat(DM.zeros(lplusSX.shape[0]) + L_max) # upper bound on labor also doesn't bind

print(condition1(cplusSX, lplusSX, kplus, zplus, kplusplus))
nonlin_con = vertcat(condition1(cplusSX, lplusSX, kplus, zplus, kplusplus),condition3(cplusSX, lplusSX, kplus, zplus))
nlp = {'x':vertcat(cplusSX, lplusSX), 'f':objective, 'g':nonlin_con}
solver = nlpsol('solver', 'ipopt', nlp,{'ipopt.print_level':0})
solution = solver(x0=x_0,lbx=vertcat(lower_bound_C,lower_bound_L),ubx=vertcat(upper_bound_C,upper_bound_L),lbg=-1e-10,ubg=1e-10)
cplus, lplus = vertsplit(solution['x'],cplusSX.shape[0])
print(np.repeat(weights,n_eds_points))
one_for_fixed_point_iteration = (((beta/g)*((1/P)*(zplus**lambda_zeta) * alpha*kplus**(alpha-1)* lplus**(1-alpha) + 1 - delta)*(1/cplus))*cnow)
