import numpy as np
import matplotlib.pyplot as plt
from casadi import *
from scipy.optimize import fsolve
from scipy import integrate
from scipy.special import h_roots
import csv
from multiprocessing import Pool
import os

# Ppath = [(1.025,25),(1.1,25),(1.2,25),(1.4,25),(2,25),(3,25),(4.5,50)]
# Mordecaispath = range(1025,1415,15)
periods = 180
itercount = 1
# k0 = DM.ones(1)*3
k0s = range(3,4,1)
Ps = range(1250,1500,250)
nrow_sol_path = 5
sol_path = np.array([]).reshape(nrow_sol_path,0)
# for P,periods in Ppath:
# for Pval in Mordecaispath:
    # P = Pval/1000

ncol_reg = 5
reg_dat = np.array([]).reshape(0,ncol_reg)

for Pval in Ps:
    P = Pval/1000
    for k0val in k0s:
        k0 = DM.ones(1)*k0val/2


        T = 200
        
        alpha = 0.32
        g = 1.016
        delta = 0.12
        
        # theta = 1000
        # P = theta/(theta-1)
        
        beta = 0.96
        
        # taxes = [0, 0, 0]
        
        eta = 1.5
        
        Lmax = 1
        
        np.random.seed(itercount)
        sigma_rho = 0.01
        # sigma_rho=0.00000001
        rho = np.random.randn(T+1)*sigma_rho
        zeta = [1]
        lambda_zeta = 0.96
        for i in range(T+1):
            zeta.append(zeta[i]**lambda_zeta*np.exp(rho[i]))
        zeta = DM(np.array(zeta[1:T+2]))
        # zeta = DM.ones(T+1)

        
        
        # solve for steady state
        cstar = SX.sym('cstar', 1, 1)
        lstar = SX.sym('lstar', 1, 1)
        kstar = SX.sym('kstar', 1, 1)
        
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
        star_solution = star_solver(x0=star_x_0,lbg=-0.000000001,ubg=0.000000001)
        # print(star_solution['x'])
        star_sol = star_solution['x']
        print(star_sol)
        ssc, ssl, ssk = star_sol[0], star_sol[1], star_sol[2]
        
            # print(ssc, ssl, ssk)
            
            
            # print(ssk)
            # print('=')
            # print((1/g)*(ssk**alpha*ssl**(1-alpha)+ssk*(1-delta)-ssc))
            # print('rstar')
            # print((1/P)*(alpha)*ssk**(alpha-1)*ssl**(1-alpha)-delta)
            # print('=')
            # print(g/beta - 1)
            # print('ssl')
            # print(ssl)
            # print((((1-alpha)/P)*ssk**alpha*(1/ssc))**(1/(eta+alpha)))
            # print(ssk**alpha*ssl**(1-alpha) + (1-delta) * ssk - ssc - g*ssk)
            
            
            # ssc = 1
            # ssl = 1
            # ssk = 1
            
            # def R(capital,labour):
            #   return (1/P)*alpha*capital[:T]**(alpha-1)*labour**(1-alpha)
            # # R = 1 + r
            # def wage(capital,labour):
            #   return (1/P)*(1-alpha)*capital[:T]**alpha*labour**(-alpha)
            
            # def utility(consumption, labour, sigma, eta):
            #     # Standard additive CRRA utility function of consumption and labour
            #     # Input: Vector 1xT consumption, labour; 1x1 sigma, eta
            #     # Output: utility 1xT
                
            #     u_consumption = log(consumption) #**(1-sigma)/(1-sigma)
              
            #     u_labour = (1-labour)**(1+eta)/(1+eta)
                
            #     u = u_consumption + (u_labour)
            #     return u
            
            
            # def wage_fun(T, T_ret):
            #     # Return a 1xT vector of wages
            #     # input: T, and retirement age T_ret
            #     res = np.maximum(1.5,(1/2 + np.array(range(1,T+1)) * (1 - np.array(range(1,T+1))/T))) / 16
            #     res[(T_ret-1):T] = 0
            #     return res
            
            # def tax(r, capital, consumption, wage, labour, taxes, T):
            #     # Calculate the combined tax on labour income, consumption and capital
            #     # income
            
            #     tc = taxes[0]
            #     tk = taxes[1]
            #     tl = taxes[2]
              
            #     return r*capital[1:T+1] * tk + wage * labour *tl + consumption*tc
            
            
            # def asset_constraint(capital, labour, consumption, T):
            #     temp = ((1/g)*((R(capital,labour) - delta) * capital[:T] + wage(capital, labour) * labour[:T] - consumption[:T])
            #       - capital[1:T+1])
            #     return vertcat([0],temp)
        


        def l(k,z):
            return 0.7933207-0.0271536*k+0.1006549*z
        def c(k,z):
            return -0.1570917+0.1948116*k+0.6911522*z

        # def e(rho,zeta,capital):
            #     # ((1/(-0.1379+0.12963*capital+0.01188*P+0.7238591*(zeta**lambda_zeta*rho)))* #1/c
            #     #     (1.3971670-0.0204440*capital-0.5433395*P+0.1792362*(zeta**lambda_zeta*rho))* #L
            #     #     (zeta**lambda_zeta*rho)* #zeta
            #     #     (1/(rho*sigma_rho*np.sqrt(2*np.pi)))*np.exp(-(np.log(rho)**2)/(2*sigma_rho**2))) #support of rho

            #     rho_length = rho.shape[0]
            #     zeta_length = zeta.shape[0]

            #     # print((zeta.reshape(zeta_length,1)**lambda_zeta*rho.reshape(1,rho_length))*capital.reshape(zeta_length,1))
            #     # print(zeta.reshape(zeta_length,1)**lambda_zeta*rho.reshape(1,rho_length))
            #     return ((beta/g)*((1/P)*
            #         (reshape(zeta,zeta_length,1)**lambda_zeta @ rho.reshape(1,rho_length)**sigma_rho) @ 
            #         alpha*reshape(capital,zeta_length,1)**(alpha-1)*
            #         (l(reshape(capital,zeta_length,1),(reshape(zeta,zeta_length,1)**lambda_zeta @ rho.reshape(1,rho_length)**sigma_rho)))**(1-alpha) +
            #         1 - delta)*
            #         (1/(c(reshape(capital,zeta_length,1),(reshape(zeta,zeta_length,1)**lambda_zeta @ rho.reshape(1,rho_length)**sigma_rho))))*
            #         repmat((1/(rho.reshape(1,rho_length)*np.sqrt(2*np.pi)))*np.exp(-(np.log(rho.reshape(1,rho_length))**2)/(2)),zeta_length)) #support of rho)
            # def e2(rho,zeta,capital):
            #     # ((1/(-0.1379+0.12963*capital+0.01188*P+0.7238591*(zeta**lambda_zeta*rho)))* #1/c
            #         # (1.3971670-0.0204440*capital-0.5433395*P+0.1792362*(zeta**lambda_zeta*rho))* #L
            #         # (zeta**lambda_zeta*rho)* #zeta
            #         # (1/(rho*sigma_rho*np.sqrt(2*np.pi)))*np.exp(-(np.log(rho)**2)/(2*sigma_rho**2))) #support of rho
            #     # print((zeta.reshape(zeta_length,1)**lambda_zeta*rho.reshape(1,rho_length))*capital.reshape(zeta_length,1))
            #     # print(zeta.reshape(zeta_length,1)**lambda_zeta*rho.reshape(1,rho_length))
            #     return ((beta/g)*((1/P)*
            #         (zeta**lambda_zeta*rho)*
            #         alpha*capital**(alpha-1)*
            #         (0.7933207-0.0271536*capital+0.1006549*(zeta**lambda_zeta*rho))**(1-alpha) +
            #         1 - delta)*
            #         (1/(-0.1570917+0.1948116*capital+0.6911522*(zeta**lambda_zeta*rho)))*
            #         (1/(rho*sigma_rho*np.sqrt(2*np.pi)))*np.exp(-(np.log(rho)**2)/(2*sigma_rho**2))) #support of rho)
            # # zeta_t**lambda_zeta*sigma_rho/2
            # # print(e(np.array([1,1.01,1.02,1.03]),np.array([1,2]),np.array([1,2])))
            # # def reimann(zeta,capital,n,fineness):
            #     # samplepoints = np.array(range(n))/fineness+0.5/fineness
            #     # sampleweights = np.ones(n)/fineness
            #     # return np.dot(e(samplepoints,zeta,capital),sampleweights)
            # # def trapezoid(zeta,capital,n,fineness):
            #     # samplepoints = np.array(range(n))/fineness+1e-13
            #     # samplepoints2 = np.array(range(n))/fineness+1/fineness+1e-13
            #     # sampleweights = np.ones(n)/fineness
            #     # return np.dot((e(samplepoints,zeta,capital)+e(samplepoints2,zeta,capital))/2,sampleweights)
            # # def simpson13(zeta,capital,n,fineness):
            # #     samplepoints = np.array(range(n))/fineness+1e-13
            # #     sampleweights = np.append(np.tile(np.array([2/3,4/3]),n//2),1/3)
            # #     sampleweights[0]=1/3
            # #     sampleweights = sampleweights/fineness
            # #     return e(samplepoints,zeta,capital)@sampleweights
            # # def boole(zeta,capital,n,fineness):
            # #     samplepoints = np.array(range(n))/fineness+1e-13
            # #     sampleweights = np.append(np.tile(np.array([14,32,12,32]),n//4),7)
            # #     sampleweights[0]=7
            # #     sampleweights = sampleweights*2/(45*fineness)
            # #     print(e(samplepoints,zeta,capital).shape)
            # #     print(reshape(sampleweights,n,1).shape)
            # #     return e(samplepoints,zeta,capital)@sampleweights

        def transformede(rho,zeta,capital):
            rho_length = rho.shape[0]
            zeta_length = zeta.shape[0]
            return ((beta/g)*((1/P)*
                (reshape(zeta,zeta_length,1)**lambda_zeta @ rho.reshape(1,rho_length)**sigma_rho) @ 
                alpha*reshape(capital,zeta_length,1)**(alpha-1)*
                (l(reshape(capital,zeta_length,1),(reshape(zeta,zeta_length,1)**lambda_zeta @ rho.reshape(1,rho_length)**sigma_rho)))**(1-alpha) +
                1 - delta)*
                (1/(c(reshape(capital,zeta_length,1),(reshape(zeta,zeta_length,1)**lambda_zeta @ rho.reshape(1,rho_length)**sigma_rho)))))

        def logquad(zeta,capital,n):
            samplepoints,sampleweights = h_roots(n)
            #print(e(samplepoints,zeta,capital).shape)
            print(reshape(sampleweights,n,1).shape)
            return (transformede(np.exp(np.sqrt(2)*samplepoints),zeta,capital)/np.sqrt(np.pi))@sampleweights
        

        print(1/ssc)
        # print(boole(DM.ones(2),DM.ones(2)*ssk,80001,305))
        print('logquad')
        print(logquad(DM.ones(2),DM.ones(2)*ssk,20))

        def condition1(consumption,labour,capital,mu):
            temp = ((1/g)*(zeta[:T]*capital[:T]**alpha*labour[:T]**(1-alpha) + (1-delta) * capital[:T] - consumption[:T]) - capital[1:T+1] - mu)
            return temp
        
        def condition2(consumption,labour,capital,mu):
            # temp = ((beta/g)*((1/P)*zeta[1:T+1]*alpha*capital[1:T+1]**(alpha-1)*labour[1:T+1]**(1-alpha) + 1 - delta)*(1/consumption[1:T+1])
            # intvals = np.zeros(T)
            # for t in range(T):
            #     intvals[t] = integrate.quad(e, 0, np.inf,args = (zeta[t],capital[t+1]))
            # temp = (intvals
            temp = (logquad(zeta[:T],capital[:T],20)
                 - 1/consumption[:T] - mu)
            # samplepoints = np.array(range(0,100000))/100+0.005
            # sampleweights = np.ones((100000,1))/100
            # for t in range(T):
            #     print(integrate.quad(e, 0, np.inf,args = (zeta[t],capital[t+1])))
            #     print(zeta[t])
            #     print(capital[t+1])
            #     print(1/ssc)
            return temp #length = T
        
        def condition3(consumption, labour, capital,mu):
            temp = ((1/P)*(1-alpha)*zeta[:T]*capital[:T]**alpha*labour[:T]**(-alpha)
                - consumption[:T]*labour[:T]**eta - mu)
            return temp #length = T
        
        # def condition4(consumption,labour,capital,mu):
        #     temp = ((1/g)*(capital[:T]**alpha*labour[:T]**(1-alpha) + (1-delta) * capital[:T] - consumption[:T])
        #         - capital[1:T+1] + mu)
        #     return temp #length = T
        
        # def condition5(consumption,labour,capital,mu):
        #     temp = ((beta/g)*((1/P)*alpha*capital[:T]**(alpha-1)*labour[:T]**(1-alpha) + 1 - delta)*consumption[:T]
        #         - vertcat(consumption[1:T],ssc) + mu)
        #     return temp #length = T
        
        # def condition6(consumption, labour, capital,mu):
        #     temp = ((1/P)*(1-alpha)*capital[:T]**alpha*labour[:T]**(-alpha)
        #         - consumption[:T]*labour[:T]**eta + mu)
        #     return temp #length = T
        
        # def FOCs(consumption, labour, capital, mu):
            # return vertcat(condition1(consumption, labour, capital, mu),condition2(consumption, labour, capital, mu),
                # condition3(consumption, labour, capital, mu),condition4(consumption, labour, capital, mu),
                # condition5(consumption, labour, capital, mu),condition6(consumption, labour, capital, mu))
        
        def FOCs(consumption, labour, capital, mu):
            return vertcat(condition1(consumption, labour, capital, mu),condition2(consumption, labour, capital, mu),
                condition3(consumption, labour, capital, mu))
        
        
        consumption = SX.sym('consumption', T+1, 1)
        labour = SX.sym('labour', T+1, 1)
        capital = SX.sym('capital', T+1, 1)
        # mu = SX.sym('mu', T, 1)
        mu = DM.zeros(T)
        
        objective = -sum1(consumption)
        
        lower_bound_C = vertcat(DM.ones(T)*0.001, ssc)    # lower bound on the consumption -> not binding anyway
        lower_bound_L = vertcat(DM.zeros(T), ssl)
        lower_bound_K = vertcat(k0, -1000 * DM.ones(T-1), ssk)
        
        upper_bound_C = vertcat(DM.ones(T)*np.inf,ssc)
        upper_bound_L = vertcat(DM.ones(T)*Lmax - 0.001, ssl)
        upper_bound_K = vertcat(k0, DM.ones(T-1)*np.inf, ssk)
        
        
        lb_x = vertcat(lower_bound_C, lower_bound_L, lower_bound_K) #,DM.zeros(T))
        ub_x = vertcat(upper_bound_C, upper_bound_L, upper_bound_K) #,DM.ones(T)*np.inf)
        
        
        
        # Define the start point
        x_0 = vertcat(DM.ones(T+1), DM.zeros(T+1)+0.5,DM.ones(T+1)) #,DM.ones(T))
        
        
        nonlin_con = FOCs(consumption, labour, capital, mu)
        
        
        nlp = {'x':vertcat(consumption,labour,capital), 'f':objective, 'g':nonlin_con}
        solver = nlpsol('solver', 'ipopt', nlp,{'ipopt.print_level':0})
        # solution = solver(x0=x_0,lbx=lb_x,ubx=ub_x,lbg=vertcat(DM.ones(3*T)*-np.inf,DM.zeros(3*T)-0.0001),ubg=vertcat(DM.zeros(3*T)+0.0001,DM.ones(3*T)*np.inf))
        solution = solver(x0=x_0,lbx=lb_x,ubx=ub_x,lbg=vertcat(DM.zeros(3*T)-0.0001),ubg=vertcat(DM.zeros(3*T)+0.0001))
        sol = solution['x']
        # print(sol[0:T+1],sol[T+1:2*T+2],sol[2*T+1:])
        c, l, k = sol[0:T+1],sol[T+1:2*T+2],sol[2*T+2:]
        # print(integrate.quad(e, 0, np.inf,args = (zeta,k)))
    #    cstar, lstar, kstar = sol[2*T//3],sol[(5*T)//3+1],sol[(8*T)//3+2]
    #     print('**************')
    #     print(cstar,lstar,kstar)
    #     print(ssc, ssl, ssk)
        r = (1/P)*(alpha)*sol[2*T+2:3*T+3]**(alpha-1)*sol[T+1:2*T+2]**(1-alpha)-delta
    #     # rstar = R(sol[2*T:],sol[T:2*T])
    #     # print(rstar)
    #     # print(g/beta - delta)
    #     # print(wage(sol[2*T:],sol[T:2*T]))
    #     # print((1/P)*(1-alpha)*sol[2*T:3*T]**alpha*sol[T:2*T]**(-alpha))
    #     # print(kstar)
    #     # print((1/g)*(kstar**alpha*lstar**(1-alpha)+kstar*(1-delta)+cstar))
    #     print('kstar')
    #     print(kstar)
    #     print('=')
    #     print((1/g)*(kstar**alpha*lstar**(1-alpha)+kstar*(1-delta)-cstar))
    #     print('rstar')
    #     print((1/P)*(alpha)*kstar**(alpha-1)*lstar**(1-alpha)-delta)
    #     print('=')
    #     print(g/beta - 1)
    #     print('lstar')
    #     print(lstar)
    #     print((((1-alpha)/P)*kstar**alpha*(1/cstar))**(1/(eta+alpha)))
    #     print((1/g)*(kstar**alpha*lstar**(1-alpha) + (1-delta) * kstar - cstar) - kstar)
    #     print(np.array(vertcat(c[:periods],l[:periods],k[:periods],np.ones(periods)*P,r[:periods])).reshape(nrow_sol_path,periods))
        sol_path = np.concatenate([sol_path,np.array(vertcat(c[:periods],l[:periods],k[:periods],np.ones(periods)*P,r[:periods])).reshape(nrow_sol_path,periods)],axis=1)
        print(sol_path)
    #     k0 = k[periods]
        itercount = itercount +1
        reg_dat = np.concatenate([reg_dat,np.array(horzcat(c[:periods],l[:periods],k[:periods],zeta[:periods],np.ones(periods)*P)).reshape(periods,ncol_reg)],axis=0)


# print(reg_dat)
np.savetxt("Model1ConsumptionRegressionData.csv", reg_dat, delimiter=",")
# print(sol_path)

# def plot_solution(solution,T):
#     plt.figure()
#     plt.plot(solution[0:T+1],'.')
#     plt.title('Consumption')
#     plt.figure()
#     plt.plot(solution[T+1:2*T+2],'.')
#     plt.title('Labour')
#     plt.figure()
#     plt.plot(solution[2*T+2:3*T+3],'.')
#     plt.title('Capital')
#     plt.figure()
#     plt.plot(r,'.')
#     plt.title('Interest Rate')
#     plt.show()

# # def plot_solution2(solution,T):
# #     plt.figure()
# #     plt.plot(solution[0:T],'.')
# #     plt.title('Consumption')
# #     plt.figure()
# #     plt.plot(solution[T:2*T],'.')
# #     plt.title('Labour')
# #     plt.figure()
# #     plt.plot(solution[2*T:3*T],'.')
# #     plt.title('Capital')
# #     plt.figure()
# #     plt.plot(solution[4*T:5*T],'.')
# #     plt.title('Interest Rate')
# #     plt.figure()
# #     plt.plot(solution[3*T:4*T],'.')
# #     plt.title('P')
# #     plt.show()

def plot_solution2(solution,T):
    plt.clf()
    fig = plt.figure()
    plt.subplot(2, 2, 1)
    plt.plot(solution[0:T],'.')
    plt.title('Consumption')
    plt.subplot(2, 2, 2)
    plt.plot(solution[2*T:3*T],'.')
    plt.title('Capital')
    plt.subplot(2, 2, 3)
    plt.plot(solution[4*T:5*T],'.')
    plt.title('Interest Rate')
    plt.subplot(2, 2, 4)
    plt.plot(solution[3*T:4*T],'.')
    plt.title('P')
    plt.tight_layout()
    fileloc = 'C:/Users/Nathan/Downloads/flder.png'
    if os.path.isfile(fileloc):
        os.remove(fileloc)
    fig.savefig(fileloc)
    # plt.show()

total_periods_simmed = sol_path.shape[1]
# print(total_periods_simmed)
plot_solution2(sol_path.reshape(nrow_sol_path*total_periods_simmed,1), total_periods_simmed)
print(sol_path)