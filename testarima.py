from statsmodels.tsa.arima_process import arma_generate_sample as AR_gen
import numpy as np
np.random.seed(15)
lambda_zeta = 0.92
T = 1000000
x = AR_gen([1,-lambda_zeta],[1],T,burnin=0,scale = 0.0072)
np.savetxt(
    "ar1data.csv",
    np.exp(x),
    delimiter=",",
    comments="",
    header="vals",
)
y = AR_gen([1,lambda_zeta],[1],T,burnin=0,scale = 0.0072)
np.savetxt(
    "ar1dataAlt2.csv",
    np.exp(y),
    delimiter=",",
    comments="",
    header="vals",
)
rho = np.random.randn(T+1)*0.0072
toZeta = [1]
for i in range(T+1):
    toZeta.append(toZeta[i]**lambda_zeta*np.exp(rho[i]))
zeta = np.array(toZeta[1:T+2])
np.savetxt(
    "ar1dataAlt.csv",
    zeta,
    delimiter=",",
    comments="",
    header="vals",
)