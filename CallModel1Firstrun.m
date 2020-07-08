DELTA   = 0.08;  %depreciation rate
ALFA    = 0.32;  %capital share
BETTA   = 0.98; %discount rate
G       = 1.014;
SIGM    = 0.9;
LAMBDAP = 0.95;
LAMBDAZ = 0.9;
sigma_Z = 0.0072;
sigma_P = 0.005;
MU      = 0.33;
FRISCHELAS = 0.5;
STEADYSTATEL = 0.3;
T=100;


addpath('C:/Users/Nathan/Downloads/casadi-windows-matlabR2016a-v3.5.1')
import casadi.*


Loop_P = 0;

shocks = [0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015,0.015];
mus = [0.3,0.3,0.33,0.33,0.38,0.38,0.43,0.43,0.48,0.48];
k0_mult=1;
startopposite = 1;
multU = 0;
    
N = length(shocks);
output_vars=9;
output = NaN([T output_vars N]);
tic
parfor i = 1:N
    MU = mus(i)
    shock = shocks(i)
    if mod(i,2) == 0
        [output(:,:,i),~] = model1_firstrun(DELTA,ALFA,BETTA,G,SIGM,LAMBDAP,LAMBDAZ,sigma_Z,sigma_P,MU,FRISCHELAS,STEADYSTATEL,T,shock,k0_mult,multU,Loop_P,startopposite);
    end
end
toc