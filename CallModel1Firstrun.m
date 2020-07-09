clear all

DELTA   = 0.08;  %depreciation rate
ALFA    = 0.32;  %capital share
BETTA   = 0.98; %discount rate
G       = 1.014;
SIGM    = 0.9;
LAMBDAZ = 0.92;
sigma_Z = 0.017;
sigma_P = 0.03;
FRISCHELAS = 0.5;
STEADYSTATEL = 0.3;
MU = 0.086;
T=100;
k0_mult=1;
startopposite = 1;
regimechanges = 0;
regime_change_frequency = 50;
LAMBDAP = 0.95;
order = 4;
randomseq = 2;

MultiplicativeU = 1;


import casadi.*

shocks = ["historical","historical_endogenous_P","historical_endogenous_P_postwar_trend"];
    
N = length(shocks);
output_vars=9;
output = NaN([T output_vars N]);
shock_character_vector = char(shocks(1));
if string(shock_character_vector(1:10)) == "historical"
    output = NaN([37 output_vars N]);
end
tic
parfor i = 1:N
    shock = shocks(i);
    [output(:,:,i)] = model1_firstrun(DELTA,ALFA,BETTA,G,SIGM,LAMBDAP,LAMBDAZ,sigma_Z,sigma_P,MU,FRISCHELAS,STEADYSTATEL,T,shock,k0_mult,MultiplicativeU,startopposite,regimechanges,regime_change_frequency,randomseq,order);
end
toc

splitA = num2cell(output, [1 2]);
stacked = vertcat(splitA{:});
writematrix(stacked,"C:/Users/Nathan/Downloads/PerturbationMethods/Model1/EncounteredPs/HistoricalSims.csv")