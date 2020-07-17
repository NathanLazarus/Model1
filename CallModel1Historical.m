clear all

DELTA   = 0.08;  %depreciation rate
ALFA    = 0.32;  %capital share
BETTA   = 0.98; %discount rate
G       = 1.014;
SIGM    = 0.9;
LAMBDAP = 0.95;
LAMBDAZ = 0.92;
sigma_Z = 0.017;
sigma_P = 0.03;
FRISCHELAS = 0.5;
STEADYSTATEL = 0.3;
MU = 0.086;
T=1000;
k0_mult=1;
startopposite = 1;
regimechanges = 0;
regime_change_frequency = 50;
order = 4;
randomseq = 2;


MultiplicativeU = 1;

import casadi.*

if MultiplicativeU
    utility_function_str = "Multiplicative";
else
    utility_function_str = "Additive";
end

shocks = ["historical_postwar_trend","historical","historical_endogenous_P_postwar_trend","historical_endogenous_P"];

    
N = length(shocks);
output_vars=10;
second_run_output = NaN([T output_vars N]);
shock_character_vector = char(shocks(1));
if string(shock_character_vector(1:10)) == "historical"
    second_run_output = NaN([33 output_vars N]);
end

first_run_output_vars=9;
first_run_output = NaN([T first_run_output_vars N]);
shock_character_vector = char(shocks(1));
if string(shock_character_vector(1:10)) == "historical"
    first_run_output = NaN([37 first_run_output_vars N]);
end

tic
parfor i = 1:N
    shock = shocks(i);
    [first_run_output(:,:,i),second_run_output(:,:,i)] = model1_run(DELTA,ALFA,BETTA,G,SIGM,LAMBDAP,LAMBDAZ,...
        sigma_Z,sigma_P,MU,FRISCHELAS,STEADYSTATEL,T,shock,k0_mult,MultiplicativeU,...
        startopposite,regimechanges,regime_change_frequency,randomseq,order);
end
toc

writematrix(second_run_output(:,:,1:2),"C:/Users/Nathan/Downloads/PerturbationMethods/Historical_Model_Outputs_K_Dynamics.xlsx",'Sheet',1,'Range','B4')
writematrix(second_run_output(:,:,3:4),"C:/Users/Nathan/Downloads/PerturbationMethods/Historical_Model_Outputs_K_Dynamics.xlsx",'Sheet',2,'Range','B4')

split_to_cell = num2cell(first_run_output, [1 2]);
stacked = vertcat(split_to_cell{:});
% writematrix(stacked,"C:/Users/Nathan/Downloads/PerturbationMethods/Model1/EncounteredPs/HistoricalSims.csv")