% model1_run.M
% Calls: model1.m num_eval.m  model1_ss_numeric.m gx_hx.m gxx_hxx.m gss_hss.m
function [exogenous_var_results,rgwkcl_mat] = model1_run(DELTA,ALFA,BETTA,G,SIGM,LAMBDAP,LAMBDAZ,...
    sigma_Z,sigma_P,MU,FRISCHELAS,STEADYSTATEL,T,shock,k0_mult,MultiplicativeU,...
    startopposite,regimechanges,regime_change_frequency,randomseq,order)

defaults = {0.08,0.32,0.98,1.014,...
    0.9,0.95,0.92,0.017,...
    0.03,0.086,0.5,0.3,...
    10,"historical",1,0,...
    1,1,50,...
    2,4};

var = ["DELTA","ALFA","BETTA","G",...
    "SIGM","LAMBDAP","LAMBDAZ","sigma_Z",...
    "sigma_P","MU","FRISCHELAS","STEADYSTATEL",...
    "T","shock","k0_mult","MultiplicativeU",...
    "startopposite","regimechanges","regime_change_frequency",...
    "randomseq","order"];

for i = 1:length(defaults)
    if ~exist(var{i},"var")
        if class(defaults{i}) == "string"
            eval(sprintf('%s = "%s";',var(i),defaults{i}))
        else
            eval(sprintf("%s = %g;",var(i),defaults{i}))
        end
    end
end

folder = fileparts(which(mfilename));
addpath(join([folder,'\','MyHelperFunctions'],""))
addpath(join([folder,'\','SchmittGroheUribeHelperFunctions'],""))

if MultiplicativeU
    u = multiplicative_u;
    utility_function_str = "Multiplicative";
else
    u = additive_u;
    utility_function_str = "Additive";
end


ZSTAR   = 1; %steady-state value of technology shock 
PSTAR   = G^(1/(1-LAMBDAP));
eta     = [0 1]'; %Matrix defining driving force
[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f] = model1(u);
[fx_no_stock,fxp_no_stock,fy_no_stock,fyp_no_stock,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,f_no_stock] = model1_no_stock(u);


LAMBDAPhigh = 0.95;
LAMBDAPlow = 0.8;
if LAMBDAP == LAMBDAPhigh
    LAMBDAPopposite = LAMBDAPlow;
end
if LAMBDAP == LAMBDAPlow
    LAMBDAPopposite = LAMBDAPhigh;
end
value_of_P_where_LSTAR_equals_STEADYSTATEL = G^(1/(1-LAMBDAPlow));

import casadi.*

sym_labor_supply = laborsupply(u);
intertemporal_euler_ss = dupdcp_over_dudc(u,1);

[KSTAR,~,~,~,~,GAMA,ETA]=model1_ss_numericsetGAMAandETA(DELTA,ALFA,BETTA,G,value_of_P_where_LSTAR_equals_STEADYSTATEL,FRISCHELAS,STEADYSTATEL,SIGM,ZSTAR,sym_labor_supply,intertemporal_euler_ss,u);
[~,~,~,~,~,~,multiplicativeETA]=model1_ss_numericsetGAMAandETA(DELTA,ALFA,BETTA,G,value_of_P_where_LSTAR_equals_STEADYSTATEL,FRISCHELAS,STEADYSTATEL,SIGM,ZSTAR,laborsupply(multiplicative_u),dupdcp_over_dudc(multiplicative_u,1),multiplicative_u);
[~,~,~,~,~,~,additiveETA]=model1_ss_numericsetGAMAandETA(DELTA,ALFA,BETTA,G,value_of_P_where_LSTAR_equals_STEADYSTATEL,FRISCHELAS,STEADYSTATEL,SIGM,ZSTAR,laborsupply(additive_u),dupdcp_over_dudc(additive_u,1),additive_u);

if startopposite
    PSTARopposite = G^(1/(1-LAMBDAPopposite));
    [KSTARopposite,CSTARopposite,LSTARopposite,WSTARopposite,RSTARopposite]=model1_ss_numeric(1,0.3,0.3,DELTA,ALFA,BETTA,G,PSTARopposite,ETA,GAMA,SIGM,ZSTAR,sym_labor_supply,intertemporal_euler_ss);
%     stockSTARopposite = G*((PSTARopposite-1)/PSTARopposite)*y_func(KSTARopposite,LSTARopposite,ZSTAR,ALFA)/((1+RSTARopposite) - G);
    k0_mult = KSTARopposite/KSTAR;
end

if string(shock) == "historical_postwar_trend"
    tfpfile = "C:/Users/Nathan/Downloads/PerturbationMethods/Parameterizations/TFPshocks_relative_to_postwar_trend.csv";
    shock = "historical";
elseif string(shock) == "historical_endogenous_P_postwar_trend"
    tfpfile = "C:/Users/Nathan/Downloads/PerturbationMethods/Parameterizations/TFPshocks_relative_to_postwar_trend.csv";
    shock = "historical_endogenous_P";
else
    tfpfile = "C:/Users/Nathan/Downloads/PerturbationMethods/Parameterizations/TFPshocks.csv";
end

rng(13466910+randomseq,'twister');
rho_Z = normrnd(0,sigma_Z,[1 T]);
if ~(string(shock) == "none"||string(shock) == "historical"|| string(shock) == "historical_endogenous_P")
    rho_Z(1:5)=shock;
end
rng(123140+randomseq,"twister");
rho_P = normrnd(0,sigma_P,[1 T]);

if (string(shock) == "historical" || string(shock) == "historical_endogenous_P")
    T = 37;
    realTFPshocks = readmatrix(tfpfile);
    rho_Z = zeros([1 T]) + realTFPshocks(realTFPshocks(:,1)>1980.5&realTFPshocks(:,1)<2017.5,2);
    if string(shock) == "historical"
        realProfits = readmatrix("C:/Users/Nathan/Downloads/PerturbationMethods/Parameterizations/ProfitShare.csv");
        True_P_Path = (1./(1-realProfits(realProfits(:,1)>1980.5&realProfits(:,1)<2017.5,2)))';
        historical_Z_path = zeros([1 T]) + ZSTAR;
        rho_P = NaN([1 T]);
        for i=1:T
            historical_Z_path(i)=historical_Z_path(max(i-1,1))^LAMBDAZ*exp(rho_Z(i));
            if i == 1
                if startopposite
                    historical_P_implied_by_Z = P_func(G,PSTARopposite,LAMBDAP,historical_Z_path(i),historical_Z_path(max(i-1,1)),historical_Z_path(max(i-2,1)),historical_Z_path(max(i-3,1)),historical_Z_path(max(i-4,1)),MU,0);
                else
                    historical_P_implied_by_Z = P_func(G,PSTAR,LAMBDAP,historical_Z_path(i),historical_Z_path(max(i-1,1)),historical_Z_path(max(i-2,1)),historical_Z_path(max(i-3,1)),historical_Z_path(max(i-4,1)),MU,0);
                end
            else
                historical_P_implied_by_Z = P_func(G,True_P_Path(max(i-1,1)),LAMBDAP,historical_Z_path(i),historical_Z_path(max(i-1,1)),historical_Z_path(max(i-2,1)),historical_Z_path(max(i-3,1)),historical_Z_path(max(i-4,1)),MU,0);
            end
            rho_P(i) = log(True_P_Path(i)/historical_P_implied_by_Z);
        end
    end
    if string(shock) == "historical_endogenous_P"
        rho_P = zeros([1 T]);
    end
end
  
try
    MathematicaOutputs = readmatrix("C:/Users/Nathan/Downloads/PerturbationMethods/Model1/coefsbyP.csv");
catch
    MathematicaOutputs = [];
end

Z_sim = zeros([1 T]) + ZSTAR;
P_sim = zeros([1 T]) + PSTAR;
k_sim = zeros([1 T+1]) + KSTAR;
if startopposite
    P_sim = zeros([1 T]) + PSTARopposite;
    k_sim = zeros([1 T+1]) + KSTARopposite;
end
l_sim(1:T) = NaN([1 T]);
c_sim(1:T) = NaN([1 T]);
stock_sim(1:T) = NaN([1 T]);
ssk(1:T) = NaN([1 T]);
ssl(1:T) = NaN([1 T]);
ssc(1:T) = NaN([1 T]);
ssstock(1:T) = NaN([1 T]);
stable_dcdk(1:T) = NaN([1 T]);
        
for i = 1:T
    Z_sim(i)=Z_sim(max(i-1,1))^LAMBDAZ*exp(rho_Z(i));
    if regimechanges && mod(floor((i-1)/regime_change_frequency),2) == 1
        P_sim(i) = P_func_greater_than_1(G,P_sim(max(i-1,1)),LAMBDAPopposite,Z_sim(i),Z_sim(max(i-1,1)),Z_sim(max(i-2,1)),Z_sim(max(i-3,1)),Z_sim(max(i-4,1)),MU,rho_P(i));
    else
        P_sim(i) = P_func_greater_than_1(G,P_sim(max(i-1,1)),LAMBDAP,Z_sim(i),Z_sim(max(i-1,1)),Z_sim(max(i-2,1)),Z_sim(max(i-3,1)),Z_sim(max(i-4,1)),MU,rho_P(i));
%         P_func_greater_than_1(G,P_sim(max(i-1,1)),LAMBDAP,Z_sim(i),Z_sim(max(i-1,1)),Z_sim(max(i-2,1)),Z_sim(max(i-3,1)),Z_sim(max(i-4,1)),MU,rho_P(i))
    end
    if i == 4 && string(shock) == "historical_endogenous_P" && startopposite
        P_sim(i) = PSTARopposite;
    end
    if i == 4 && string(shock) == "historical_endogenous_P" && ~startopposite
        P_sim(i) = PSTAR;
    end
    if i == 5 && (string(shock) == "historical" || string(shock) == "historical_endogenous_P")
        k_sim(i) = KSTAR*k0_mult;
    end
    if abs(P_sim(i) - 1) > 1e-9
        [ssk(i),ssl(i),ssc(i),ssstock(i),stable_dcdk(i),k_sim(i+1),c_sim(i),l_sim(i),stock_sim(i)]=model1_P(P_sim(i),k_sim(i),Z_sim(i),LAMBDAZ,sigma_Z,BETTA,DELTA,ALFA,ETA,GAMA,SIGM,G,eta,ZSTAR,1,0.3,0.3,...
            order,sym_labor_supply,intertemporal_euler_ss,fx,fxp,fy,fyp,f,MathematicaOutputs);
    else
        [ssk(i),ssl(i),ssc(i),ssstock(i),stable_dcdk(i),k_sim(i+1),c_sim(i),l_sim(i),stock_sim(i)]=model1_P(P_sim(i),k_sim(i),Z_sim(i),LAMBDAZ,sigma_Z,BETTA,DELTA,ALFA,ETA,GAMA,SIGM,G,eta,ZSTAR,1,0.3,0.3,...
            order,sym_labor_supply,intertemporal_euler_ss,fx_no_stock,fxp_no_stock,fy_no_stock,fyp_no_stock,f_no_stock,MathematicaOutputs);
    end
    k_sim = k_sim(1:T);
end
    
dim = size(P_sim);
exogenous_var_results = [P_sim',ssk',ssl',ssc',ssstock',stable_dcdk',ones([dim(2) 1])*GAMA,ones([dim(2) 1]) * multiplicativeETA,ones([dim(2) 1]) * additiveETA];


w_sim = w_func(k_sim,l_sim,P_sim,Z_sim,ALFA);
r_sim = little_r(k_sim,l_sim,P_sim,Z_sim,ALFA,DELTA);
y_sim = y_func(k_sim,l_sim,Z_sim,ALFA);
g_sim = [NaN,(G*y_sim(2:T)-y_sim(1:T-1))./y_sim(1:T-1)];
profits_sim = ((P_sim-1)./P_sim).*y_sim;

rgwkcl_mat = [r_sim',g_sim',w_sim',k_sim',c_sim',l_sim',stock_sim',y_sim',P_sim',Z_sim'];

if (string(shock) == "historical" || string(shock) == "historical_endogenous_P")
    rgwkcl_mat = rgwkcl_mat(5:37,:);
end

