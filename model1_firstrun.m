% model1_run.M
% Calls: model1.m num_eval.m  model1_ss_numeric.m gx_hx.m gxx_hxx.m gss_hss.m
function [firstrunresults] = model1_firstrun(DELTA,ALFA,BETTA,G,SIGM,LAMBDAP,LAMBDAZ,sigma_Z,sigma_P,MU,FRISCHELAS,STEADYSTATEL,T,shock,k0_mult,MultiplicativeU,Loop_P,startopposite,randomseq)

defaults = {0.08,0.32,0.98,1.014,...
    0.9,0.95,0.92,0.017,...
    0.03,0.086,0.5,0.3,...
    200,"none",1,1,...
    1,1,50,2,4};

var = ["DELTA","ALFA","BETTA","G",...
    "SIGM","LAMBDAP","LAMBDAZ","sigma_Z",...
    "sigma_P","MU","FRISCHELAS","STEADYSTATEL",...
    "T","shock","k0_mult","MultiplicativeU",...
    "startopposite","regimechanges","regime_change_frequency","randomseq","order"];

for i = 1:length(defaults)
    if ~exist(var{i},"var")
        if class(defaults{i}) == "string"
            eval(sprintf('%s = "%s";',var(i),defaults{i}))
        else
            eval(sprintf("%s = %g;",var(i),defaults{i}))
        end
    end
end

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

[~,~,~,~,~,GAMA,ETA]=model1_ss_numericsetGAMAandETA(DELTA,ALFA,BETTA,G,value_of_P_where_LSTAR_equals_STEADYSTATEL,FRISCHELAS,STEADYSTATEL,SIGM,ZSTAR,sym_labor_supply,intertemporal_euler_ss,u);

if startopposite
    PSTARopposite = G^(1/(1-LAMBDAPopposite));
%     [KSTARopposite,CSTARopposite,LSTARopposite,WSTARopposite,RSTARopposite]=model1_ss_numeric(1,0.3,0.3,DELTA,ALFA,BETTA,G,PSTARopposite,ETA,GAMA,SIGM,ZSTAR,sym_labor_supply,intertemporal_euler_ss);
%     stockSTARopposite = G*((PSTARopposite-1)/PSTARopposite)*y_func(KSTARopposite,LSTARopposite,ZSTAR,ALFA)/((1+RSTARopposite) - G);
%     k0_mult = KSTARopposite/KSTAR;
end

rng(13466910+randomseq,'twister');
rho_zeta = normrnd(0,sigma_Z,[1 T]);
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
        historical_Z_path = zeros([1 T]) + ZSTAR^LAMBDAZ*exp(rho_Z(1));
        rho_P = zeros([1 T]) + log(True_P_Path(1)/P_func(G,PSTAR,LAMBDAP,historical_Z_path(1),historical_Z_path(1),historical_Z_path(1),historical_Z_path(1),historical_Z_path(1),MU,0));
        for i=2:T
            historical_Z_path(i)=historical_Z_path(i-1)^LAMBDAZ*exp(rho_Z(i));
            historical_P_implied_by_Z = P_func(G,True_P_Path(max(i-1,1)),LAMBDAP,historical_Z_path(i),historical_Z_path(max(i-1,1)),historical_Z_path(max(i-2,1)),historical_Z_path(max(i-3,1)),historical_Z_path(max(i-4,1)),MU,0);
            rho_P(i) = log(True_P_Path(i)/historical_P_implied_by_Z);
        end
    end
    if string(shock) == "historical_endogenous_P"
        rho_P = zeros([1 T]);
    end
end
        
        
Z_sim = zeros([1 T]) + ZSTAR;
P_sim = zeros([1 T]) + PSTAR;
if startopposite
    P_sim = zeros([1 T]) + PSTARopposite;
end
ssk(1:T) = NaN([1 T]);
ssl(1:T) = NaN([1 T]);
ssc(1:T) = NaN([1 T]);
stable_dkplusdk(1:T) = NaN([1 T]);
        
for i = 1:T
    Z_sim(i)=Z_sim(max(i-1,1))^LAMBDAZ*exp(rho_zeta(i));
    if regimechanges && mod(floor((i-1)/regime_change_frequency),2) == 1
        P_sim(i) = P_func_greater_than_1(G,P_sim(max(i-1,1)),LAMBDAPopposite,Z_sim(i),Z_sim(max(i-1,1)),Z_sim(max(i-2,1)),Z_sim(max(i-3,1)),Z_sim(max(i-4,1)),MU,rho_P(i));
    else
        P_sim(i) = P_func_greater_than_1(G,P_sim(max(i-1,1)),LAMBDAP,Z_sim(i),Z_sim(max(i-1,1)),Z_sim(max(i-2,1)),Z_sim(max(i-3,1)),Z_sim(max(i-4,1)),MU,rho_P(i));
    end
    if i == 4 && string(shock) == "historical_endogenous_P" && startopposite
        P_sim(i) = PSTARopposite;
    end
    if i == 4 && string(shock) == "historical_endogenous_P" && ~startopposite
        P_sim(i) = PSTAR;
    end
    if abs(P_sim(i) - 1) > 1e-9
        [ssk(i),ssl(i),ssc(i),stable_dkplusdk(i)]=model1_P_firstrun(P_sim(i),LAMBDAZ,sigma_Z,BETTA,DELTA,ALFA,ETA,GAMA,SIGM,G,eta,ZSTAR,1,0.3,0.3,...
            sym_labor_supply,intertemporal_euler_ss,fx,fxp,fy,fyp,f);
    else
        [ssk(i),ssl(i),ssc(i),stable_dkplusdk(i)]=model1_P_firstrun(P_sim(i),LAMBDAZ,sigma_Z,BETTA,DELTA,ALFA,ETA,GAMA,SIGM,G,eta,ZSTAR,1,0.3,0.3,...
            sym_labor_supply,intertemporal_euler_ss,fx_no_stock,fxp_no_stock,fy_no_stock,fyp_no_stock,f_no_stock);
    end
        
end
    

dim = size(P_sim);
firstrunresults = [P_sim',ssk',ssl',ssc',stable_dkplusdk',ones([dim(2) 1])*GAMA];