% model1_run.M
% Calls: model1.m num_eval.m  model1_ss_numeric.m gx_hx.m gxx_hxx.m gss_hss.m
function [rgwkcl_mat] = model1_secondrun(DELTA,ALFA,BETTA,G,SIGM,LAMBDAP,LAMBDAZ,sigma_Z,sigma_P,MU,FRISCHELAS,STEADYSTATEL,T,shock,k0_mult,MultiplicativeU,Loop_P,startopposite,randomseq)

defaults = [0.08,0.32,0.98,1.014,0.9,0.95,0.9,0.0072,0.005,0.33,0.5,0.3,250,77,1,1,1,1,0];
var={'DELTA','ALFA','BETTA','G','SIGM',...
    'LAMBDAP','LAMBDAZ','sigma_Z','sigma_P','MU',...
    'FRISCHELAS','STEADYSTATEL','T','shock','k0_mult',...
    'MultiplicativeU','Loop_P','startopposite','randomseq'};

regimechanges = 1;

for i = 1:length(defaults)
    if ~exist(var{i},'var')
        eval(sprintf('%s = %g;',var{i},defaults(i)))
    end
end
shock

if MultiplicativeU
    u = multiplicative_u;
else
    u = additive_u;
end

P       = G^(1/(1-LAMBDAP));
Loop_P_Path = 1.03:0.01:1.35;
eta     = [0 1]'; %Matrix defining driving force
T_Psims = T;

%lower the bounds on the casadi lbg's

% impulse response functions setup
irf=0;
Z_shock=sigma_Z;
T_irf = T;

% simulations setup
simulations=0;
% T=100;
% T_Psims = 40;
stats=0;
Euler_Error = 0;
Sim_P=1;

ZSTAR = 1; %steady-state value of technology shock 

addpath('C:/Users/Nathan/Downloads/casadi-windows-matlabR2016a-v3.5.1')
import casadi.*

sym_labor_supply = laborsupply(u);
intertemporal_euler_ss = dupdcp_over_dudc(u,1);
intertemporal_euler_sym = dupdcp_over_dudc(u,0);

[KSTAR,CSTAR,LSTAR,WSTAR,RSTAR,GAMA,ETA]=model1_ss_numericsetGAMAandETA(DELTA,ALFA,BETTA,G,P,FRISCHELAS,STEADYSTATEL,SIGM,ZSTAR,sym_labor_supply,intertemporal_euler_ss,u)
PSTAR = G^(1/(1-LAMBDAP));
if startopposite
    if LAMBDAP == 0.95
    Popposite = G^(1/(1-0.77));
    end
    if LAMBDAP == 0.77
    Popposite = G^(1/(1-0.95));
    end
    KSTARopposite = model1_ss_numeric(KSTAR,CSTAR,LSTAR,DELTA,ALFA,BETTA,G,Popposite,ETA,GAMA,SIGM,ZSTAR,sym_labor_supply,intertemporal_euler_ss);
    k0_mult = KSTARopposite/KSTAR;
end
    
k=KSTAR; c=CSTAR; l=LSTAR; Z=ZSTAR;
kp=k; cp=c; lp=l; Zp=Z; riskless_r_ = RSTAR;


rng(13466910+randomseq,'twister');
rho_zeta = normrnd(0,sigma_Z,[1 max(T,T_Psims)]); %should be max(T,T_Psims,LP once I add the loop start and stop as parameters)
rho_zeta(1:5)=shock;
positiveshock = 0.015;
if shock == 99
    realTFP = readmatrix('TFPshocks.csv');
    %realTFP(realTFP(:,1)>(max(realTFP(:,1))-33),:)
    rho_zeta(1:33) = realTFP(realTFP(:,1)>1984.5&realTFP(:,1)<2017.5,2);
    real_rho_zeta_back_to_1980(1:33+4) = realTFP(realTFP(:,1)>1980.5&realTFP(:,1)<2017.5,2);
end
if shock == 88
    realTFP = readmatrix('TFPshocks.csv');
    rho_zeta(1:33+4) = realTFP(realTFP(:,1)>1980.5&realTFP(:,1)<2017.5,2);
    Loop_P = 0;
    T_Psims = 33+4;
end
if shock == 77
    realTFP = readmatrix('TFPshocks.csv');
    rho_zeta(1:33) = realTFP(realTFP(:,1)>1984.5&realTFP(:,1)<2017.5,2);
    real_rho_zeta_back_to_1980(1:33+4) = realTFP(realTFP(:,1)>1980.5&realTFP(:,1)<2017.5,2);
    realProfits = readmatrix('ProfitShare.csv');
    Loop_P_Path = (1./(1-realProfits(realProfits(:,1)>1984.5&realProfits(:,1)<2017.5,2)))';
end
if shock == 66
    rho_zeta(1:5) = positiveshock;
    realProfits = readmatrix('ProfitShare.csv');
    Loop_P_Path = (1./(1-realProfits(realProfits(:,1)>1984.5&realProfits(:,1)<2017.5,2)))';
end
if shock == 55
    rho_zeta(1:5) = positiveshock;
    Loop_P = 0;
    T_Psims = 33;
end

rng(20,'twister');
rho_P = normrnd(0,sigma_P,[1 T_Psims]);
        

    if Sim_P
        if Loop_P
            P_P=Loop_P_Path;
            LP=length(P_P);
            k_P = ones([1 LP]) * KSTAR * k0_mult;
            Z_P = ones([1 LP]) * ZSTAR;
            if shock == 99 || shock == 77
                Z_P_back_to_1980 = ones([1 4]) * ZSTAR;
                for i = 1:4
                    Z_P_back_to_1980(i)=Z_P_back_to_1980(max(i-1,1))^LAMBDAZ*exp(real_rho_zeta_back_to_1980(i));
                end
                Z_P = ones([1 LP]) * Z_P_back_to_1980(4);
            end
        else
            LP = T_Psims;
            k_P = ones([1 LP]) * KSTAR * k0_mult;
            Z_P = ones([1 LP]) * ZSTAR;
            % this is from when the loop went from 2:LP instead of 1:LP, so P(1) needed to be P(1) and not PSTAR, I think P_P=ones([1 LP]) * P_func(G,PSTAR,LAMBDAP,Z_P(1),Z_P(1),Z_P(1),Z_P(1),Z_P(1),MU,rho_P(1));
            P_P=ones([1 LP])*PSTAR;
            if startopposite
                P_P = ones([1 LP])*Popposite;
            end
        end
        
        c_P(1:LP) = ones([1 LP]);
        l_P(1:LP) = ones([1 LP]);
        
        ssk(1:LP) = ones([1 LP]);
        ssl(1:LP) = ones([1 LP]);
        ssc(1:LP) = ones([1 LP]);
        
        MathematicaOutputs = readmatrix("C:/Users/Nathan/Downloads/PerturbationMethods/Model1/coefsbyP.csv");
              
        for i = 1:LP
            Z_P(i)=Z_P(max(i-1,1))^LAMBDAZ*exp(rho_zeta(i));
            if ~Loop_P
                if regimechanges&&(i>50&&i<=100)||(i>150&&i<=200)
                   P_P(i) = P_func(G,P_P(max(i-1,1)),0.77,Z_P(i),Z_P(max(i-1,1)),Z_P(max(i-2,1)),Z_P(max(i-3,1)),Z_P(max(i-4,1)),MU,rho_P(i));
                else
                   P_P(i) = P_func(G,P_P(max(i-1,1)),LAMBDAP,Z_P(i),Z_P(max(i-1,1)),Z_P(max(i-2,1)),Z_P(max(i-3,1)),Z_P(max(i-4,1)),MU,rho_P(i));
                end
            end
            if i == 4 && shock == 88
                P_P(i) = Popposite;
            end
            [k_P(i+1),c_P(i),l_P(i)]=model1_P_secondrun(P_P(i),k_P(i),Z_P(i),LAMBDAZ,sigma_Z,BETTA,DELTA,ALFA,ETA,GAMA,SIGM,G,eta,ZSTAR,KSTAR,CSTAR,LSTAR,...
                sym_labor_supply,intertemporal_euler_ss,MathematicaOutputs);
            if i == 4 && shock == 88
                k_P(i+1) = KSTAR * k0_mult;
            end
        end
        k_P = k_P(1:LP);
        w_P = w_func(k_P,l_P,P_P,Z_P,ALFA);
        r_P = big_R(k_P,l_P,P_P,Z_P,ALFA,DELTA)-1;
        y_P = y_func(k_P,l_P,Z_P,ALFA);
        g_P = [NaN,(G*y_P(2:LP)-y_P(1:LP-1))./y_P(1:LP-1)];
    
    end
    
if shock == 88
     g_P(5) = NaN;
end
rgwkcl_mat = [r_P',g_P',w_P',k_P',c_P',l_P',y_P',P_P',Z_P'];
if shock == 88
     rgwkcl_mat = rgwkcl_mat(5:37,:);
end
% filename = "C:/Users/Nathan/Downloads/PerturbationMethods/Comparing6Models.xlsx";
% sheetloc = 'J6'; %J6, AB6, AT6
% if exist('filename','var')   
%     writematrix(rgwkcl_mat,filename,'Sheet',1,'Range',sheetloc)
% end
% 
filename = 'C:/Users/Nathan/Downloads/PerturbationMethods/DifferentMUs.xlsx';
sheetloc = 'S7';
writematrix(rgwkcl_mat,filename,'Sheet','Historical Zetas','Range',sheetloc)

end