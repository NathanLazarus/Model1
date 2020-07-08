% model1_run.M
% Calls: model1.m num_eval.m  model1_ss_numeric.m gx_hx.m gxx_hxx.m gss_hss.m
function [rgwkcl_mat,firstrunresults] = model1_firstrun(DELTA,ALFA,BETTA,G,SIGM,LAMBDAP,LAMBDAZ,sigma_Z,sigma_P,MU,FRISCHELAS,STEADYSTATEL,T,shock,k0_mult,MultiplicativeU,Loop_P,startopposite,randomseq)

defaults = [0.08,0.32,0.98,1.014,0.9,0.95,0.9,0.0072,0.005,0.33,0.5,0.3,250,77,1,0,1,1,0];
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

[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f] = model1(u);

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

%Order of approximation desired 
approx = 1;

%Obtain numerical derivatives of f
num_eval

%First-order approximation
[gx,hx] = gx_hx(nfy,nfx,nfyp,nfxp);

flatten = @(A) A(:);

if approx == 2
    %Second-order approximation
    [gxx,hxx] = gxx_hxx(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx)

    [gss,hss] = gss_hss(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx,gxx,eta)
    dec_k=[KSTAR,hx(1,:),1/2*flatten(hxx(1,:,:))',1/2*hss(1)];
    dec_l=[LSTAR,gx(1,:),1/2*flatten(gxx(1,:,:))',1/2*gss(1)];
    dec_c=[CSTAR,gx(2,:),1/2*flatten(gxx(2,:,:))',1/2*gss(2)];
else
    dec_k=[KSTAR,hx(1,:),0,0,0,0,0]; 
    dec_l=[LSTAR,gx(1,:),0,0,0,0,0];
    dec_c=[CSTAR,gx(2,:),0,0,0,0,0];
end

if irf
    k(1:T_irf)=KSTAR;
    c(1:T_irf)=CSTAR;
    Z(1:T_irf)=ZSTAR;
    l(1:T_irf)=LSTAR;
    w(1:T_irf)=WSTAR;
    r(1:T_irf)=RSTAR;
    
    for i=2:T_irf
        k(i)=decision_func(dec_k,[k(i-1) Z(i-1)],[KSTAR ZSTAR],sigma_Z);
        if i==2
            Z(i) = Z(i-1)^LAMBDAZ*exp(Z_shock);
        else
            Z(i) = Z(i-1)^LAMBDAZ;
        end
        c(i)=decision_func(dec_c,[k(i) Z(i)],[KSTAR ZSTAR],sigma_Z);
        l(i)=decision_func(dec_l,[k(i) Z(i)],[KSTAR ZSTAR],sigma_Z);
    end

    w = w_func(k,l,P,Z,ALFA);
    r = little_r(k,l,P,Z,ALFA,DELTA);

    % figure(1)
    % subplot(2,3,1)
    % plot(k)
    % title('Capital ($k$)','Interpreter','latex')
    % subplot(2,3,2)
    % plot(c)
    % title('Consumption ($c$)','Interpreter','latex')
    % subplot(2,3,3)
    % plot(l)
    % title('Labor ($l$)','Interpreter','latex')
    % subplot(2,3,4)
    % plot(Z)
    % title('Productivity shock ($\zeta$)','Interpreter','latex')
    % subplot(2,3,5)
    % plot(r)
    % title('Interest rate ($r$)','Interpreter','latex')
    % subplot(2,3,6)
    % plot(w)
    % title('Wage rate ($w$)','Interpreter','latex')
    % print(['IRF_comp_1_approx_',num2str(approx)],'-djpeg','-r150')
    % close(1)
end

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
        
if simulations
    % fprintf('\n mean(rho_zeta)=%g, estimated=%g\n',0, mean(rho_zeta))
    % fprintf('sigma_Z=%g, estimated=%g\n\n', sigma_Z, std(rho_zeta))
    
    % Start from the non-stochastic steady state
    k_sim(1:T)=KSTAR;
    c_sim(1:T)=CSTAR;
    Z_sim(1:T)=ZSTAR;
    l_sim(1:T)=LSTAR;
    w_sim(1:T)=WSTAR;
    r_sim(1:T)=RSTAR;

    for i=2:T
        k_sim(i)=decision_func(dec_k,[k_sim(i-1) Z_sim(i-1)],[KSTAR ZSTAR],sigma_Z);
        Z_sim(i)=Z_sim(i-1)^LAMBDAZ*exp(rho_zeta(i));

        c_sim(i)=decision_func(dec_c,[k_sim(i) Z_sim(i)],[KSTAR ZSTAR],sigma_Z);
        l_sim(i)=decision_func(dec_c,[k_sim(i) Z_sim(i)],[KSTAR ZSTAR],sigma_Z);       
    end

    w_sim=w_func(k_sim,l_sim,P,Z_sim,ALFA);
    r_sim=little_r(k_sim,l_sim,P,Z_sim,ALFA,DELTA);
    % figure(2)
    % subplot(2,3,1)
    % plot(k_sim)
    % title('Capital ($k$)','Interpreter','latex')
    % subplot(2,3,2)
    % plot(c_sim)
    % title('Consumption ($c$)','Interpreter','latex')
    % subplot(2,3,3)
    % plot(l_sim)
    % title('Labor ($l$)','Interpreter','latex')
    % subplot(2,3,4)
    % plot(Z_sim)
    % title('Productivity shock ($\zeta$)','Interpreter','latex')    
    % subplot(2,3,5)
    % plot(r_sim)
    % title('Interest rate ($r$)','Interpreter','latex')
    % subplot(2,3,6)
    % plot(w_sim)
    % title('Wage rate ($w$)','Interpreter','latex')
    % print(['SIM_comp_1_approx_',num2str(approx)],'-djpeg','-r150')
    % close(2)
    if stats
        M=mean([k_sim;c_sim;l_sim;Z_sim;r_sim;w_sim],2);
        V=var([k_sim;c_sim;l_sim;Z_sim;r_sim;w_sim],0,2);
        Max=max([k_sim;c_sim;l_sim;Z_sim;r_sim;w_sim],[],2);
        Min=min([k_sim;c_sim;l_sim;Z_sim;r_sim;w_sim],[],2);
        fprintf('\n mean(k)=%g, mean(c)=%g, mean(l)=%g, mean(Z)=%g, mean(r)=%g, mean(w)=%g\n',M)
        fprintf('\n var(k)=%g, var(c)=%g, var(l)=%g, var(Z)=%g, var(r)=%g, var(w)=%g\n',V)
        fprintf('\n max(k)=%g, max(c)=%g, max(l)=%g, max(Z)=%g, max(r)=%g, max(w)=%g\n',Max)
        fprintf('\n min(k)=%g, min(c)=%g, min(l)=%g, min(Z)=%g, min(r)=%g, min(w)=%g\n',Min)
    end
    
    if Euler_Error
        error = ones([T-1 5]);
        Marg_Util = c_sim.^-1 + l_sim.^ETA;
        
        euler_eqs = subs(f,{sym('cp'),sym('lp')},...
                {dec_c*[1,(sym('kp')-KSTAR),(sym('Zp')-ZSTAR),(sym('kp')-KSTAR)^2,(sym('kp')-KSTAR)*(sym('Zp')-ZSTAR),(sym('kp')-KSTAR)*(sym('Zp')-ZSTAR),(sym('Zp2')-2*sym('Zp')*ZSTAR+ZSTAR^2),sigma_Z^2]',...
                dec_l*[1,(sym('kp')-KSTAR),(sym('Zp')-ZSTAR),(sym('kp')-KSTAR)^2,(sym('kp')-KSTAR)*(sym('Zp')-ZSTAR),(sym('kp')-KSTAR)*(sym('Zp')-ZSTAR),(sym('Zp2')-2*sym('Zp')*ZSTAR+ZSTAR^2),sigma_Z^2]'...
                });
            
        for t = 1:T-1
            euler_eqs_t = subs(euler_eqs,{sym('kp'),sym('lp'),sym('k'),sym('Z'),sym('l'),sym('c')},...
                {k_sim(t+1),l_sim(t+1),k_sim(t),Z_sim(t),l_sim(t),c_sim(t)});
            error(t,:) = eval(subs(euler_eqs_t,{sym('Zp'),sym('Zp2')},{Z_sim(t)^LAMBDAZ*exp(sigma_Z^2/2),(Z_sim(t)^LAMBDAZ)^2*exp(2*sigma_Z^2)}))';
        end
        
        unit_free_error = error./Marg_Util(1:T-1)';
        
        mean_relative_error = mean(log10(abs(unit_free_error)+1e-15));
        fprintf('\nOrder of magnitude of the mean error for equation 1: %g, equation 2: %g, equation 3: %g\n\n', mean_relative_error(1),mean_relative_error(2),mean_relative_error(3))
        max_relative_error = max(log10(abs(unit_free_error)+1e-15));
        fprintf('\nOrder of magnitude of the max error for equation 1: %g, equation 2: %g, equation 3: %g\n\n', max_relative_error(1),max_relative_error(2),max_relative_error(3))
    end
end


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
        E_l(1:LP) = ones([1 LP]);
        E_c(1:LP) = ones([1 LP]);
        
        ssk(1:LP) = ones([1 LP]);
        ssl(1:LP) = ones([1 LP]);
        ssc(1:LP) = ones([1 LP]);
        stable_dkplusdk(1:LP) = ones([1 LP]);
        
        riskless_r_P(1:LP) = ones([1 LP]);
        riskless_r_alt(1:LP) = ones([1 LP]);
              
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
            [k_P(i+1),c_P(i),l_P(i),ssk(i),ssl(i),ssc(i),stable_dkplusdk(i)]=model1_P_firstrun(P_P(i),approx,k_P(i),Z_P(i),LAMBDAZ,sigma_Z,BETTA,DELTA,ALFA,ETA,GAMA,SIGM,G,eta,ZSTAR,KSTAR,CSTAR,LSTAR,...
                sym_labor_supply,intertemporal_euler_ss,intertemporal_euler_sym,fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f,rho_zeta(min(i+1,T)));
            if i == 4 && shock == 88
                k_P(i+1) = KSTAR * k0_mult;
            end
        end
        k_P = k_P(1:LP);
        w_P = w_func(k_P,l_P,P_P,Z_P,ALFA);
        E_Z = Z_P;
        for i = 2:LP
            E_Z(i) = Z_P(i-1)^LAMBDAZ*exp(sigma_Z^2/2);
        end
        r_P = big_R(k_P,l_P,P_P,Z_P,ALFA,DELTA)-1;
        y_P = y_func(k_P,l_P,Z_P,ALFA);
        g_P = [NaN,(G*y_P(2:LP)-y_P(1:LP-1))./y_P(1:LP-1)];
        another_riskless = 1./(BETTA*(c_P./(G*E_c)).^SIGM)-1;
        % figure(3)
        % subplot(2,3,1)
        % plot(k_P)
        % title('Capital ($k$)','Interpreter','latex')
        % subplot(2,3,2)
        % plot(c_P)
        % title('Consumption ($c$)','Interpreter','latex')
        % subplot(2,3,3)
        % plot(l_P)
        % title('Labor ($l$)','Interpreter','latex')
        % subplot(2,3,4)
        % plot(Z_P)
        % title('Productivity shock ($\zeta$)','Interpreter','latex')
        % subplot(2,3,5)
        % plot(r_P)
        % title('Interest rate ($r$)','Interpreter','latex')
        % subplot(2,3,6)
        % plot(P_P)
        % title('Markup ($P$)','Interpreter','latex')
        % print(['SIM_P_comp_5_approx_',num2str(approx)],'-djpeg','-r200')
        % close(3)
    
    end
    

rgwkcl_mat = [r_P',g_P',w_P',k_P',c_P',l_P',y_P',P_P',Z_P']; % ,riskless_r_alt',another_riskless',E_c'
if shock == 88
     rgwkcl_mat = rgwkcl_mat(5:37,:);
     rgwkcl_mat(1,3) = NaN;
end
if exist('filename','var')   
    writematrix(rgwkcl_mat,filename,'Sheet',1,'Range',sheetloc)
end

dim = size(P_P);
firstrunresults = [P_P',ssk',ssl',ssc',stable_dkplusdk',ones([dim(2) 1])*GAMA];
% writematrix(firstrunresults,join(["C:/Users/Nathan/Downloads/PerturbationMethods/Model1/EncounteredPs/Mu",string(round(MU*100)),"positiveshock.csv"],""));
writematrix(firstrunresults,join(["C:/Users/Nathan/Downloads/PerturbationMethods/Model1/EncounteredPs/RegimeChangesPositiveShock",randomseq,".csv"],""))
end