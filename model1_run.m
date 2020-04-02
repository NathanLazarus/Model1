% model1_run.M
% Calls: model1.m num_eval.m  model1_ss_numeric.m gx_hx.m gxx_hxx.m gss_hss.m
function [rgwkcl_mat] = model1_run(DELTA,ALFA,BETTA,G,SIGM,LAMBDAP,LAMBDAZ,sigma_Z,sigma_P,MU,FRISCHELAS,STEADYSTATEL,T,shock,k0_mult,MultiplicativeU)

defaults = [0.12,0.32,0.96,1.018,0.9,0.95,0.896,0.0072,0.005,0.35,0.5,0.3,10,0,1,0];
var={'DELTA','ALFA','BETTA','G','SIGM','LAMBDAP','LAMBDAZ','sigma_Z','sigma_P','MU','FRISCHELAS','STEADYSTATEL','T','shock','k0_mult','MultiplicativeU'};

for i = 1:length(defaults)
    if ~exist(var{i},'var')
        eval(sprintf('%s = %g',var{i},defaults(i)))
    end
end

if MultiplicativeU
    u = multiplicative_u;
else
    u = additive_u;
end

[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f] = model1(u);

P       = G^(1/(1-LAMBDAP));
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
Loop_P=0;
Sim_P=1;
Euler_Error = 0;

ZSTAR = 1; %steady-state value of technology shock 

[KSTAR,CSTAR,LSTAR,WSTAR,RSTAR,GAMA,ETA]=model1_ss_numericsetGAMAandETA(DELTA,ALFA,BETTA,G,P,FRISCHELAS,STEADYSTATEL,SIGM,ZSTAR,u);


k=KSTAR; c=CSTAR; l=LSTAR; Z=ZSTAR;
kp=k; cp=c; lp=l; Zp=Z;

%Order of approximation desired 
approx = 2;

%Obtain numerical derivatives of f
num_eval

%First-order approximation
[gx,hx] = gx_hx(nfy,nfx,nfyp,nfxp);

flatten = @(A) A(:);

if approx == 2
    %Second-order approximation
    [gxx,hxx] = gxx_hxx(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx); 

    [gss,hss] = gss_hss(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx,gxx,eta);
    dec_k=[KSTAR,hx(1,:),1/2*flatten(hxx(1,:,:))',1/2*hss(1)];
    dec_l=[LSTAR,gx(1,:),1/2*flatten(gxx(1,:,:))',1/2*gss(1)];
    dec_c=[LSTAR,gx(2,:),1/2*flatten(gxx(2,:,:))',1/2*gss(2)];
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
        l(i)=decision_func(dec_c,[k(i) Z(i)],[KSTAR ZSTAR],sigma_Z);
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

rng(13466910,'twister');
rho_zeta = normrnd(0,sigma_Z,[1 max(T,T_Psims)]);
rho_zeta(1:5)=shock;

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
        k_sim(i)=decision_func(dec_k,[k(i-1) Z(i-1)],[KSTAR ZSTAR],sigma_Z);
        Z_sim(i)=Z_sim(i-1)^LAMBDAZ*exp(rho_zeta(i));

        c_sim(i)=decision_func(dec_c,[k(i) Z(i)],[KSTAR ZSTAR],sigma_Z);
        l_sim(i)=decision_func(dec_c,[k(i) Z(i)],[KSTAR ZSTAR],sigma_Z);       
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
        error = ones([T-1 4]);
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
    
    
    if Loop_P
        Pl=1.03:0.01:1.35;
        LP=length(Pl);
        
        k_P(1:LP)=KSTAR*k0_mult;
        Z_P=ones([1 LP])*exp(rho_zeta(1));
        
        [~,cstart,lstart]=model1_P(Pl(1),approx,k_P(1),Z_P(1),Z_P(1),LAMBDAZ,sigma_Z,BETTA,DELTA,ALFA,ETA,GAMA,SIGM,G,eta,ZSTAR,u,1);
        c_P(1:LP)=cstart;
        l_P(1:LP)=lstart;
        
      
        for i=2:LP
            Z_P(i)=Z_P(i-1)^LAMBDAZ*exp(rho_zeta(i));
            [k_P(i),c_P(i),l_P(i)]=model1_P(Pl(i),approx,k_P(i-1),Z_P(i-1),Z_P(i),LAMBDAZ,sigma_Z,BETTA,DELTA,ALFA,ETA,GAMA,SIGM,G,eta,ZSTAR,u,0);
        end

        w_P=w_func(k_P,l_P,P,Z_P,ALFA);
        r_P=little_r(k_P,l_P,P,Z_P,ALFA,DELTA);
        y_P=y_func(k_P,l_P,Z_P,ALFA);
        g_P = [NaN,(G*y_P(2:LP)-y_P(1:LP-1))./y_P(1:LP-1)];

        % figure(3)
        % subplot(2,3,1)
        % plot(Pl,k_P)
        % title('Capital ($k$)','Interpreter','latex')
        % subplot(2,3,2)
        % plot(Pl,c_P)
        % title('Consumption ($c$)','Interpreter','latex')
        % subplot(2,3,3)
        % plot(Pl,l_P)
        % title('Labor ($l$)','Interpreter','latex')
        % subplot(2,3,4)
        % plot(Pl,Z_P)
        % title('Productivity shock ($\zeta$)','Interpreter','latex')
        % subplot(2,3,5)
        % plot(Pl,r_P)
        % title('Interest rate ($r$)','Interpreter','latex')
        % subplot(2,3,6)
        % plot(Pl,w_P)
        % title('Wage rate ($w$)','Interpreter','latex')
        % print(['SIM_P_comp_1_approx_',num2str(approx)],'-djpeg','-r150')
        %close(3)
    
    end
        
    if Sim_P
        LP = T_Psims;
        
        k_P = ones([1 LP]) * KSTAR * k0_mult;
        Z_P = ones([1 LP]) * exp(rho_zeta(1));
                
        PSTAR = G^(1/(1-LAMBDAP));
        P_P=ones([1 LP]) * G * PSTAR^LAMBDAP * Z_P(1)^(MU) * Z_P(1)^(MU^2) * Z_P(1)^(MU^3) * Z_P(1)^(MU^4) * Z_P(1)^(MU^5)*exp(rho_P(1));
        
        [~,cstart,lstart,rstart,wstart,ystart] = model1_P(P_P(1),approx,k_P(1),Z_P(1),Z_P(1),LAMBDAZ,sigma_Z,BETTA,DELTA,ALFA,ETA,GAMA,SIGM,G,eta,ZSTAR,u,1);
        c_P(1:LP) = cstart;
        l_P(1:LP) = lstart;
        
              
        for i = 2:LP
            Z_P(i)=Z_P(i-1)^LAMBDAZ*exp(rho_zeta(i));
            P_P(i) = G * P_P(i-1)^LAMBDAP * Z_P(i-1)^(MU) * Z_P(max(i-2,1))^(MU^2) * Z_P(max(i-3,1))^(MU^3) * Z_P(max(i-4,1))^(MU^4) * Z_P(max(i-5,1))^(MU^5)*exp(rho_P(i));
            [k_P(i),c_P(i),l_P(i)]=model1_P(P_P(i),approx,k_P(i-1),Z_P(i-1),Z_P(i),LAMBDAZ,sigma_Z,BETTA,DELTA,ALFA,ETA,GAMA,SIGM,G,eta,ZSTAR,u,0);
        end
        w_P = w_func(k_P,l_P,P,Z_P,ALFA);
        r_P = little_r(k_P,l_P,P,Z_P,ALFA,DELTA);
        y_P = y_func(k_P,l_P,Z_P,ALFA);
        g_P = [NaN,(G*y_P(2:LP)-y_P(1:LP-1))./y_P(1:LP-1)];
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
    
rgwkcl_mat = [r_P',g_P',w_P',k_P',c_P',P_P',Z_P',l_P'];
if exist('filename','var')   
    writematrix(rgwkcl_mat,filename,'Sheet',1,'Range',sheetloc)
end

end
