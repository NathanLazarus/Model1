% model1_run.M
% Calls: model1.m num_eval.m  model1_ss_numeric.m gx_hx.m gxx_hxx.m gss_hss.m

clear all

[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f] = model1;

% Numerical Evaluation
% Steady State and Parameter Values
% choose omega/1+eta such that lstar = .2
BETTA   = 0.96; %discount rate
DELTA   = 0.12;  %depreciation rate
ALFA    = 0.32;  %capital share
ETA     = 1.5;
P       = 1.2;
G       = 1.016;
LAMBDAZ = 0.96;
LAMBDAP = 0.785; %reg 1931-1985 and 1889-1901 & 1985- to estimate lambdap
MU = 0.5; %reg 1931-1985 and 1889-1901 & 1985- to estimate lambdap
eta     = [0 1]'; %Matrix defining driving force
sigma_Z = 0.0072;

% impulse response functions setup
irf=1;
Z_shock=sigma_Z;

% simulations setup
simulations=1;
T=50;
T_Psims = 30;
stats=1;
Loop_P=1;
Sim_P=1;
Euler_Error = 1;

ZSTAR = 1; %steady-state value of technology shock 

[KSTAR,CSTAR,LSTAR,WSTAR,RSTAR]=model1_ss_numeric(DELTA,ALFA,BETTA,G,P,ETA,ZSTAR);

k=KSTAR; c=CSTAR; l=LSTAR; Z=ZSTAR;
kp=k; cp=c; lp=l; Zp=Z;

%Order of approximation desired 
approx = 2;

%Obtain numerical derivatives of f
num_eval

%First-order approximation
[gx,hx] = gx_hx(nfy,nfx,nfyp,nfxp);

if approx == 2
    %Second-order approximation
    [gxx,hxx] = gxx_hxx(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx); 

    [gss,hss] = gss_hss(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx,gxx,eta);
    dec_k=[KSTAR,hx(1,1),hx(1,2),1/2*hxx(1,1,1),hxx(1,2,1),1/2*hxx(1,2,2),1/2*hss(1)]; 
    dec_l=[LSTAR,gx(1,1),gx(1,2),1/2*gxx(1,1,1),gxx(1,2,1),1/2*gxx(1,2,2),1/2*gss(1)];
    dec_c=[CSTAR,gx(2,1),gx(2,2),1/2*gxx(2,1,1),gxx(2,2,1),1/2*gxx(2,2,2),1/2*gss(2)];
else
    dec_k=[KSTAR,hx(1,1),hx(1,2),0,0,0,0]; 
    dec_l=[LSTAR,gx(1,1),gx(1,2),0,0,0,0];
    dec_c=[CSTAR,gx(2,1),gx(2,2),0,0,0,0];
end

if irf
    k(1)=KSTAR;
    c(1)=CSTAR;
    Z(1)=ZSTAR;
    l(1)=LSTAR;
    w(1)=WSTAR;
    r(1)=RSTAR;
    
    if approx == 1
        k(2)=dec_k*[1,k(1)-KSTAR,Z(1)-1,0,0,0,0]';
        Z(2)=Z(1)^LAMBDAZ*exp(Z_shock);
        c(2)=dec_c*[1,k(2)-KSTAR,Z(2)-1,0,0,0,0]';
        l(2)=dec_l*[1,k(2)-KSTAR,Z(2)-1,0,0,0,0]';
    else
        k(2)=dec_k*[1,k(1)-KSTAR,Z(1)-1,(k(1)-KSTAR)^2,(k(1)-KSTAR)*(Z(1)-1),(Z(1)-1)^2,sigma_Z^2]';
        Z(2)=Z(1)^LAMBDAZ*exp(Z_shock);
        c(2)=dec_c*[1,k(2)-KSTAR,Z(2)-1,(k(2)-KSTAR)^2,(k(2)-KSTAR)*(Z(2)-1),(Z(2)-1)^2,sigma_Z^2]';
        l(2)=dec_l*[1,k(2)-KSTAR,Z(2)-1,(k(2)-KSTAR)^2,(k(2)-KSTAR)*(Z(2)-1),(Z(2)-1)^2,sigma_Z^2]';
    end
    w(2)=Z(2)/P*(1-ALFA)*k(2)^ALFA*l(2)^(-ALFA);
    r(2)=Z(2)/P*ALFA*k(2)^(ALFA-1)*l(2)^(1-ALFA)-DELTA;

    for i=3:40
        k(i)=dec_k*[1,k(i-1)-KSTAR,Z(i-1)-1,(k(i-1)-KSTAR)^2,(k(i-1)-KSTAR)*(Z(i-1)-1),(Z(i-1)-1)^2,sigma_Z^2]';
        Z(i)=Z(i-1)^LAMBDAZ;
        c(i)=dec_c*[1,k(i)-KSTAR,Z(i)-1,(k(i)-KSTAR)^2,(k(i)-KSTAR)*(Z(i)-1),(Z(i)-1)^2,sigma_Z^2]';
        l(i)=dec_l*[1,k(i)-KSTAR,Z(i)-1,(k(i)-KSTAR)^2,(k(i)-KSTAR)*(Z(i)-1),(Z(i)-1)^2,sigma_Z^2]';
        w(i)=Z(i)/P*(1-ALFA)*k(i)^ALFA*l(i)^(-ALFA);
        r(i)=Z(i)/P*ALFA*k(i)^(ALFA-1)*l(i)^(1-ALFA)-DELTA;
    end
    figure(1)
    subplot(2,3,1)
    plot(k)
    title('Capital ($k$)','Interpreter','latex')
    subplot(2,3,2)
    plot(c)
    title('Consumption ($c$)','Interpreter','latex')
    subplot(2,3,3)
    plot(l)
    title('Labor ($l$)','Interpreter','latex')
    subplot(2,3,4)
    plot(Z)
    title('Productivity shock ($\zeta$)','Interpreter','latex')
    subplot(2,3,5)
    plot(r)
    title('Interest rate ($r$)','Interpreter','latex')
    subplot(2,3,6)
    plot(w)
    title('Wage rate ($w$)','Interpreter','latex')
    print(['IRF_comp_1_approx_',num2str(approx)],'-djpeg','-r150')
    close(1)
end

if simulations
    rng(13466917,'twister');
    rho_zeta = normrnd(0,sigma_Z,[1 T]);
    fprintf('\n mean(rho_zeta)=%g, estimated=%g\n',0, mean(rho_zeta))
    fprintf('sigma_Z=%g, estimated=%g\n\n', sigma_Z, std(rho_zeta))
    
    % Start from the non-stochastic steady state
    k_sim(1:T)=KSTAR;
    c_sim(1:T)=CSTAR;
    Z_sim(1:T)=ZSTAR;
    l_sim(1:T)=LSTAR;
    w_sim(1:T)=WSTAR;
    r_sim(1:T)=RSTAR;
    for i=2:T
        k_sim(i)=dec_k*[1,k_sim(i-1)-KSTAR,Z_sim(i-1)-1,(k_sim(i-1)-KSTAR)^2,(k_sim(i-1)-KSTAR)*(Z_sim(i-1)-1),(Z_sim(i-1)-1)^2,sigma_Z^2]';
        Z_sim(i)=Z_sim(i-1)^LAMBDAZ*exp(rho_zeta(i));
        c_sim(i)=dec_c*[1,k_sim(i)-KSTAR,Z_sim(i)-1,(k_sim(i)-KSTAR)^2,(k_sim(i)-KSTAR)*(Z_sim(i)-1),(Z_sim(i)-1)^2,sigma_Z^2]';
        l_sim(i)=dec_l*[1,k_sim(i)-KSTAR,Z_sim(i)-1,(k_sim(i)-KSTAR)^2,(k_sim(i)-KSTAR)*(Z_sim(i)-1),(Z_sim(i)-1)^2,sigma_Z^2]';
        w_sim(i)=Z_sim(i)/P*(1-ALFA)*k_sim(i)^ALFA*l_sim(i)^(-ALFA);
        r_sim(i)=Z_sim(i)/P*ALFA*k_sim(i)^(ALFA-1)*l_sim(i)^(1-ALFA)-DELTA;
    end
    figure(2)
    subplot(2,3,1)
    plot(k_sim)
    title('Capital ($k$)','Interpreter','latex')
    subplot(2,3,2)
    plot(c_sim)
    title('Consumption ($c$)','Interpreter','latex')
    subplot(2,3,3)
    plot(l_sim)
    title('Labor ($l$)','Interpreter','latex')
    subplot(2,3,4)
    plot(Z_sim)
    title('Productivity shock ($\zeta$)','Interpreter','latex')    
    subplot(2,3,5)
    plot(r_sim)
    title('Interest rate ($r$)','Interpreter','latex')
    subplot(2,3,6)
    plot(w_sim)
    title('Wage rate ($w$)','Interpreter','latex')
    print(['SIM_comp_1_approx_',num2str(approx)],'-djpeg','-r150')
    close(2)
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
                {dec_c*[1,(sym('kp')-KSTAR),(sym('Zp')-ZSTAR),(sym('kp')-KSTAR)^2,(sym('kp')-KSTAR)*(sym('Zp')-ZSTAR),(sym('Zp2')-2*sym('Zp')*ZSTAR+ZSTAR^2),sigma_Z^2]',...
                dec_l*[1,(sym('kp')-KSTAR),(sym('Zp')-ZSTAR),(sym('kp')-KSTAR)^2,(sym('kp')-KSTAR)*(sym('Zp')-ZSTAR),(sym('Zp2')-2*sym('Zp')*ZSTAR+ZSTAR^2),sigma_Z^2]'...
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
    
    
    if Loop_P
        Pl=1:0.005:1.4;
        LP=length(Pl);
        
        k_P(1:LP)=k_sim(2);
        c_P(1:LP)=c_sim(2);
        l_P(1:LP)=l_sim(2);
        Z_P(1:LP)=Z_sim(2);
        r_P(1:LP)=r_sim(2);
        w_P(1:LP)=w_sim(2);
      
        for i=2:LP
            Z_P(i)=Z_sim(i+1);
            [k_P(i),c_P(i),l_P(i),r_P(i),w_P(i)]=model1_P(Pl(i),approx,k_P(i-1),Z_P(i-1),Z_P(i),LAMBDAZ,sigma_Z);
        end
        figure(3)
        subplot(2,3,1)
        plot(Pl,k_P)
        title('Capital ($k$)','Interpreter','latex')
        subplot(2,3,2)
        plot(Pl,c_P)
        title('Consumption ($c$)','Interpreter','latex')
        subplot(2,3,3)
        plot(Pl,l_P)
        title('Labor ($l$)','Interpreter','latex')
        subplot(2,3,4)
        plot(Pl,Z_P)
        title('Productivity shock ($\zeta$)','Interpreter','latex')
        subplot(2,3,5)
        plot(Pl,r_P)
        title('Interest rate ($r$)','Interpreter','latex')
        subplot(2,3,6)
        plot(Pl,w_P)
        title('Wage rate ($w$)','Interpreter','latex')
        print(['SIM_P_comp_1_approx_',num2str(approx)],'-djpeg','-r150')
        close(3)
    
    end
    
    if Sim_P
        LP=min(T_Psims,T-1);
        
        k_P=ones([1 LP])*k_sim(2);
        c_P=ones([1 LP])*c_sim(2);
        l_P=ones([1 LP])*l_sim(2);
        Z_P=ones([1 LP])*Z_sim(2);
        r_P=ones([1 LP])*r_sim(2);
        w_P=ones([1 LP])*w_sim(2);
        sigma_P = 0.005;
        
        rng(13,'twister');
        rho_P = normrnd(0,sigma_P,[1 LP]);
        
        Pl=ones([1 LP])*G^(1/(1-LAMBDAP))*exp(rho_P(1));
        
        mu_mat = [MU,MU^2,MU^3,MU^4,MU^5];
      
        for i=2:LP
            Z_P(i)=Z_sim(i+1);
            first_Z = max(1,i-5);
            Pl(i) = Pl(i-1)^LAMBDAP*G*prod(Z_sim(first_Z:(i-1)).^mu_mat(1:min(i-1,5)))*exp(rho_P(i));
            Pl(i)
            i
            [k_P(i),c_P(i),l_P(i),r_P(i),w_P(i)]=model1_P(Pl(i),approx,k_P(i-1),Z_P(i-1),Z_P(i),LAMBDAZ,sigma_Z);
        end
        figure(3)
        subplot(2,3,1)
        plot(k_P)
        title('Capital ($k$)','Interpreter','latex')
        subplot(2,3,2)
        plot(c_P)
        title('Consumption ($c$)','Interpreter','latex')
        subplot(2,3,3)
        plot(l_P)
        title('Labor ($l$)','Interpreter','latex')
        subplot(2,3,4)
        plot(Z_P)
        title('Productivity shock ($\zeta$)','Interpreter','latex')
        subplot(2,3,5)
        plot(r_P)
        title('Interest rate ($r$)','Interpreter','latex')
        subplot(2,3,6)
        plot(Pl)
        title('Markup ($P$)','Interpreter','latex')
        print(['SIM_P_comp_5_approx_',num2str(approx)],'-djpeg','-r200')
%         close(3)
    
    end
    
end