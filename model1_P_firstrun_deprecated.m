% Calls: model1.m num_eval.m  model1_ss_numeric.m gx_hx.m gxx_hxx.m gss_hss.m

function [k_sim_plus,c_sim,l_sim,KSTAR,LSTAR,CSTAR,stable_dkplusdk]=model1_P_firstrun(P,approx,k_sim,Z_sim,LAMBDAZ,sigma_Z,BETTA,DELTA,ALFA,ETA,GAMA,SIGM,G,eta,ZSTAR,k_guess,c_guess,l_guess,...
    sym_labor_supply,intertemporal_euler_ss,intertemporal_euler_sym,fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f,nextshock)

[KSTAR,CSTAR,LSTAR,WSTAR,RSTAR]=model1_ss_numeric(k_guess,c_guess,l_guess,DELTA,ALFA,BETTA,G,P,ETA,GAMA,SIGM,ZSTAR,sym_labor_supply,intertemporal_euler_ss);

k=KSTAR; c=CSTAR; l=LSTAR; Z=ZSTAR;
kp=k; cp=c; lp=l; Zp=Z; riskless_r_ = RSTAR;

%Obtain numerical derivatives of f
num_eval

%First-order approximation
[gx,hx] = gx_hx(nfy,nfx,nfyp,nfxp);
stable_dkplusdk = hx(1,1);

flatten = @(A) A(:);
if approx == 2
    %Second-order approximation
    [gxx,hxx] = gxx_hxx(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx); 
    [gss,hss] = gss_hss(nfx,nfxp,nfy,nfyp,nfypyp,nfypy,nfypxp,nfypx,nfyyp,nfyy,nfyxp,nfyx,nfxpyp,nfxpy,nfxpxp,nfxpx,nfxyp,nfxy,nfxxp,nfxx,hx,gx,gxx,eta);
    dec_k=[KSTAR,hx(1,:),1/2*flatten(hxx(1,:,:))',1/2*hss(1)]; 
    dec_l=[LSTAR,gx(1,:),1/2*flatten(gxx(1,:,:))',1/2*gss(1)];
    dec_c=[CSTAR,gx(2,:),1/2*flatten(gxx(2,:,:))',1/2*gss(2)];
    dec_riskless_r=[RSTAR,gx(3,:),1/2*flatten(gxx(3,:,:))',1/2*gss(3)];
else
    dec_k=[KSTAR,hx(1,:),0,0,0,0,0]; 
    dec_l=[LSTAR,gx(1,:),0,0,0,0,0];
    dec_c=[CSTAR,gx(2,:),0,0,0,0,0];   
    dec_riskless_r=[RSTAR,gx(3,:),0,0,0,0,0];
end

k_sim_plus=decision_func(dec_k,[k_sim Z_sim],[KSTAR ZSTAR],sigma_Z);
    
c_sim=decision_func(dec_c,[k_sim Z_sim],[KSTAR ZSTAR],sigma_Z);
l_sim=decision_func(dec_l,[k_sim Z_sim],[KSTAR ZSTAR],sigma_Z);
riskless_r_sim_alt = decision_func(dec_riskless_r,[k_sim Z_sim],[KSTAR ZSTAR],sigma_Z);
% c_sim = (1-DELTA) * k_sim + y_func(k_sim,l_sim,Z_sim,ALFA) - G*k_sim_plus
% fprintf('%.25g Z_sim\n',Z_sim)
E_Z = Z_sim^LAMBDAZ*exp(sigma_Z^2/2);
E_Z2 = Z_sim^LAMBDAZ*exp(2*sigma_Z^2);
E_l_sim = decision_func(dec_l,[k_sim Z_sim^LAMBDAZ*exp(sigma_Z^2/2)],[KSTAR ZSTAR],sigma_Z);

E_c_sim = dec_c * [1,k_sim-KSTAR,E_Z-ZSTAR,(k_sim-KSTAR)^2,(k_sim-KSTAR)*(E_Z-ZSTAR),(k_sim-KSTAR)*(E_Z-ZSTAR),E_Z2-2*E_Z*ZSTAR+ZSTAR^2,sigma_Z^2]';

riskless_r_sim = riskless_r_model1(k_sim_plus,Z_sim,c_sim,l_sim,dec_c,dec_l,LAMBDAZ,KSTAR,ZSTAR,sigma_Z,BETTA,G,GAMA,ETA,SIGM,intertemporal_euler_sym,DELTA,ALFA,dec_k,nextshock);
end

        


 

