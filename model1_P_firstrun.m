% Calls: model1.m num_eval.m  model1_ss_numeric.m gx_hx.m gxx_hxx.m gss_hss.m

function [KSTAR,LSTAR,CSTAR,stable_dkplusdk]=model1_P_firstrun(P,LAMBDAZ,sigma_Z,BETTA,DELTA,ALFA,ETA,GAMA,SIGM,G,eta,ZSTAR,k_guess,c_guess,l_guess,...
    sym_labor_supply,intertemporal_euler_ss,fx,fxp,fy,fyp,f)

[KSTAR,CSTAR,LSTAR,WSTAR,RSTAR]=model1_ss_numeric(k_guess,c_guess,l_guess,DELTA,ALFA,BETTA,G,P,ETA,GAMA,SIGM,ZSTAR,sym_labor_supply,intertemporal_euler_ss);

k=KSTAR; c=CSTAR; l=LSTAR; Z=ZSTAR;
kp=k; cp=c; lp=l; Zp=Z;
stockSTAR = G*((P-1)/P)*y_func(KSTAR,LSTAR,ZSTAR,ALFA)/((1+RSTAR) - G); stock = stockSTAR; stockp = stockSTAR;


%Obtain numerical derivatives of f
approx = 1;
num_eval

%First-order approximation
[~,hx] = gx_hx(nfy,nfx,nfyp,nfxp);
stable_dkplusdk = hx(1,1);