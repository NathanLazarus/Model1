% Calls: model1.m num_eval.m  model1_ss_numeric.m gx_hx.m gxx_hxx.m gss_hss.m

function [k_sim,c_sim,l_sim]=model1_P(P,approx,k_sim_,Z_sim_,Z_sim,LAMBDAZ,sigma_Z,BETTA,DELTA,ALFA,ETA,GAMA,SIGM,G,eta,ZSTAR,u,force)

[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f] = model1(u);

ZSTAR = 1; %steady-state value of technology shock 

[KSTAR,CSTAR,LSTAR]=model1_ss_numeric(DELTA,ALFA,BETTA,G,P,ETA,GAMA,SIGM,ZSTAR,u);

k=KSTAR; c=CSTAR; l=LSTAR; Z=ZSTAR;
kp=k; cp=c; lp=l; Zp=Z;

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
    dec_c=[CSTAR,gx(2,:),1/2*flatten(gxx(2,:,:))',1/2*gss(2)];
else
    dec_k=[KSTAR,hx(1,:),0,0,0,0,0]; 
    dec_l=[LSTAR,gx(1,:),0,0,0,0,0];
    dec_c=[CSTAR,gx(2,:),0,0,0,0,0];   
end

k_sim=decision_func(dec_k,[k(i-1) Z(i-1)],[KSTAR ZSTAR],sigma_Z);

if force
    k_sim = k_sim_; %allows exogenous k for period 1
end
    
c_sim=decision_func(dec_c,[k(i) Z(i)],[KSTAR ZSTAR],sigma_Z);
l_sim=decision_func(dec_l,[k(i) Z(i)],[KSTAR ZSTAR],sigma_Z);

end

        


 

