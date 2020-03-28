% Calls: model1.m num_eval.m  model1_ss_numeric.m gx_hx.m gxx_hxx.m gss_hss.m

function [k_sim,c_sim,l_sim,r_sim,w_sim,y_sim]=model1_P(P,approx,k_sim_,Z_sim_,Z_sim,LAMBDAZ,sigma_Z,BETTA,DELTA,ALFA,ETA,GAMA,G,eta,ZSTAR,force)

[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f] = model1;

ZSTAR = 1; %steady-state value of technology shock 

[KSTAR,CSTAR,LSTAR]=model1_ss_numeric(DELTA,ALFA,BETTA,G,P,ETA,GAMA,ZSTAR);

k=KSTAR; c=CSTAR; l=LSTAR; Z=ZSTAR;
kp=k; cp=c; lp=l; Zp=Z;

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

k_sim=dec_k*[1,k_sim_-KSTAR,Z_sim_-1,(k_sim_-KSTAR)^2,(k_sim_-KSTAR)*(Z_sim_-1),(Z_sim_-1)^2,sigma_Z^2]';
if force
    k_sim = k_sim_; %allows exogenous k for period 1
end
    
c_sim=dec_c*[1,k_sim-KSTAR,Z_sim-1,(k_sim-KSTAR)^2,(k_sim-KSTAR)*(Z_sim-1),(Z_sim-1)^2,sigma_Z^2]';
l_sim=dec_l*[1,k_sim-KSTAR,Z_sim-1,(k_sim-KSTAR)^2,(k_sim-KSTAR)*(Z_sim-1),(Z_sim-1)^2,sigma_Z^2]';
w_sim=Z_sim/P*(1-ALFA)*k_sim^ALFA*l_sim^(-ALFA);
r_sim=Z_sim/P*ALFA*k_sim^(ALFA-1)*l_sim^(1-ALFA)-DELTA;
y_sim=Z_sim*k_sim^ALFA*l_sim^(1-ALFA);

end

        


 

