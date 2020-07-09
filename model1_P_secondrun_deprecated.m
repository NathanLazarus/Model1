% Calls: model1.m num_eval.m  model1_ss_numeric.m gx_hx.m gxx_hxx.m gss_hss.m

function [k_sim_plus,c_sim,l_sim]=model1_P_secondrun(P,k_sim,Z_sim,LAMBDAZ,sigma_Z,BETTA,DELTA,ALFA,ETA,GAMA,SIGM,G,eta,ZSTAR,k_guess,c_guess,l_guess,...
    sym_labor_supply,intertemporal_euler_ss,MathematicaOutputs)

[KSTAR,CSTAR,LSTAR,WSTAR,RSTAR]=model1_ss_numeric(k_guess,c_guess,l_guess,DELTA,ALFA,BETTA,G,P,ETA,GAMA,SIGM,ZSTAR,sym_labor_supply,intertemporal_euler_ss);

k=KSTAR; c=CSTAR; l=LSTAR; Z=ZSTAR;
kp=k; cp=c; lp=l; Zp=Z; riskless_r_ = RSTAR;
[~,currentPloc]=min(abs(MathematicaOutputs(:,1)-P));
[~,width] = size(MathematicaOutputs);
currentPcoefs = mat2cell(MathematicaOutputs(currentPloc,2:width),1,[(width-1)/3,(width-1)/3,(width-1)/3]);
h1 = currentPcoefs{1};
g1 = currentPcoefs{2};
g2 = currentPcoefs{3};

dec_k=[KSTAR,h1(1:3),1/2*h1(4:12),1/3*h1(13:39)];
dec_l=[LSTAR,g1(1:3),1/2*g1(4:12),1/3*g1(13:39)];
dec_c=[CSTAR,g2(1:3),1/2*g2(4:12),1/3*g2(13:39)];
if currentPloc == 1
    dec_l
end

k_sim_plus=decision_func_thirdorder_model1(dec_k,[k_sim Z_sim],[KSTAR ZSTAR],sigma_Z);
    
c_sim=decision_func_thirdorder_model1(dec_c,[k_sim Z_sim],[KSTAR ZSTAR],sigma_Z);
l_sim=decision_func_thirdorder_model1(dec_l,[k_sim Z_sim],[KSTAR ZSTAR],sigma_Z);

end

        


 

