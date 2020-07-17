function [KSTAR,LSTAR,CSTAR,stockSTAR,stable_dcdk,k_sim_plus,c_sim,l_sim,stock_sim]=model1_P_firstrun(P,k_sim,Z_sim,LAMBDAZ,sigma_Z,BETTA,DELTA,ALFA,ETA,GAMA,SIGM,G,eta,ZSTAR,k_guess,c_guess,l_guess,...
    order,sym_labor_supply,intertemporal_euler_ss,fx,fxp,fy,fyp,f,MathematicaOutputs)

[KSTAR,CSTAR,LSTAR,WSTAR,RSTAR]=model1_ss_numeric(k_guess,c_guess,l_guess,DELTA,ALFA,BETTA,G,P,ETA,GAMA,SIGM,ZSTAR,sym_labor_supply,intertemporal_euler_ss);
stockSTAR = G*((P-1)/P)*y_func(KSTAR,LSTAR,ZSTAR,ALFA)/((1+RSTAR) - G);

k=KSTAR; c=CSTAR; l=LSTAR; Z=ZSTAR;
kp=k; cp=c; lp=l; Zp=Z;
stock = stockSTAR; stockp = stockSTAR;


%Obtain numerical derivatives of f
approx = 1;
num_eval

%First-order approximation
[gx,~] = gx_hx(nfy,nfx,nfyp,nfxp);
stable_dcdk = gx(2,1);
[k_sim_plus,c_sim,l_sim,stock_sim] = deal(NaN);


if ~isempty(MathematicaOutputs)
    [deviation,currentPloc]=min(abs(MathematicaOutputs(:,1)-P));
    if deviation > 1e-8
        warning('No coefficients found for P = %f, using coefficients for P = %f',P,MathematicaOutputs(currentPloc,1))
    end
    [~,width] = size(MathematicaOutputs);
    currentPcoefs = mat2cell(MathematicaOutputs(currentPloc,2:width),1,[(width-1)/4,(width-1)/4,(width-1)/4,(width-1)/4]);
    h1 = currentPcoefs{1};
    g1 = currentPcoefs{2};
    g2 = currentPcoefs{3};
    g3 = currentPcoefs{4};
    
    max_order_I_can_handle = 4;
    nstates_plus_shocks = 3;
    decision_func_to_use = str2func(join(["decision_func_","order",max_order_I_can_handle,"_model1"],""));
    
    for order_to_check_for_coefs=1:max_order_I_can_handle
        if order >= order_to_check_for_coefs
            index_lb = sum(nstates_plus_shocks.^(1:(order_to_check_for_coefs-1)))+1;
            index_ub = sum(nstates_plus_shocks.^(1:order_to_check_for_coefs));
            eval(sprintf('order%1$d = [h1(index_lb:index_ub);g1(index_lb:index_ub);g2(index_lb:index_ub);g3(index_lb:index_ub)];',order_to_check_for_coefs))
        else
            eval(sprintf('order%1$d = zeros([4 nstates_plus_shocks^%1$d]);',order_to_check_for_coefs))
        end
    end

    dec_k=[KSTAR,order1(1,:),1/2*order2(1,:),1/6*order3(1,:),1/24*order4(1,:)];
    dec_stock=[stockSTAR,order1(2,:),1/2*order2(2,:),1/6*order3(2,:),1/24*order4(2,:)];
    dec_l=[LSTAR,order1(3,:),1/2*order2(3,:),1/6*order3(3,:),1/24*order4(3,:)];
    dec_c=[CSTAR,order1(4,:),1/2*order2(4,:),1/6*order3(4,:),1/24*order4(4,:)];

    k_sim_plus=decision_func_to_use(dec_k,[k_sim Z_sim],[KSTAR ZSTAR],sigma_Z);

    c_sim=decision_func_to_use(dec_c,[k_sim Z_sim],[KSTAR ZSTAR],sigma_Z);
    l_sim=decision_func_to_use(dec_l,[k_sim Z_sim],[KSTAR ZSTAR],sigma_Z);
    stock_sim=decision_func_to_use(dec_stock,[k_sim Z_sim],[KSTAR ZSTAR],sigma_Z);
end