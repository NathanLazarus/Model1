function riskless_r = riskless_r(k,Z,c_input,l_input,dec_c,dec_l,LAMBDAZ,KSTAR,ZSTAR,sigma_Z,BETTA,G,GAMA,ETA,SIGM,u)
riskless_r = NaN([1 length(k)]);
%intertemporal_euler = dupdcp_over_dudc(u,0);
sigma_Z = 0.000072
quadpoints = [-5.38748089, -4.60368245, -3.94476404, -3.34785457, -2.78880606,...
       -2.254974  , -1.73853771, -1.23407622, -0.73747373, -0.24534071,...
        0.24534071,  0.73747373,  1.23407622,  1.73853771,  2.254974  ,...
        2.78880606,  3.34785457,  3.94476404,  4.60368245,  5.38748089];
    
quadweights = [2.22939365e-13, 4.39934099e-10, 1.08606937e-07, 7.80255648e-06,...
       2.28338636e-04, 3.24377334e-03, 2.48105209e-02, 1.09017206e-01,...
       2.86675505e-01, 4.62243670e-01, 4.62243670e-01, 2.86675505e-01,...
       1.09017206e-01, 2.48105209e-02, 3.24377334e-03, 2.28338636e-04,...
       7.80255648e-06, 1.08606937e-07, 4.39934099e-10, 2.22939365e-13];
   
for i=2:length(k)
    E_Zp = Z(i-1)^LAMBDAZ;
    dec_c
    E_cp = decision_func(dec_c,[k(i) E_Zp],[KSTAR ZSTAR],sigma_Z)
    E_lp = decision_func(dec_l,[k(i) E_Zp],[KSTAR ZSTAR],sigma_Z)
    wrong_e_cp = 1/(BETTA*(G*E_cp/c_input(i-1))^-SIGM)
    quad = 0;
    c = c_input(i-1);
    l = l_input(i-1);
    cp = E_cp
    lp = E_lp
    for j=1:length(quadpoints)
        cp = decision_func(dec_c,[k(i) Z(i-1)^LAMBDAZ*exp(sqrt(2)*quadpoints(j)).^(sigma_Z)],[KSTAR ZSTAR],sigma_Z);
        lp = decision_func(dec_l,[k(i) Z(i-1)^LAMBDAZ*exp(sqrt(2)*quadpoints(j)).^(sigma_Z)],[KSTAR ZSTAR],sigma_Z);
        quad = quad + quadweights(j)*(1/sqrt(pi))*BETTA/eval(dupdcp_over_dudc(u,0));
    end
    right_e_cp = quad
    riskless_r(i) = quad - 1;
    fprintf('%.15g\n',wrong_e_cp)
    fprintf('%.15g\n',right_e_cp)
%     E_lp = decision_func(dec_l,[k(i) E_Zp],[KSTAR ZSTAR],sigma_Z);
%     cp = E_cp;
%     lp = E_lp;
%     c = c_input(i-1);
%     l = l_input(i-1);
%     dupdcp_over_dudc(u,0)
%     riskless_r(i) = BETTA/subs(dupdcp_over_dudc(u,0)) - 1; %,[c cp l lp G GAMA ETA SIGM],[c E_cp l E_lp G GAMA ETA SIGM]
end

dupdcp_over_dudc(u,0)