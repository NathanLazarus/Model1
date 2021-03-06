% model1.M
function [fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx,f] = model1(u)

syms DELTA ALFA BETTA G P LAMBDAZ ETA GAMA SIGM
syms c cp l lp k kp Z Zp riskless_r_ riskless_rp_ stock stockp

up = subs(u,[c l],[cp lp]);
dupdcp = jacobian(up,cp);
dudc = jacobian(u,c);

f1 = c + G*kp - (1-DELTA) * k - y_func(k,l,Z,ALFA);
f2 = dudc - BETTA * subs(dupdcp,cp,G*cp) * big_R(kp,lp,P,Zp,ALFA,DELTA);
f3 = laborsupply(u) + w_func(k,l,P,Z,ALFA);
f4 = Zp - Z^LAMBDAZ;
f5 = dudc - BETTA * subs(dupdcp,cp,G*cp) * G * (stockp + ((P-1)/P)*y_func(kp,lp,Zp,ALFA))/stock;

f = [f1;f2;f3;f4;f5];

x = [k Z];
y = [l c stock];
xp = [kp Zp];
yp = [lp cp stockp];

[fx,fxp,fy,fyp,fypyp,fypy,fypxp,fypx,fyyp,fyy,fyxp,fyx,fxpyp,fxpy,fxpxp,fxpx,fxyp,fxy,fxxp,fxx]=anal_deriv(f,x,y,xp,yp);