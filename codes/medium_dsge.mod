// Chapter 4 baseline: risk premium, GDP measurement and GDP Taylor rule. Nonlinear equilibrium E1-E30, Dynare 7.1.
@#ifndef RUN_ESTIMATION
 @#define RUN_ESTIMATION = 0
@#endif
@#ifndef MH_DRAWS
 @#define MH_DRAWS = 100000
@#endif
@#ifndef RESUME_MODE
 @#define RESUME_MODE = 0
@#endif
@#ifndef EXTEND_CHAINS
 @#define EXTEND_CHAINS = 0
@#endif
@#ifndef MH_CHAINS
 @#define MH_CHAINS = 4
@#endif
var C I K u q lam N L wh w Y Yw mc Pi R G A
    ps ws Fp Hp Fw Hw Dp Dw rp zeta mp mw
    GDP GDPobs Cobs Iobs N_obs wobs Piobs Robs Lobs Kobs uobs qobs;
varexo eA eG eR eb eI epmark ewmark;
parameters beta eta chi alpha delta h phiI a1 a2 ep ew tp tw
    rhoR phipi phiy rhoA rhoG sigA sigG sigR A_ss G_ss Pi_ss Y_ss
    rhob rhoI rhop rhow sigb sigI sigp sigw;

load calibration_input.mat;
beta=data_beta; eta=1; delta=0.025; ep=6; ew=6;
h=0.60; phiI=4; tp=0.75; tw=0.75;
rhoR=0.75; phipi=1.50; phiy=0.125;
rhoA=0.5; rhoG=0.5; rhob=0.5; rhoI=0.5; rhop=0.5; rhow=0.5;
sigA=0.01; sigG=0.05; sigR=0.0025; sigb=0.0025;
sigI=0.01; sigp=0.01; sigw=0.01;
Y_ss=1; Pi_ss=data_Pi; G_ss=data_sG;
a1=1/beta-1+delta; a2=0.1;
alpha=data_sI*a1/(delta*((ep-1)/ep));
A_ss=1/((data_sI/delta)^alpha*(1/3)^(1-alpha));
chi=((ew-1)/ew)*(1-alpha)*((ep-1)/ep)/(1/3)/((1-h)*data_sC)/(1/3)^eta;

model;
 # S=phiI/2*(I/I(-1)-1)^2;
 # dS=phiI*(I/I(-1)-1);
 # dSp=phiI*(I(+1)/I-1);
 # au=a1*(u-1)+a2/2*(u-1)^2;
 # aup=a1*(u(+1)-1)+a2/2*(u(+1)-1)^2;
 # disc=beta*lam(+1)/lam;
 # R_ss=Pi_ss/beta;
 lam=1/(C-h*C(-1));
 chi*N^eta=lam*wh;
 1=disc*R*exp(rp)/Pi(+1);
 Yw=A*(u*K(-1))^alpha*L^(1-alpha);
 w*L=(1-alpha)*mc*Yw;
 alpha*mc*Yw=u*(a1+a2*(u-1))*K(-1);
 q=disc*(alpha*mc(+1)*Yw(+1)/K-aup+(1-delta)*q(+1));
 1=q*zeta*(1-S-(I/I(-1))*dS)+disc*q(+1)*zeta(+1)*dSp*(I(+1)/I)^2;
 K=(1-delta)*K(-1)+zeta*(1-S)*I;
 ps=ep/(ep-1)*Fp/Hp;
 Fp=exp(mp)*mc*Y+tp*disc*(Pi(+1)/Pi_ss)^ep*Fp(+1);
 Hp=Y+tp*disc*(Pi(+1)/Pi_ss)^(ep-1)*Hp(+1);
 1=(1-tp)*ps^(1-ep)+tp*(Pi_ss/Pi)^(1-ep);
 Dp=(1-tp)*ps^(-ep)+tp*(Pi_ss/Pi)^(-ep)*Dp(-1);
 ws=ew/(ew-1)*Fw/Hw;
 Fw=exp(mw)*wh*w^ew*L+tw*disc*(Pi(+1)/Pi_ss)^ew*Fw(+1);
 Hw=w^ew*L+tw*disc*(Pi(+1)/Pi_ss)^(ew-1)*Hw(+1);
 w^(1-ew)=(1-tw)*ws^(1-ew)+tw*(Pi_ss/Pi*w(-1))^(1-ew);
 Dw=(1-tw)*(ws/w)^(-ew)+tw*(Pi_ss/Pi*w(-1)/w)^(-ew)*Dw(-1);
 N=Dw*L;
 Yw=Dp*Y;
 Y=C+I+G+au*K(-1);
 log(G/G_ss)=rhoG*log(G(-1)/G_ss)+sigG*eG;
 log(A/A_ss)=rhoA*log(A(-1)/A_ss)+sigA*eA;
 log(R/R_ss)=rhoR*log(R(-1)/R_ss)+(1-rhoR)*(phipi*log(Pi/Pi_ss)+phiy*log(GDP/steady_state(GDP)))+sigR*eR;
 rp=rhob*rp(-1)+sigb*eb;
 log(zeta)=rhoI*log(zeta(-1))+sigI*eI;
 mp=rhop*mp(-1)+sigp*epmark;
 mw=rhow*mw(-1)+sigw*ewmark;
 GDP=C+I+G;
 GDPobs=100*log(GDP/steady_state(GDP));
 Cobs=100*log(C/steady_state(C));
 Iobs=100*log(I/steady_state(I));
 N_obs=100*log(N/steady_state(N));
 wobs=100*log(w/steady_state(w));
 Piobs=100*log(Pi/Pi_ss);
 Robs=100*log(R/R_ss);
 Lobs=100*log(L/steady_state(L));
 Kobs=100*log(K/steady_state(K));
 uobs=100*log(u);
 qobs=100*log(q);
end;

steady_state_model;
 mc=(ep-1)/ep;
 a1=1/beta-1+delta;
 Y=Y_ss; Yw=Y; G=G_ss;
 K=alpha*mc*Y/a1; I=delta*K; C=Y-I-G;
 A=A_ss; L=(Y/(A*K^alpha))^(1/(1-alpha)); N=L;
 w=(1-alpha)*mc*Y/L; wh=(ew-1)/ew*w;
 lam=1/((1-h)*C); chi=lam*wh/N^eta;
 Pi=Pi_ss; R=Pi_ss/beta; u=1; q=1; ps=1; ws=w; Dp=1; Dw=1;
 rp=0; zeta=1; mp=0; mw=0;
 Fp=mc*Y/(1-beta*tp); Hp=Y/(1-beta*tp);
 Fw=wh*w^ew*L/(1-beta*tw); Hw=w^ew*L/(1-beta*tw);
 GDP=C+I+G; GDPobs=0; Cobs=0; Iobs=0; N_obs=0; wobs=0; Piobs=0; Robs=0;
 Lobs=0; Kobs=0; uobs=0; qobs=0;
end;
shocks;
 var eA; stderr 1;
 var eG; stderr 1;
 var eR; stderr 1;
 var eb; stderr 1;
 var eI; stderr 1;
 var epmark; stderr 1;
 var ewmark; stderr 1;
end;
steady;
resid;
check;
model_diagnostics;

@#if RUN_ESTIMATION
estimated_params;
 h, beta_pdf, 0.60, 0.15;
 phiI, gamma_pdf, 4, 2;
 a2, gamma_pdf, 0.1, 0.1;
 tp, beta_pdf, 0.75, 0.10;
 tw, beta_pdf, 0.75, 0.10;
 rhoR, beta_pdf, 0.75, 0.10;
 phipi, gamma_pdf, 1.50, 0.25;
 phiy, gamma_pdf, 0.125, 0.075;
 rhoA, beta_pdf, 0.50, 0.20;
 rhoG, beta_pdf, 0.50, 0.20;
 rhob, beta_pdf, 0.50, 0.20;
 rhoI, beta_pdf, 0.50, 0.20;
 rhop, beta_pdf, 0.50, 0.20;
 rhow, beta_pdf, 0.50, 0.20;
 sigA, inv_gamma1_pdf, 0.010, 0.005;
 sigG, inv_gamma1_pdf, 0.050, 0.025;
 sigR, inv_gamma1_pdf, 0.0025, 0.00125;
 sigb, inv_gamma1_pdf, 0.0025, 0.00125;
 sigI, inv_gamma1_pdf, 0.010, 0.005;
 sigp, inv_gamma1_pdf, 0.010, 0.005;
 sigw, inv_gamma1_pdf, 0.010, 0.005;
end;
varobs GDPobs Cobs Iobs N_obs wobs Piobs Robs;
set_dynare_seed(20261110);
estimation(datafile='dsge_observables.mat', first_obs=1, nobs=94,
 presample=0, order=1, prefilter=0, lik_init=1,
@#if RESUME_MODE
 mode_compute=0, mode_file='medium_dsge/Output/medium_dsge_mode.mat',
 posterior_sampler_options=('scale_file','medium_dsge/Output/medium_dsge_optimal_mh_scale_parameter.mat'),
@#else
 mode_compute=6,
@#endif
 mh_replic=@{MH_DRAWS}, mh_nblocks=@{MH_CHAINS}, mh_drop=0.5,
@#if EXTEND_CHAINS
 load_mh_file,
@#endif
 mh_conf_sig=0.90, nograph, nodisplay)
 GDPobs Cobs Iobs N_obs wobs Piobs Robs Lobs Kobs uobs qobs;
@#if MH_DRAWS > 0
shock_decomposition(parameter_set=posterior_mean, nograph)
 GDPobs Cobs Iobs N_obs Piobs Robs;
stoch_simul(order=1, irf=40, nograph, noprint)
 GDPobs Cobs Iobs N_obs Piobs Robs wobs Kobs uobs qobs;
@#endif
@#else
stoch_simul(order=1, irf=40, nograph, noprint)
 GDPobs Cobs Iobs N_obs Piobs Robs wobs Kobs uobs qobs;
@#endif
