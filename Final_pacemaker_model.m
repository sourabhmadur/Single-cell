

  clear;            
  clc;                 
  %close all;
  %------------------------------------------------------------------------
  
  Ts=8;               % Total time of stimulation - 60 seconds
  dt=0.0001;           % Step duration for integration using Eulers Method
  
  %------------------------------------------------------------------------
  
  T = 310.0;           % Kelvin
  T_exp=297;           % Kelvin, temperature of experiments
  F = 96.4846;         % microcoulomb_per_nanomole
  R = 8.3144;          % nanojoule_per_nanomole_kelvin
  FoRT = F/(R*T);      % per_millivolt
  RToF = (R*T)/F;      % millivolt
  Cm = 0.025;          % nanofarad, cell capacitance

  %------------------------------------------------------------------------
  %                   Cell Geometry and buffering
  %------------------------------------------------------------------------
 
  Vol = 1.0e-12;      % in liters. 1000 cubic microns=10^-12 litre
  P_PU = 0.001;       % dimensionless, fraction of caveolar space
  Pmito = 0.12871;    % dimensionless, fraction of mitochondrial, double of F&K because of big mitochondria in ICC
  Per = 0.10;         % dimensionless, fraction of Er space
  P_cyto=0.7;         % dimensionless, fraction of cytoplasmic area available for Ca
  V_PU = Vol*P_PU;    % l
  Vmito = Vol*Pmito;  % l
  Ver = Vol*Per;      % l
  Vcyto = Vol*P_cyto; % l
  fc = 0.01;          % dimensionless, buffering in cytoplasm
  fm = 0.0003;        % dimensionless, buffering in Mitochondria
  fe = 0.01;          % dimensionless, buffering in ER

 %-------------------------------------------------------------------------
 %                   Ionic concentrations (constants)
 %-------------------------------------------------------------------------
  
  K_i  = 120;
  Na_i = 30;
  Ca_o = 2.5;
  K_o = 7;
  Na_o = 137;
  Cl_i= 88;
  Cl_o= 134;
 IP3 =0.0006;        % millimolar : baseline IP3 concentration
 %IP3 =0.000;        % millimolar : baseline IP3 concentration

  %------------------------------------------------------------------------
  %                     Constants from M&K 1997
  %------------------------------------------------------------------------
  
  % M&K 1997 Table 1
  deltapH = -0.4;       % dimensionless
  Cmito = 0.006995;     % millifarads
  
  % M&K 1997 Table 2
  K_res = 1.35e018;     % dimensionless
  r1 = 2.077e-018;      % dimensionless
  r2 = 1.728e-009;      % dimensionless
  r3 = 1.059e-026;      % dimensionless
  ra = 6.394e-010;      % per_second
  rb = 1.762e-013;      % per_second
  rc1 = 2.656e-019;     % per_second
  rc2 = 8.632e-027;     % per_second
  deltaPsi_B = 50.0;    % millivolt
  g = 0.85;             % dimensionless
  
  % M&K 1997 Table 3
  K_F1 = 1.71e009;      % millimolar :
                                     
  Pi_m = 20;            % millimolar
  p1 = 1.346e-008;      % dimensionless
  p2 = 7.739e-007;      % dimensionless
  p3 = 6.65e-015;       % dimensionless
  pa = 1.656e-005;      % per_second
  pb = 3.373e-007;      % per_second
  pc1 = 9.651e-014;     % per_second
  pc2 = 4.845e-019;     % per_second
  
  % M&K 1997 Table 4
  frac   = 0.5;         % dimensionless : (=f)
 
  % M&K 1997 Table 5
  K_act = 0.00038;      % millimolar
  na = 2.8;             % dimensionless
  deltaPsi_star = 91;   % millivolt
  K_Na = 9.4;           % millimolar
  K_Ca = 0.003;         % millimolar

  %------------------------------------------------------------------------
  %                  Constants from M&K 1998a
  %------------------------------------------------------------------------
  % M&K 1998a Table 1
  K_trans = 0.006;      % millimolar
  L = 50;               % dimensionless
  b = 0.5;              % dimensionless

  % M&K 1998a Table 3
  beta_max = 2.055;     % per_second
  beta1 = 1.66;         % per_millimolar
  beta2 = 0.0249;       % per_millimolar :
                        % Missing from M&K 1998a Equation (10)
  beta3 = 4.0;          % per_millimolar
  beta4 = 2.83;         % per_millimolar
  beta5 = 1.3;          % per_millimolar
  beta6 = 2.66;         % per_millimolar
  beta7 = 0.16;         % per_millimolar
  KCa_PDH = 0.00005;    % millimolar
  u1 = 15;              % dimensionless
  u2 = 1.1;             % dimensionless

  % Other
  n = 2;                % dimensionless : (Na:Ca)
  K_Glc = 8.7;          % millimolar
  nhyd = 2.7;           % dimensionless
  K_hyd = 0.05125;      % per_second : Missing from text
  
  %------------------------------------------------------------------------
  %             ER Ca2+ handling parameters
  %------------------------------------------------------------------------
  
  J_ERleak = 1.666667;  % per_second
  Jmax_IP3 = 50000;     % per_second
  d_IP3 = 0.00025;      % millimolar
  d_ACT = 0.001;        % millimolar
  d_INH = 0.0014;       % millimolar
  tauh = 4.0;           % seconds
  Jmax_serca = 1.8333;  % millimolar_per_second
  k_serca = 0.00042;    % millimolar

  % Mitochondrial Ca2+ handling parameters
  conc = 0.001;        % millimolar : correcting units of Juni
  Jmax_uni = 5000;     % per_second
  Jmax_NaCa = 0.05;    % millimolar_per_second
  J_max_leak = 0.01;   % per second

  % Mitochondrial respiration parameters
  rho_res = 0.4;        % millimolar
  rho_F1 = 0.7;         % millimolar
  g_H = 0.0033333;      % millimolar_per_second_per_millivolt
  J_red_basal = 0.3333; % millimolar_per_second
  Jmax_ANT = 15;        % millimolar_per_second
  J_hyd_max = 0.037625; % millimolar_per_second
  Glc = 1.0;            % millimolar : Varies ~1-10

  %------------------------------------------------------------------------
  % PMCA parameters
  
  J_max_PMCA = 0.088464; % mM/s,
  k_PMCA = 0.000298;     % millimolar

  % Temperature Corrections
  Q10Ca=2.1;
  Q10K=1.5;
  Q10Na=2.45;
  
  Tcorrection_Ca=(Q10Ca^((T-T_exp)/10));
  Tcorrection_K=(Q10K^((T-T_exp)/10));
  Tcorrection_Na=(Q10Na^((T-T_exp)/10));
  T_correction_BK=1.1*(T-T_exp);

%   % Maximal conductances in nS/cm^2
%   G_ERG =   2.5;    %2.5;
%   G_Ca_VDDR = 3;      %3;
%   G_Ca_Ltype= 2 ;    % 2;
%   G_Kv11 =   6.3  ;      %6.3;
%   G_bk =      23+T_correction_BK  ;    %23+T_correction_BK;
%   %G_bk = 0;
%   G_Kb =  0.15  ;    %0.15;
%   G_Na =  20 ;    %20;
%   G_Cacl =  10.1;    %10.1;
%   G_NSCC = 12.15;     %12.15;
%   g_Ano1 =  20;      %20;           % nS

  G_ERG =   2.5;    %2.5;
  G_Ca_VDDR = 3;      %3;
  G_Ca_Ltype= 2 ;    % 2;
  G_Kv11 =   6.3  ;      %6.3;
  G_bk =  23+T_correction_BK;    %23+T_correction_BK;
  G_Kb = 0.15  ;    %0.15;
  G_Na =  20 ;    %20;
  G_Cacl = 10.1;    %10.1;
  G_NSCC = 12.15;     %12.15;
  g_Ano1 =  20;      %20;           % nS

  
  
  % Other parameters of ion channels
  hBK = 2;      % hill coefficient for BK channels
  kbk = -17;    % slope factor for the BK current
  Caset= 0.001;

  KClCa=140e-6;
  hClca=3;
  K_NSCC=0.0745*1.0e-3;
  hNSCC=-85;

  %Time constants
  tau_d_VDDR=(6)*1.0e-3*Tcorrection_Ca;
  tau_f_VDDR=(40)*1.0e-3*Tcorrection_Ca;
  tau_d_Ltype=1*1.0e-3*Tcorrection_Ca;
  tau_f_Ltype=86*1.0e-3*Tcorrection_Ca;
  tau_f_Ca_Ltype=2.0*1.0e-3*Tcorrection_Ca;
  tau_d_kv11=5.0*1.0e-3*Tcorrection_K;
  tau_f_kv11=5*1.0e-3*Tcorrection_K;
  tau_ERG=3*1.0e-3*Tcorrection_K;
  tau_d_Na=3*1.0e-3*Tcorrection_Na;
  tau_f_Na=1.6*1.0e-3*Tcorrection_Na;
  tau_NSCC=350*1.0e-3;
  tau_act_Cacl=30*1.0e-3;

  %------------------------------------------------------------------------
  % Nernst potentials
  %------------------------------------------------------------------------
  E_Na=(R*T/F)*log(Na_o/Na_i);
  E_k=(R*T/F)*log (K_o/K_i);

  E_cl=(R*T/F)*log (Cl_i/Cl_o);
  E_NSCC=0;
  
  % Fixed concentrations
  total_NAD_m = 8.0;   % millimolar : total = NADH + NAD+
  total_ANP_m = 12.0;  % millimolar : total = ATP_m + ADP_m
  total_ANP_i = 2.0;   % millimolar : total = ATP_i + ADP_i
 
     
  Ca_cyto = zeros(1,Ts/dt);               % millimolar
  
  vclamp = -70;
%   Vm = vclamp*ones(1,Ts/dt);
  Vm = zeros(1,Ts/dt);                    % mV
  ADP_i = zeros(1,Ts/dt);                 % millimolar
  ADP_m = zeros(1,Ts/dt);                 % millimolar
  Ca_ER = zeros(1,Ts/dt);                 % millimolar
  Ca_PU = zeros(1,Ts/dt);                 % millimolar
  Ca_m = zeros(1,Ts/dt);                  % millimolar
  NADH_m = zeros(1,Ts/dt);                % millimolar
  deltaPsi = zeros(1,Ts/dt);              % mV
  h = zeros(1,Ts/dt);                     % dimensionless
  act_Cacl = zeros(1,Ts/dt);              % dimensionless
  d_ERG = zeros(1,Ts/dt);                 % dimensionless
  d_Ltype = zeros(1,Ts/dt);               % dimensionless
  P_openNSCC = zeros(1,Ts/dt);            % dimensionless
  d_Na = zeros(1,Ts/dt);                  % dimensionless
  d_VDDR = zeros(1,Ts/dt);                % dimensionless
  d_kv11 = zeros(1,Ts/dt);                % dimensionless
  f_Ltype = zeros(1,Ts/dt);               % dimensionless
  f_Na = zeros(1,Ts/dt);                  % dimensionless
  f_VDDR = zeros(1,Ts/dt);                % dimensionless
  f_ca_Ltype = zeros(1,Ts/dt);            % dimensionless
  f_kv11 = zeros(1,Ts/dt);                % dimensionless
  
  I_ca_VDDR  = zeros(1,Ts/dt); 
  I_ca_Ltype = zeros(1,Ts/dt); 
  I_NSCC     = zeros(1,Ts/dt); 
  I_kv11     = zeros(1,Ts/dt); 
  I_ERG      = zeros(1,Ts/dt); 
  I_BK       = zeros(1,Ts/dt); 
  I_na       = zeros(1,Ts/dt); 
  I_Kb       = zeros(1,Ts/dt); 
  I_CaCl     = zeros(1,Ts/dt);
  I_pmca     = zeros(1,Ts/dt);
  
  Ca_cyto(1) = 0.00000993087;   % millimolar
  Vm(1) = vclamp;
  %Vm(1) = -70.0;                % mV
  ADP_i(1) = 0.0077282;         % millimolar
  ADP_m(1) = 2.60093454;        % millimolar
  Ca_ER(1) = 0.007299;          % millimolar
  Ca_PU(1) = 0.0000902;         % millimolar
  Ca_m(1) = 0.000136;           % millimolar
  NADH_m(1) = 0.101476;         % millimolar
  deltaPsi(1) = 164.000044;     % mV
  h(1) = 0.9397;                % dimensionless
  act_Cacl(1) = 0.0;            % dimensionless
  d_ERG(1) = 0.0;               % dimensionless
  d_Ltype(1) = 0.0;             % dimensionless
  P_openNSCC(1) = 0.0;          % dimensionless
  d_Na(1) = 0.0;                % dimensionless
  d_VDDR(1) = 0.0;              % dimensionless
  d_kv11(1) = 0.0;              % dimensionless
  f_Ltype(1) = 1.0;             % dimensionless
  f_Na(1) = 1.0;                % dimensionless
  f_VDDR(1) = 1.0;              % dimensionless
  f_ca_Ltype(1) = 1.0;          % dimensionless
  f_kv11(1) = 1.0;              % dimensionless
  
  
  %inserting ISOC and ANO1 mechanism
  
  % Ano1 model parameters

%   g_Ano1      =   20          ;           % nS
  EC_50_i     =   1.39e-3     ;           % mM
  k_C         =   0.01248     ;           % mV^-1
  V_h         =   -100        ;           % mV
  k_V         =   0.0156      ;           % mV^-1
  
  % Cytosolic Calcium flux parameters

  K_SOC       =   200         ;           % 200 uM 

  % Ano1-SOC localisation parameters

  D_c         =   250         ;           % uM^2 s^-1
  D_m         =   75          ;           % um^2 s^-1
  K_m         =   1           ;           % uM
  B_m         =   50          ;           % uM
  r           =   0.05        ;           % uM
  v_scale     =   10e15       ;           % uM^3L^-1
  n_SOC       =   50          ;
  
  g_SOC       =   0.1         ;           % nS

  V_SOC       =   5000        ;           % ???????????? check
 
  O_Ano1      =   zeros(1,Ts/dt) ;
  I_SOC       =   zeros(1,Ts/dt) ;
  I_Ano1      =   zeros(1,Ts/dt) ;
  tau_Ano1    =   zeros(1,Ts/dt) ;
  E_ca=(R*T/(2.0*F))*log (Ca_o/Ca_cyto(1));
  
  v           =   1.0e-12       ;           % L   volume of cell
  v_CYTO      =   0.7e-12       ;           % L   volume of cytosol
  O_SOC       =   (K_SOC^8)/((K_SOC^8)+((Ca_ER(1)/1000)^8));
  I_SOC(1)    =   g_SOC * O_SOC * (Vm(1) - E_ca);         %pA
  J_SOC       =   (1e-9) * (-I_SOC(1)/(2*F*v_CYTO));
  
%   sigma       =   (v_scale*J_SOC*v_CYTO)/n_SOC;
%   phi_m       =   D_m*B_m*K_m ;
%   C2          =   (D_c*(Ca_cyto(1)/1000))-(phi_m/(K_m+(Ca_cyto(1)/1000)));
%   Ca_Ano1     =   ((-D_c*K_m)+(sigma/2*3.14*r)+C2+sqrt((((D_c*K_m)+(sigma/(2*3.14*r))+(C2))^2)+(4*D_c*phi_m)))/2*D_c;

  % ANO1 model
  EC_50       =   EC_50_i*exp(-k_C*Vm(1)); 
  O_Ano1(1)   =   1/((1+exp((V_h-Vm(1))*k_V))*(1+((EC_50/Ca_PU(1))^2)));

  
  for i=(1:(Ts/dt-1))
  
  E_ca=(R*T/(2.0*F))*log (Ca_o/Ca_cyto(i));
  
  %inserting I_SOC and ANO1
  O_SOC       =   (K_SOC^8)/((K_SOC^8)+(Ca_ER(i)^8));
  I_SOC(i)    =   g_SOC * O_SOC * (Vm(i) - E_ca);         %pA
  J_SOC       =   (1e-9) * (-I_SOC(i)/(2*F*v_CYTO));
%   
%   sigma       =   (v_scale*J_SOC*v_CYTO)/n_SOC;
%   phi_m       =   D_m*B_m*K_m ;
%   C2          =   (D_c*(Ca_cyto(i)/1000))-(phi_m/(K_m+(Ca_cyto(i)/1000)));
%   Ca_Ano1     =   ((-D_c*K_m)+(sigma/2*3.14*r)+C2+sqrt((((D_c*K_m)+(sigma/(2*3.14*r))+(C2))^2)+(4*D_c*phi_m)))/2*D_c;

  % ANO1 model
  EC_50       =   EC_50_i*exp(-k_C*Vm(i)); 
  O_Ano1_inf  =   1/((1+exp((V_h-Vm(i))*k_V))*(1+((EC_50/Ca_PU(i))^2)));
  
  t1          =   0.08163*exp(-0.57*Ca_PU(i));
  t2          =   0.07617*exp(-0.05374*Ca_PU(i));
  t3          =   70.3*exp(0.153*Ca_PU(i));
  tau_Ano1(i) =   t1 +(t2*exp(Vm(i)/t3));
  dO_Ano1     =   dt*((O_Ano1_inf-O_Ano1(i))/tau_Ano1(i));
  O_Ano1(i+1) =   dO_Ano1+O_Ano1(i);
  
  % Ano1 Ca2+ activated Cl channel
  I_Ano1(i)   =   g_Ano1*O_Ano1(i)*(Vm(i)-E_cl);                    % done
 
  % Initial concentrations (all in millimolar)
  
  NAD_m = total_NAD_m-NADH_m(i);
  ATP_m = total_ANP_m - ADP_m(i);
  ADP_mfree = 0.8*ADP_m(i);
  ADP3_m = 0.45*ADP_mfree;
  ATP4_m = 0.05*ATP_m;
  ATP_i = total_ANP_i - ADP_i(i);
  ADP_ifree = 0.3*ADP_i(i);
  ADP3_i = 0.45*ADP_ifree;
  MgADP_i = 0.55*ADP_ifree;
  ATP4_i = 0.05*ATP_i;
  
  %------------------------------------------------------------------------
  %                 Gates steady states
  %------------------------------------------------------------------------
  
  d_inf_Ltype=1.0/(1.0+exp(-(17.0+Vm(i))/4.3));
  f_inf_Ltype=1.0/(1.0+exp((43.0+Vm(i))/8.9));
  f_Ca_Ltype_INF=1.0-(1.0/(1.0+exp(-((Ca_cyto(i)-0.0001)-0.000214)/0.0000131)));
  
  d_inf_VDDR=1.0/(1.0+exp(-(26.0+Vm(i))/6.0));
  f_inf_VDDR=1.0/(1.0+exp((66.0+Vm(i))/6.0));
  
  p_open_BK=1.0/(1.0+exp((Vm(i)/kbk)-hBK*log(Ca_cyto(i)/Caset)));
  
  d_ERG_inf=0.8/(1+exp(-(20+Vm(i))/1.8))+0.2;
  
  d_Na_inf=1/(1+exp(-(47+Vm(i))/4.8));
  f_Na_inf=1/(1+exp((78+Vm(i))/(7)));
  
  d_kv11_inf=1/(1+exp(-(25+Vm(i))/7.7));
  f_kv11_inf=0.5/(1+exp((44.8+Vm(i))/4.4))+0.5;
  
  act_Cacl_inf=1.0/(1.0+((KClCa/Ca_cyto(i))^hClca));
  
  P_openNSCC_inf=1.0/(1.0+((K_NSCC/Ca_PU(i))^hNSCC));


  % ER Ca2+ handling components--------------------------------------------
  % Ca release from the SR via IP3-----------------------------------------
  J_ERout = (Jmax_IP3*((IP3/(IP3+d_IP3))^3)*((Ca_PU(i)/(Ca_PU(i)+d_ACT))^3)*(h(i)^3)+J_ERleak)*(Ca_ER(i)-Ca_PU(i));
  %J_ERout = 0;

  % Ca uptake into the ER via SERCA pumps------------
  J_SERCA = Jmax_serca*((Ca_PU(i)^2)/((k_serca^2)+(Ca_PU(i)^2)));
  %J_SERCA = 0;
  % Mitochondrial Ca2+ handling components---------------------------------
  % Ca influx through uniporter--------------------------------------------
  MWC = conc*((Ca_PU(i)/K_trans)*((1+Ca_PU(i)/K_trans)^3))/(((1+Ca_PU(i)/K_trans)^4)+(L/((1+Ca_PU(i)/K_act)^na)));
  J_uni = Jmax_uni*((MWC-(Ca_m(i)*exp(-2*FoRT*(deltaPsi(i)-deltaPsi_star))))*(2*FoRT*(deltaPsi(i)-deltaPsi_star))/(1-exp(-2*FoRT*(deltaPsi(i)-deltaPsi_star))));
  %J_uni =0;
  % Ca efflux through NaCa exchanger------------------
  J_NaCa = Jmax_NaCa*(exp(b*FoRT*(deltaPsi(i)-deltaPsi_star))/(((1+((K_Na/Na_i)^n)))*(1+(K_Ca/Ca_m(i)))));
  %J_NaCa = 0;
  % PMCA ---------------------------------
  J_PMCA = J_max_PMCA*(1.0/(1.0+(k_PMCA/Ca_cyto(i))));

  % leakage from PU to cyto
  J_leak=J_max_leak*(Ca_PU(i)-Ca_cyto(i));
  %J_leak = 100*J_leak;
  %------------------------------------------------------------------------
  %         Mitochondrial components that do not transport Ca2+
  %------------------------------------------------------------------------

  % Oxidation of NADH to NAD+------------------------
  % M&K 1997 Equation (4)
  A_res = RToF*log(K_res*(sqrt(NADH_m(i))/sqrt(NAD_m)));

  % Respiration Driven NADH Oxidation (Oxygen consumption rate)
  % M&K 1997 Equation (5)
  J_o = 0.5*rho_res*((ra*(10.0^(6.0*deltapH))+rc1*exp(6.0*FoRT*deltaPsi_B))*exp(FoRT*A_res)-ra*exp(g*6.0*FoRT*deltaPsi(i))+rc2*exp(FoRT*A_res)*exp(g*6.0*FoRT*deltaPsi(i)))/((1.0+r1*exp(FoRT*A_res))*exp(6.0*FoRT*deltaPsi_B)+(r2+r3*exp(FoRT*A_res))*exp(g*6.0*FoRT*deltaPsi(i)));

  % Respiration Driven Proton Efflux
  % M&K 1997 Equation (6)
  J_Hres = 0.661*6.0*rho_res*(ra*(10.0^(6.0*deltapH))*exp(FoRT*A_res)+rb*(10.0^(6.0*deltapH))-(ra+rb)*exp(g*6.0*FoRT*deltaPsi(i)))/((1.0+r1*exp(FoRT*A_res))*exp(6.0*FoRT*deltaPsi_B)+(r2+r3*exp(FoRT*A_res))*exp(g*6.0*FoRT*deltaPsi(i)));

  % Reduction of NAD+ to NADH----------------------------------------------
  % D-glucose utilisation rate = rate of glucokinase
  % M&K 1998a Equation (10)
  J_glyTotal = (beta_max*(1+(beta1*Glc))*beta2*Glc*ATP_i)/(1+(beta3*ATP_i)+((1+beta4*ATP_i)*beta5*Glc)+((1+beta6*ATP_i)*beta7*Glc));

  % Fraction of activated pyruvate dehydrogenase
  % M&K 1998a Equation (18)
  f_PDHa = 1/(1+(u2*(1+(u1/((1+(Ca_m(i)/KCa_PDH))^2)))));

  % Mitochondrial reducing equivalents from D-glucose
  % M&K 1998a Equation (21) : Coefficient miscalculated in paper = 7.36
  J_red = J_red_basal+6.3944*f_PDHa*J_glyTotal;

  % ATP production via TCA phosphorylation---------------------------------
  % M&K 1998a Equation (21)
  J_pTCA = (J_red_basal/3)+(0.84*f_PDHa*J_glyTotal);

  % ATP production via oxidative phosphorylation---------------------------
  % M&K 1997 Equation (12)
  A_F1 = RToF*log(K_F1*ATP_m/(ADP_mfree*Pi_m));

  % ATP Production via F1F0-ATPase
  % M&K 1997 Equation (13)
  J_pF1 = -rho_F1*(((pa*(10^(3*deltapH))+pc1*exp(3*FoRT*deltaPsi_B))*exp(FoRT*A_F1)-pa*exp(3*FoRT*deltaPsi(i))+pc2*exp(FoRT*A_F1)*exp(3*FoRT*deltaPsi(i)))/((1+p1*exp(FoRT*A_F1))*exp(3*FoRT*deltaPsi_B)+(p2+p3*exp(FoRT*A_F1))*exp(3*FoRT*deltaPsi(i))));

  % Proton Uptake due to F1F0-ATPase activity
  % M&K 1997 Equation (14)
  J_HF1 = -3.0*rho_F1*((pa*(10^(3*deltapH))*exp(FoRT*A_F1)+pb*(10^(3*deltapH))-(pa+pb)*exp(3*FoRT*deltaPsi(i)))/((1+p1*exp(FoRT*A_F1))*exp(3*FoRT*deltaPsi_B)+(p2+p3*exp(FoRT*A_F1))*exp(3*FoRT*deltaPsi(i))));

  %-ADP(in)/ATP(out) Exchange through the ANT------------------------------
  % M&K 1997 EquationICC_model (16)
  J_ANT = Jmax_ANT*((1-(((ATP4_i*ADP3_m)/(ADP3_i*ATP4_m))*exp(-FoRT*deltaPsi(i))))/((1+((ATP4_i/ADP3_i)*exp(-frac*FoRT*deltaPsi(i))))*(1+(ADP3_m/ATP4_m))));

  %-Proton leak into the mitochondria--------------------------------------
  % Proton Motive Force
  PMF = deltaPsi(i)-2.303*RToF*deltapH;

  % Proton Leakage
  % M&K 1997 Equation (7)
  J_Hleak = g_H*PMF;

  %-Cytosolic components---------------------------------------------------
  %-Cytosolic ATP phosphorylation------------------------------------------
  % M&K 1998a near Equation (34)
  J_pGly = 0.15*J_glyTotal;

  %-Cytosolic ATP hydrolysis-----------------------------------------------
  % Steady state hydrolysis rate
  % M&K 1998a Equation (37)
  J_hydSS = J_hyd_max/(1+((K_Glc/Glc)^nhyd));

  % Rate of cytosolic ATP hydrolysis
  % M&K 1998a Equation (35)
  J_hyd = (K_hyd*ATP_i)+J_hydSS;


  %------------------------------------------------------------------------
  % Computation of currents, in pA since they are nS*mV
  %------------------------------------------------------------------------
  I_ca_VDDR(i)=G_Ca_VDDR*d_VDDR(i)*f_VDDR(i)*(Vm(i)-E_ca);
  I_ca_Ltype(i)=G_Ca_Ltype*d_Ltype(i)*f_Ltype(i)*f_ca_Ltype(i)*(Vm(i)-E_ca);
  I_NSCC(i)=G_NSCC*P_openNSCC(i)*(Vm(i)-E_NSCC);
  I_kv11(i)=G_Kv11*d_kv11(i)*f_kv11(i)*(Vm(i)-E_k);
  I_ERG(i)=G_ERG*d_ERG(i)*(Vm(i)-E_k);
  I_BK(i)=G_bk*p_open_BK*(Vm(i)-E_k);
  I_na(i)=G_Na*d_Na(i)*f_Na(i)*(Vm(i)-E_Na);
  I_Kb(i)=G_Kb*(Vm(i)-E_k);
  I_CaCl(i)=(G_Cacl*act_Cacl(i)*(Vm(i)-E_cl))+I_Ano1(i);
  I_pmca(i)=(2*F*1.0e12*Vcyto)*J_PMCA;

  I_ion = I_kv11(i)+I_NSCC(i)+I_ca_VDDR(i)+I_ca_Ltype(i)+I_na(i)+I_ERG(i)+I_CaCl(i)+I_Kb(i)+I_BK(i)+(2*F*1.0e12*Vcyto)*J_PMCA;

  % Initialise the dx vector to be zero
  

  d_Ca_cyto = fc*((((-I_ca_VDDR(i)-I_ca_Ltype(i))/(2*F*1.0e12*Vcyto)))-J_PMCA+J_leak);   % millimolar
  d_Vm = -I_ion/Cm;             % mV+
%   d_Vm = 0;
  
  % M&K 1998a Equation (34) : F&K 2001 removed need for gamma_1
  d_ADP_i = ((-J_ANT*(Vmito/Vcyto)) + (J_hyd - J_pGly));      % millimolar
  
  %-Update ADPm & ATPm-----------------------------------------------------
  % M&K 1998a Equation (25)
  d_ADP_m = (J_ANT - J_pTCA - J_pF1);                          % millimolar
  d_Ca_ER = fe*(J_SERCA - J_ERout);                            % millimolar
  d_Ca_PU = fc*( ((J_NaCa-J_uni)*(Vmito/V_PU))+ ((J_ERout-J_SERCA)*(Ver/V_PU))- J_leak*(Vcyto/V_PU)); % millimolar
  d_Ca_m = fm*(J_uni - J_NaCa);                                % millimolar
  
  % M&K 1998a Equation (22)
  d_NADH_m = (J_red - J_o);                                    % millimolar
  
  % M&K 1998a Equation (28)------------------------------------------------
  % Note on units. 10 power 6 instead of 10 power 12 because Cmito is given by
  % F&K in millifarad and not nanofarad.
  d_deltaPsi = - (1.0/Cmito)*(-J_Hres + J_HF1 + J_ANT + J_Hleak + (2*J_uni))*Vmito*F*1.0e6;  % mV
  
  %-Update h  F&K 2001 Equation p157---------------------------------------
  d_h = ((d_INH-((Ca_PU(i)+d_INH)*h(i)))/tauh);                 % dimensionless
  d_act_Cacl = (act_Cacl_inf-act_Cacl(i))/tau_act_Cacl;         % dimensionless
  d_d_ERG = (d_ERG_inf-d_ERG(i))/tau_ERG;                       % dimensionless
  d_d_Ltype = (d_inf_Ltype-d_Ltype(i))/tau_d_Ltype;             % dimensionless
  d_P_openNSCC = (P_openNSCC_inf-P_openNSCC(i))/tau_NSCC;       % dimensionless
  d_d_Na = (d_Na_inf-d_Na(i))/tau_d_Na;                         % dimensionless
  d_d_VDDR = (d_inf_VDDR-d_VDDR(i))/tau_d_VDDR;                 % dimensionless
  d_d_kv11 = (d_kv11_inf-d_kv11(i))/tau_d_kv11;                 % dimensionless
  d_f_Ltype = (f_inf_Ltype-f_Ltype(i))/tau_f_Ltype;             % dimensionless
  d_f_Na = (f_Na_inf-f_Na(i))/tau_f_Na;                         % dimensionless
  d_f_VDDR = (f_inf_VDDR-f_VDDR(i))/tau_f_VDDR;                 % dimensionless
  d_f_ca_Ltype = (f_Ca_Ltype_INF-f_ca_Ltype(i))/tau_f_Ca_Ltype; % dimensionless
  d_f_kv11 = (f_kv11_inf-f_kv11(i))/tau_f_kv11;                 % dimensionless
  
  Ca_cyto(i+1) = Ca_cyto(i)+(dt*d_Ca_cyto);             % millimolar
  Vm(i+1) = Vm(i)+(dt*d_Vm);                            % mV
  ADP_i(i+1) = ADP_i(i)+((dt*d_ADP_i));                 % millimolar
  ADP_m(i+1) = ADP_m(i)+((dt*d_ADP_m));                 % millimolar
  Ca_ER(i+1) = Ca_ER(i)+(dt*d_Ca_ER);                   % millimolar
  Ca_PU(i+1) = Ca_PU(i)+(dt*d_Ca_PU);                   % millimolar
  Ca_m(i+1) = Ca_m(i)+(dt*d_Ca_m);                      % millimolar
  NADH_m(i+1) = NADH_m(i)+(dt*d_NADH_m);                % millimolar
  %deltaPsi(i+1) = deltaPsi(i)+(dt*d_deltaPsi);          % mV
   deltaPsi(i+1) = deltaPsi(i)+(dt*0);          % mV

  h(i+1) = h(i)+(dt*d_h);                               % dimensionless
  act_Cacl(i+1) = act_Cacl(i)+(dt*d_act_Cacl);          % dimensionless
  d_ERG(i+1) = d_ERG(i)+(dt*d_d_ERG);                   % dimensionless
  d_Ltype(i+1) = d_Ltype(i)+(dt*d_d_Ltype);             % dimensionless
  P_openNSCC(i+1) = P_openNSCC(i)+(dt*d_P_openNSCC);    % dimensionless
  d_Na(i+1) = d_Na(i)+(dt*d_d_Na);                      % dimensionless
  d_VDDR(i+1) = d_VDDR(i)+(dt*d_d_VDDR);                % dimensionless
  d_kv11(i+1) = d_kv11(i)+(dt*d_d_kv11);                % dimensionless
  f_Ltype(i+1) = f_Ltype(i)+(dt*d_f_Ltype);             % dimensionless
  f_Na(i+1) = f_Na(i)+(dt*d_f_Na);                      % dimensionless
  f_VDDR(i+1) = f_VDDR(i)+(dt*d_f_VDDR);                % dimensionless
  f_ca_Ltype(i+1) = f_ca_Ltype(i)+(dt*d_f_ca_Ltype);    % dimensionless
  f_kv11(i+1) = f_kv11(i)+(dt*d_f_kv11);                % dimensionless
  
  end;

figure(1);
plot((1:Ts/dt),Vm);
% hold on;
% plot((1:Ts/dt),Vm,'r');
% xlabel('Time');
% ylabel('Membrane Voltage (mV)');
%plot((1:Ts/dt),deltaPsi);

% plotting ionic currents (part 1)
% figure(5);
% subplot(4,1,1);
% plot((1:Ts/dt),I_ca_VDDR);
% xlabel('Time');
% ylabel('I_ca_VDDR (pA)');
% 
% subplot(4,1,2);
% plot((1:Ts/dt),I_ca_Ltype);
% xlabel('Time');
% ylabel('I_ca_Ltype (pA)');
% 
% subplot(4,1,3);
% plot((1:Ts/dt),I_NSCC);
% xlabel('Time');
% ylabel('I_NSCC (pA)');
% 
% subplot(4,1,4);
% plot((1:Ts/dt),I_kv11);
% xlabel('Time');
% ylabel('I_kv11 (pA)');

% plotting ionic currents (part 2)
% figure(6);
% subplot(5,1,1);
% plot((1:Ts/dt),I_ERG);
% xlabel('Time');
% ylabel('I_ERG (pA)');
% 
% subplot(5,1,2);
% plot((1:Ts/dt),I_BK);
% xlabel('Time');
% ylabel('I_BK (pA)');
% 
% subplot(5,1,3);
% plot((1:Ts/dt),I_na);
% xlabel('Time');
% ylabel('I_na (pA)');
% 
% subplot(5,1,4);
% plot((1:Ts/dt),I_Kb);
% xlabel('Time');
% ylabel('I_Kb (pA)');
% 
% subplot(5,1,5);
% plot((1:Ts/dt), I_CaCl);
% xlabel('Time');
% ylabel(' I_CaCl (pA)');

% plotting calcium concentrations
% figure(7);
% subplot(4,1,1);
% plot((1:Ts/dt),Ca_cyto);
% xlabel('Time');
% ylabel('Ca_cyto(mM)');
% 
% subplot(4,1,2);
% plot((1:Ts/dt),Ca_ER);
% xlabel('Time');
% ylabel('Ca_ER (mM)');
% 
% subplot(4,1,3);
% plot((1:Ts/dt),I_na);
% xlabel('Time');
% ylabel('Ca_PU (mM)');
% 
% subplot(4,1,4);
% plot((1:Ts/dt),Ca_PU);
% xlabel('Time');
% ylabel('Ca_mitochondria (mM)');

% plotting ANO1 current
% figure(8);
% subplot(4,1,1);
% plot((1:Ts/dt),I_SOC);
% xlabel('Time');
% ylabel('I_SOC(pA)');
% 
% subplot(4,1,2);
% plot((1:Ts/dt),I_Ano1);
% xlabel('Time');
% ylabel('I_Ano1(pA)');
% 
% subplot(4,1,3);
% plot((1:Ts/dt),tau_Ano1);
% xlabel('Time');
% ylabel('tau_Ano1 (S)');
% 
% subplot(4,1,4);
% plot((1:Ts/dt),O_Ano1);
% xlabel('Time');
% ylabel('O_Ano1');

% figure(1);
% plot((1:Ts/dt),Vm,'r');
% % plot((1:Ts/dt),Ca_PU);
% xlabel('Time (s)');
% ylabel('Membrane Voltage (mV)');
% a = I_CaCl - I_Ano1;        %this is the icl_cacl current
% figure(2);
%I_kv11(i)+0*I_NSCC(i)+0*I_ca_VDDR(i)+0*I_ca_Ltype(i)+0*I_na(i)+I_ERG(i)+0*I_CaCl(i)+0*I_Kb(i)+0*I_BK(i)
% plot((1:Ts/dt),I_Ano1);
% title('icl_cacl');
