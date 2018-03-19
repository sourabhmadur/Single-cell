 

T=310
T_exp=297           # Kelvin, temperature of experiments      
F = 96.4846         # microcoulomb_per_nanomole
R = 8.3144


 # Temperature Corrections
Q10Ca=2.1
Q10K=1.5
Q10Na=2.45
Tcorrection_Ca=(Q10Ca**((T-T_exp)/10))
Tcorrection_K=(Q10K**((T-T_exp)/10))
Tcorrection_Na=(Q10Na**((T-T_exp)/10))
T_correction_BK=1.1*(T-T_exp)

t_corr = 1000

tau_d_VDDR=(6)*1.0e-3*Tcorrection_Ca*t_corr
tau_f_VDDR=(40)*1.0e-3*Tcorrection_Ca*t_corr
tau_d_Ltype=1*1.0e-3*Tcorrection_Ca*t_corr
tau_f_Ltype=86*1.0e-3*Tcorrection_Ca*t_corr
tau_f_Ca_Ltype=2.0*1.0e-3*Tcorrection_Ca*t_corr
tau_d_kv11=5.0*1.0e-3*Tcorrection_K*t_corr
tau_f_kv11=5*1.0e-3*Tcorrection_K*t_corr
tau_ERG=3*1.0e-3*Tcorrection_K*t_corr
tau_d_Na=3*1.0e-3*Tcorrection_Na*t_corr
tau_f_Na=1.6*1.0e-3*Tcorrection_Na*t_corr
tau_NSCC=350*1.0e-3*t_corr
tau_act_Cacl=30*1.0e-3*t_corr

