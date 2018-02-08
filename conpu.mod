TITLE conpu channel by madur

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
}

NEURON {
    SUFFIX conpu
    USEION ca READ cai 
	USEION caer	READ caeri WRITE caeri VALENCE 2
	USEION capu	READ capui WRITE capui VALENCE 2
	USEION cam	READ cami WRITE cami	VALENCE 2
	RANGE J_max_leak
	RANGE Jmax_serca
	RANGE Jmax_IP3
	RANGE Jmax_NaCa 
	RANGE J_ERleak
	RANGE Jmax_uni
}

PARAMETER {
    
		
	P_cyto=0.7
	Vol = 1.0e-12		(litre)
	fc = 0.01   
	J_max_leak = 0.01	(mM/s)
	
	:ER Ca Mechanism parameters
	P_PU = 0.001 
	Per = 0.10
	fe=0.01	
	Jmax_serca = 1.8333	(mM/s)
	k_serca = 0.00042  (mM)	
	Jmax_IP3 = 50000	(1/s)
	J_ERleak = 1.666667	(1/s)
	IP3 =0.0006       (millimolar) : baseline IP3 concentration
	d_IP3 = 0.00025   (millimolar)
	
	d_ACT = 0.001        (millimolar)
	d_INH = 0.0014       (millimolar)
	tauh = 4.0           (s)
	
	:Mitochondria mechanism parameters
	F = 96.4846         (microcoulomb/nanomole)
	R = 8.3144  
	:T = 279.6				(K)
	deltaPsi_star = 91   (mV)
	celsius
	
	K_Na = 9.4           (millimolar)
	K_Ca = 0.003         (millimolar)
	Na_i = 30			  (millimolar) 	
	b = 0.5              
	n = 2          
	Jmax_NaCa = 0.05   (mM/s)
	deltaPsi = 164.000044     (mV)
	fm = 0.0003
	Pmito = 0.12871
	:uniporter 
	Jmax_uni = 5000     (mM/s)
	naa = 2.8
	K_act = 0.00038		(mM)
	conc = 0.001        (mM) : correcting units of Juni
	K_trans = 0.006      (mM)
	L = 50        
		
}

ASSIGNED { 
	cai (mM)
	FoRT 	(1/mV):F/(R*T)       (mV)	
	
     }

STATE {
	capui (millimolar)	
	caeri (millimolar)
	cami (millimolar)	
	h
}

BREAKPOINT {
    SOLVE states METHOD cnexp
	
	}

INITIAL {
		settables()
		capui =  0.0000902
		caeri = 0.007299
		cami = 0.000136	
		h = 0.9397
		
			  }
			  
DERIVATIVE states {  

capui' = fc*(((J_NaCa(cami)-J_uni(cami))*((Vol*Pmito)/(Vol*P_PU)))+((J_ERout(capui,caeri,h)-J_SERCA(capui))*((Vol*Per)/(Vol*P_PU)))-J_leak(capui,cai)*((Vol*P_cyto)/(Vol*P_PU)))		
caeri' = fe*( J_SERCA(capui) - J_ERout(capui,caeri,h) )
cami' = fm*(J_uni(cami)- J_NaCa(cami))  

h' = ((d_INH-((capui+d_INH)*h))/tauh)	
		}
	
FUNCTION J_leak(capui(mM),cai(mM)){
	J_leak=J_max_leak*(capui - cai)	
}

FUNCTION J_SERCA(capui(mM)){
	J_SERCA= Jmax_serca*((capui^2)/((k_serca^2)+(capui^2)))
}

FUNCTION J_ERout(capui(mM),caeri(mM),h){
	J_ERout = (Jmax_IP3*((IP3/(IP3+d_IP3))^3)*((capui/(capui+d_ACT))^3)*(h^3)+J_ERleak)*(caeri-capui)
}

FUNCTION J_NaCa(cami(mM)){
	J_NaCa = Jmax_NaCa*(exp(b*FoRT*(deltaPsi-deltaPsi_star))/(((1+((K_Na/Na_i)^n)))*(1+(K_Ca/cami))))

}

FUNCTION MWC(capui(mM)){
	MWC = conc*((capui/K_trans)*((1+capui/K_trans)^3))/(((1+capui/K_trans)^4)+(L/((1+capui/K_act)^naa)))
}

FUNCTION J_uni(cami(mM)){
	 J_uni = Jmax_uni*((MWC(capui)-(cami*exp(-2*FoRT*(deltaPsi-deltaPsi_star))))*(2*FoRT*(deltaPsi-deltaPsi_star))/(1-exp(-2*FoRT*(deltaPsi-deltaPsi_star))))
}

PROCEDURE settables(){
	FoRT = F/(R*(celsius+273))	:F/(R*T)       (mV)	
}

	