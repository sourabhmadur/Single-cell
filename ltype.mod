TITLE ltype channel by madur

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
}

NEURON {
    SUFFIX ltype
    USEION ca READ eca,cai WRITE ica
    RANGE G_Ca_Ltype,tau_f_Ltype,tau_d_Ltype,tau_f_Ca_Ltype,ica
	
}

PARAMETER {
    G_Ca_Ltype= 2 (mho/cm2)
	tau_d_Ltype=1              (ms) :*1.0e-3*Tcorrection_Ca;
	tau_f_Ltype=86              (ms):*1.0e-3*Tcorrection_Ca;
	tau_f_Ca_Ltype=2.0       (ms)   :*1.0e-3*Tcorrection_Ca;
	
}

ASSIGNED { 
    v (mV)
    eca (mV)
	cai	(mM)	
    ica (mA/cm2)
	
    d_inf_Ltype
	f_inf_Ltype
	f_Ca_Ltype_INF
	}

STATE {
	d_Ltype f_Ltype f_ca_Ltype 
}

BREAKPOINT {
    SOLVE states METHOD cnexp
	ica=G_Ca_Ltype*d_Ltype*f_Ltype*f_ca_Ltype*(v-eca)
	}

INITIAL {
	f_ca_Ltype = 1.0 
	f_Ltype = 1.0
	d_Ltype = 0.0
	 
			  }

DERIVATIVE states {  
    settables(v,cai)      
	f_ca_Ltype' = (f_Ca_Ltype_INF-f_ca_Ltype)/tau_f_Ca_Ltype : dimensionless
    f_Ltype' = (f_inf_Ltype-f_Ltype)/tau_f_Ltype             : dimensionless  
    d_Ltype' = (d_inf_Ltype-d_Ltype)/tau_d_Ltype             : dimensionless
	
	}
	
UNITSOFF

PROCEDURE settables(v (mV),cai(mM)) {
   
	
	d_inf_Ltype=1.0/(1.0+exp(-(17.0+v)/4.3))
	f_inf_Ltype=1.0/(1.0+exp((43.0+v)/8.9))
	f_Ca_Ltype_INF=1.0-(1.0/(1.0+exp(-((cai-0.0001)-0.000214)/0.0000131)))
	    
}

UNITSON