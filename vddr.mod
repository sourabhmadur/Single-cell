TITLE vddr channel by madur

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
}

NEURON {
    SUFFIX vddr
    USEION ca READ eca WRITE ica
    RANGE G_Ca_VDDR,tau_d_VDDR,tau_f_VDDR,ica
}

PARAMETER {
    G_Ca_VDDR = 3(mho/cm2)
	tau_d_VDDR=6                         :*1.0e-3*Tcorrection_Ca;
	tau_f_VDDR=40                        :*1.0e-3*Tcorrection_Ca;

}

ASSIGNED { 
    v (mV)
    eca (mV)
    ica (mA/cm2)
    d_inf_VDDR
	f_inf_VDDR
}

STATE {
    d_VDDR f_VDDR 
}

BREAKPOINT {
    SOLVE states METHOD cnexp
	ica=G_Ca_VDDR*d_VDDR*f_VDDR*(v-eca)  

	}

INITIAL {
	d_VDDR= 0
	f_VDDR=1
			  }

DERIVATIVE states {  
    settables(v)      
	d_VDDR' = (d_inf_VDDR-d_VDDR)/tau_d_VDDR
	f_VDDR' = (f_inf_VDDR-f_VDDR)/tau_f_VDDR
	}
	
UNITSOFF

PROCEDURE settables(v (mV)) {
		
	d_inf_VDDR=1.0/(1.0+exp(-(26.0+v)/6.0))
	f_inf_VDDR=1.0/(1.0+exp((66.0+v)/6.0))

}

UNITSON