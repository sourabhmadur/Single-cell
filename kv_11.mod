TITLE kv11 channel by madur

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
}

NEURON {
    SUFFIX kv11
    USEION k READ ek WRITE ik
    RANGE G_Kv11, tau_d_kv11, tau_f_kv11,ik
}

PARAMETER {
    G_Kv11= 6.3 (mho/cm2)
	tau_d_kv11=  1       (ms)                           		:5.0*1.0e-3*Tcorrection_K;
	tau_f_kv11= 3          (ms)                                 	:5*1.0e-3*Tcorrection_K;

}

ASSIGNED { 
    v (mV)
    ek (mV)
    ik (mA/cm2)
    d_kv11_inf
	f_kv11_inf

}

STATE {
    d_kv11 f_kv11 	
}

BREAKPOINT {
    SOLVE states METHOD cnexp

	ik=G_Kv11*d_kv11*f_kv11*(v-ek)
	}

INITIAL {
	d_kv11= 0
	f_kv11=1

			  }

DERIVATIVE states {  
    settables(v)      
	d_kv11' = (d_kv11_inf-d_kv11)/tau_d_kv11
	f_kv11' = (f_kv11_inf-f_kv11)/tau_f_kv11

	}
	
UNITSOFF

PROCEDURE settables(v (mV)) {

	
	d_kv11_inf=1/(1+exp(-(25+v)/7.7))
	f_kv11_inf=0.5/(1+exp((44.8+v)/4.4))+0.5

	
    
}

UNITSON