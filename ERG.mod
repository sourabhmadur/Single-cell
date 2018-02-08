TITLE ERG channel by madur

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
}

NEURON {
    SUFFIX ERG
    USEION k READ ek WRITE ik
    RANGE G_ERG,tau_ERG,ik
}

PARAMETER {
    G_ERG  = 2.5 (mho/cm2)
	tau_ERG = 3	(ms)
	
}

ASSIGNED { 
    v (mV)
    ek (mV)
    ik (mA/cm2)
    d_ERG_inf
    
}

STATE {
    d_ERG
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ik  = G_ERG*d_ERG*(v-ek)
}

INITIAL {

	d_ERG = 0
			  }

DERIVATIVE states {  
    settables(v)      
	d_ERG' = (d_ERG_inf-d_ERG)/tau_ERG
    
}
	
UNITSOFF

PROCEDURE settables(v (mV)) {

    d_ERG_inf=0.8/(1+exp(-(20+v)/1.8))+0.2
	
    
}

UNITSON
