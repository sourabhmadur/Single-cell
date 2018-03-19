TITLE Na channel by madur

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
}

NEURON {
    SUFFIX Na
    USEION na READ ena WRITE ina
    RANGE G_Na,tau_d_Na,tau_f_Na,ina
}

PARAMETER {
    G_Na  = 20        (mho/cm2)
	tau_d_Na = 1       (ms)
	tau_f_Na =5        (ms)
}

ASSIGNED { 
    v (mV)
    ena (mV)
    ina (mA/cm2)
    d_Na_inf
	f_Na_inf
}

STATE {
    d_Na f_Na 
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ina  = G_Na*d_Na*f_Na*(v-ena)
}

INITIAL {
	d_Na=0
	f_Na=1
			  }

DERIVATIVE states {  
    settables(v)      
	d_Na' = (d_Na_inf-d_Na)/tau_d_Na
	f_Na' = (f_Na_inf-f_Na)/tau_f_Na
}
	
UNITSOFF

PROCEDURE settables(v (mV)) {
   

	d_Na_inf=1/(1+exp(-(47+v)/4.8))
	f_Na_inf=1/(1+exp((78+v)/(7)))
	
	
    
}

UNITSON