TITLE Kb channel by madur

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
}

NEURON {
    SUFFIX Kb
    USEION k READ ek WRITE ik
    RANGE G_Kb,ik
}

PARAMETER {
    G_Kb  = 0.15 (mho/cm2)
	
	
}

ASSIGNED { 
    v (mV)
    ek (mV)
    ik (mA/cm2)
    
    
}

BREAKPOINT {
    
    ik  = G_Kb*(v-ek)
}




	
UNITSOFF

UNITSON
