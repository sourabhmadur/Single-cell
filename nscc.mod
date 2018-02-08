TITLE nscc channel by madur

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
}

NEURON {
    SUFFIX nscc
    NONSPECIFIC_CURRENT i
    RANGE G_NSCC,i, E_NSCC ,tau_NSCC ,i
	USEION capu	READ capui
}

PARAMETER {
	G_NSCC = 12.15	(mho/cm2)
    tau_NSCC=350		(ms):e-3
	E_NSCC = 0
	Ca_PU = 0.0000902
	K_NSCC=0.0000745
	hNSCC=-85
	}

ASSIGNED { 
    v (mV)
    i (mA/cm2)
    P_openNSCC_inf
	capui (mM)
	
}

STATE {
    P_openNSCC
}

BREAKPOINT {
    SOLVE states METHOD cnexp
	i=G_NSCC*P_openNSCC*(v-E_NSCC)
	}

INITIAL {
	P_openNSCC=0
			  }

DERIVATIVE states {  
    settables(v)      
	P_openNSCC' = (P_openNSCC_inf-P_openNSCC)/tau_NSCC
	}
	
UNITSOFF

PROCEDURE settables(v (mV)) {
   
	P_openNSCC_inf=1.0/(1.0+((K_NSCC/capui)^hNSCC))
	
    
}

UNITSON