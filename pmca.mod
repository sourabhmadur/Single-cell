TITLE pmca channel by madur

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
}

NEURON {
    SUFFIX pmca
    USEION ca READ cai WRITE ica
    RANGE J_max_PMCA,ica
	
}

PARAMETER {
    J_max_PMCA= 0.088464 (millimolar/sec)
	k_PMCA = 0.000298	(millimolar)
	F = 96.4846 	 	(microcoulomb/nanomole)  		
	P_cyto=0.7
	Vol = 1.0e-12		(litre)
	fc = 0.01   
	
	
	}

ASSIGNED { 
    v (mV)
    eca (mV)
	cai	(mM)	
    ica (mA/cm2)
    }

BREAKPOINT {
    ica = (2*F*1.0e12*Vol*P_cyto)*J_max_PMCA*(1.0/(1.0+(k_PMCA/cai)))
	
	}



	
UNITSOFF


UNITSON