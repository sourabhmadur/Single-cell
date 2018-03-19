TITLE pmca channel by madur

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
}

NEURON {
    SUFFIX pmca
    USEION ca READ cai WRITE ica
    RANGE J_max_PMCA,ica,Vol,a
	
}

PARAMETER {
    J_max_PMCA= 0.088464 (millimolar/sec)
	k_PMCA = 0.000298	(millimolar)
	F = 96.4846 	 	(microcoulomb/nanomole)  		
	P_cyto=0.7
	Vol = 1.0e-12		(lit)		
	fc = 0.01   
	corrfactor = 1.0e4		:page 13 of the notebook
	a=1000			(um^2)
	
	}

ASSIGNED { 
    v (mV)
    eca (mV)
	cai	(mM)	
    ica (mA/cm2)
    }

BREAKPOINT {
    ica = ( (2*F*Vol*P_cyto)*J_max_PMCA*(1.0/(1.0+(k_PMCA/cai))) )/ ((corrfactor)*a)
	}



	
UNITSOFF


UNITSON