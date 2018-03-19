TITLE ltype channel by madur

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
}

NEURON {
    SUFFIX concyto
    USEION ca READ ica WRITE cai 
	USEION capu READ capui VALENCE 2
	RANGE J_max_leak,Jmax_serca,Vol,a

	}

PARAMETER {
    
           
	F = 96.4846 	 	(microcoulomb/nanomole)  		
	P_cyto=0.7
	Vol = 1.0e-12		(um^3)
	fc = 0.01  
	J_max_leak = 0.00001 	(1/ms)
	a=1000			(um^2)
	corrfactor = 10		:page 13 of the notebook
	
		
	
}

ASSIGNED { 
    v (mV)
    eca (mV)
    ica (mA/cm2)
	capui (mM)
	
	
	}

STATE {
	cai (millimolar)	
}

BREAKPOINT {
    SOLVE states METHOD derivimplicit
	
	}

INITIAL {
		cai =  0.00000993087 
		
			  }

DERIVATIVE states {  
       
	cai' = fc*((((-corrfactor*ica*a)/(2*F*Vol*P_cyto))) + J_leak(capui,cai))
		
	}
	
FUNCTION J_leak(capui(mM),cai(mM)){

	J_leak=J_max_leak*(capui-cai)	

}
	
