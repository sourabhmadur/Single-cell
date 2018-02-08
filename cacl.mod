TITLE cacl channel by madur

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
}

NEURON {
    SUFFIX cacl
	USEION cl READ ecl WRITE icl VALENCE -1 
    USEION ca READ cai VALENCE 2
    RANGE G_Cacl ,tau_act_Cacl,icl
}

PARAMETER {
     G_Cacl = 10.1					 (mho/cm2)
	 tau_act_Cacl=30         	:%*1.0e-3;
	 KClCa=140e-6
	 hClca=3
	
}

ASSIGNED { 
    v (mV)
    ecl (mV)
    icl (mA/cm2)
	cai	(mM)
	act_Cacl_inf
}

STATE {
    act_Cacl  
}

BREAKPOINT {  
	 				
	SOLVE states METHOD cnexp
	icl   =  G_Cacl*act_Cacl*(v-ecl)	
	}

INITIAL {
	
	act_Cacl = 0.0
	icl   =   G_Cacl*act_Cacl*(v-ecl)
				
			  }

DERIVATIVE states {  
    settables(v,cai)      
    act_Cacl' = (act_Cacl_inf-act_Cacl)/tau_act_Cacl         :% dimensionless
	}
	
UNITSOFF

PROCEDURE settables(v (mV),cai(mM)) {
	act_Cacl_inf=1.0/(1.0+((KClCa/cai)^hClca))
	
}


UNITSON