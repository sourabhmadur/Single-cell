TITLE ano1 channel by madur

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
}

NEURON {
    SUFFIX ano1
	USEION cl READ ecl WRITE icl VALENCE -1 
    USEION capu READ capui VALENCE 2
    RANGE g_Ano1,icl 
}

PARAMETER {
    g_Ano1 = 20 					 (mho/cm2)
	EC_50_i     =   1.39e-3             (mM)
	k_C         =   0.01248              (1/mV)
	V_h         =   -100                 (1/mV)
	k_V         =   0.0156               (1/mV)
}

ASSIGNED { 
    v 	(mV)
	
    ecl 	(mV)
    icl 	(mA/cm2)
	t1 t2 t3 tau_Ano1		
	EC_50
	O_Ano1_inf
	capui	(mM)
}

STATE {
    O_Ano1 
}

BREAKPOINT {  
	 				
	SOLVE states METHOD cnexp
	icl   =   g_Ano1*O_Ano1*(v-ecl)	
	}

INITIAL {
	settables(v,capui) 
	O_Ano1=1/((1+exp((V_h-v)*k_V))*(1+((EC_50/capui)^2)))
	icl   =   g_Ano1*O_Ano1*(v-ecl)	
				
			  }

DERIVATIVE states {  
    settables(v,capui)      
    O_Ano1'=((O_Ano1_inf-O_Ano1)/tau_Ano1)
	
	}
	
UNITSOFF

PROCEDURE settables(v (mV),capui(mM)) {
	t1          =   0.08163*exp(-0.57*capui)
    t2          =   0.07617*exp(-0.05374*capui)
	t3          =   70.3*exp(0.153*capui)
	tau_Ano1 =  ( t1 +(t2*exp(v/t3))  )*1.0e3		:this is the corrention from seconds to ms
	EC_50       =   EC_50_i*exp(-k_C*v)
	O_Ano1_inf  =   1/((1+exp((V_h-v)*k_V))*(1+((EC_50/capui)^2)))
}


UNITSON