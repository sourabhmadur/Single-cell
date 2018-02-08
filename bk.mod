TITLE bk channel by madur

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
}

NEURON {
    SUFFIX bk
    USEION k READ ek WRITE ik
	USEION ca READ cai
    RANGE G_bk,ik 
}

PARAMETER {
    G_bk = 23 (mho/cm2)
	kbk = -17
	hBK=2
	Caset=0.001 
	

}

ASSIGNED { 
    v (mV)
    ek (mV)
    ik (mA/cm2)
	p_open_BK
	cai (mM)
	
   }

BREAKPOINT {   
	p_open_BK=1.0/(1.0+exp((v/kbk)-hBK*log(cai/Caset)))
	ik=G_bk*p_open_BK*(v-ek)
	}


UNITSOFF


UNITSON