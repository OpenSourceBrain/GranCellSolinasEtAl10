TITLE Cerebellum Granule Cell Model

COMMENT
        Gaba A leakage
   
	Author: A. Fontana
	Last revised: 18.2.99
ENDCOMMENT
 
NEURON { 
	SUFFIX GRANULE_LKG2 
	NONSPECIFIC_CURRENT il
	RANGE Q10_diff,Q10_channel,gbar_Q10
	RANGE egaba, g , ic, gbar
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
} 
 
PARAMETER { 
	v (mV) 
	gbar = 6e-5 (mho/cm2) : Increased of 200% for Jorntell
	egaba = -65 (mV)
	Q10_diff	= 1.5
	celsius (degC)
} 

ASSIGNED { 
	il (mA/cm2) 
	ic (mA/cm2) 
	g (mho/cm2)
	gbar_Q10 (mho/cm2)
}

BREAKPOINT { 
    gbar_Q10 = gbar*(Q10_diff^((celsius-30)/10))
    g = gbar_Q10
    il = g*(v - egaba) 
    ic =il
} 
