TITLE NMDA leakage
COMMENT
	NMDA: descritto da Westbrook
	derived from the kinetic scheme from Nieus2006 and also folder "DeterministicApproximated"
	rates of the kinetic scheme from Rossi2002
	
ENDCOMMENT

NEURON {
	SUFFIX GRANULE_Nmda_leak
	NONSPECIFIC_CURRENT i
	RANGE Q10_diff
	RANGE g , ic, Erev, gmax, gmax_factor
	RANGE MgBlock,v0_block,k_block
}

UNITS {
    (nA) = (nanoamp)	
    (mV) = (millivolt)
    (umho) = (micromho)
    (mM) = (milli/liter)
    (pS) = (picosiemens)
    (nS) = (nanosiemens)
    (um) = (micrometer)
    PI	= (pi)		(1)
}

PARAMETER {
    gmax_factor = 1
    gmax	= 50e-12	(mho)
    surf        = 299.26e-8 (cm2)
    Erev	= -3.7  (mV)	: 0 (mV)
    Q10_diff	= 1.5
    v0_block = -20 (mV)	
    k_block  = 13 (mV)
    O = 10e-3 (s) : mean open time
    Or = 52 : spontaneous opening rate  (Hz)
    v		(mV)		: postsynaptic voltage
    celsius (degC)
}

ASSIGNED {
    i 		(mA/cm2)		: current = g*(v - Erev)
    ic 		(mA/cm2)		: current = g*(v - Erev)
    g 		(mho/cm2)		: actual conductance
    MgBlock
    gbar_Q10 (1)
}

INITIAL {
	rates(v)
	gbar_Q10 = Q10_diff^((celsius-30)/10)
}

BREAKPOINT {
	rates(v)
	g = gmax / surf * gbar_Q10 * O * Or * gmax_factor  * MgBlock
	i = g * (v - Erev) 
	ic = i
}

PROCEDURE rates(v(mV)) {
	: E' necessario includere DEPEND v0_block,k_block per aggiornare le tabelle!
	TABLE MgBlock DEPEND v0_block,k_block FROM -120 TO 30 WITH 150
	MgBlock = 1 / ( 1 + exp ( - ( v - v0_block ) / k_block ) )
}

