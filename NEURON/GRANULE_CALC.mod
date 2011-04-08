TITLE Cerebellum Granule Cell Model

COMMENT
        Calcium first order kinetics
   
	Author: A. Fontana
	Last revised: 12.12.98
ENDCOMMENT

NEURON {
        SUFFIX GRANULE_CALC
        USEION ca READ ica, cao WRITE cai
	RANGE Q10_diff,beta_Q10
        RANGE d, beta, cai0, ic
}

UNITS {
        (mV)    = (millivolt)
        (mA)    = (milliamp)
	(um)    = (micron)
	(molar) = (1/liter)
        (mM)    = (millimolar)
   	F      = (faraday) (coulomb)
}

PARAMETER {
        ica             (mA/cm2)
        ic             (mA/cm2)
        d = .2          (um)
        cao = 2.        (mM)         
        cai0 = 1e-4     (mM)
	Q10_diff = 3
        beta = 1.5        (/ms)
	celsius (degC)
}

ASSIGNED {
	beta_Q10 (mho/cm2)
}

STATE {
	cai (mM)
}

INITIAL {
	beta_Q10 = beta*(Q10_diff^((celsius-30)/10))
        cai = cai0 
}

BREAKPOINT {
    SOLVE conc METHOD derivimplicit
    ic = beta*(cai-cai0)
}

DERIVATIVE conc {    
	cai' = -ica/(2*F*d)*(1e4) - beta_Q10*(cai-cai0)
}

