TITLE Cerebellum Granule Cell Model

COMMENT
        KCa channel
   
	Author: E.D'Angelo, T.Nieus, A. Fontana
	Last revised: 8.5.2000
ENDCOMMENT
 
NEURON { 
	SUFFIX GRANULE_KCA
	USEION k READ ek WRITE ik
	USEION ca READ cai
	RANGE Q10_diff,Q10_channel,gbar_Q10
	RANGE gbar, ic, g, alpha_c, beta_c
	RANGE Aalpha_c, Balpha_c, Kalpha_c
	RANGE Abeta_c, Bbeta_c, Kbeta_c 
	RANGE c_inf, tau_c 
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
	(molar) = (1/liter)
	(mM) = (millimolar)
} 
 
PARAMETER { 
	Aalpha_c = 2.5 (/ms)
	Balpha_c = 1.5e-3 (mM)

	Kalpha_c =  -11.765 (mV)

	Abeta_c = 1.5 (/ms)
	Bbeta_c = 0.15e-3 (mM)

	Kbeta_c = -11.765 (mV)

	v (mV) 
	cai (mM)
	Q10_diff	= 1.5
	Q10_channel	= 3
	gbar= 0.0045 (mho/cm2) 
	ek = -84.69 (mV) 
	celsius (degC)
} 

STATE { 
	c 
} 

ASSIGNED { 
	ic (mA/cm2) 
	ik (mA/cm2)

	c_inf 
	tau_c (ms) 
	g (mho/cm2) 
	alpha_c (/ms) 
	beta_c (/ms) 
	gbar_Q10 (mho/cm2)
} 
 
INITIAL { 
	gbar_Q10 = gbar*(Q10_diff^((celsius-30)/10))
	rate(v) 
	c = c_inf 
} 
 
BREAKPOINT { 
	SOLVE states METHOD derivimplicit 
	g = gbar_Q10*c 
	ik = g*(v - ek) 
	ic = ik
	alpha_c = alp_c(v) 
	beta_c = bet_c(v) 
} 
 
DERIVATIVE states { 
	rate(v) 
	c' =(c_inf - c)/tau_c 
} 
 
FUNCTION alp_c(v(mV))(/ms) { LOCAL Q10
	Q10 = Q10_channel^((celsius-30(degC))/10(degC))
	alp_c = Q10*Aalpha_c/(1+(Balpha_c*exp(v/Kalpha_c)/cai)) 
} 
 
FUNCTION bet_c(v(mV))(/ms) { LOCAL Q10
	Q10 = Q10_channel^((celsius-30(degC))/10(degC))
	bet_c = Q10*Abeta_c/(1+cai/(Bbeta_c*exp(v/Kbeta_c))) 
} 
 
PROCEDURE rate(v (mV)) {LOCAL a_c, b_c 
	TABLE c_inf, tau_c 
	DEPEND Aalpha_c, Balpha_c, Kalpha_c, 
	       Abeta_c, Bbeta_c, Kbeta_c, celsius FROM -100 TO 30 WITH 13000 
	a_c = alp_c(v)  
	b_c = bet_c(v) 
	tau_c = 1/(a_c + b_c) 
	c_inf = a_c/(a_c + b_c) 
} 

