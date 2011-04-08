TITLE Cerebellum Granule Cell Model

COMMENT
        pNa channel
   
	Author: E.D'Angelo, T.Nieus, A. Fontana 
	Last revised: 8.5.2000
ENDCOMMENT
 
NEURON { 
	SUFFIX GRANULE_PNA 
	USEION na READ ena WRITE ina 
	RANGE Q10_diff,Q10_channel,gbar_Q10
	RANGE gbar, ina, ic, g, alpha_m, beta_m
	RANGE Aalpha_m, Kalpha_m, V0alpha_m
	RANGE Abeta_m, Kbeta_m, V0beta_m
	RANGE V0_minf, B_minf
	RANGE m_inf, tau_m
	:GLOBAL i
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
} 
 
PARAMETER { 
	Aalpha_m 	= -0.91 (/mV-ms) 	: -0.091
	Kalpha_m 	= -5 (mV)			: CONST
	V0alpha_m 	= -42 (mV)	 		: -42,-40
	Abeta_m 	= 0.62 (/mV-ms) 	: 0.062
	Kbeta_m 	= 5 (mV)			: CONST
	V0beta_m 	= -42 (mV)			: -42,-40 
	V0_minf 	= -42 (mV)			: -42,-40 
	B_minf 		= 5 (mV)			: CONST
	v (mV) 
	gbar		= 2e-5 (mho/cm2)
	Q10_diff	= 1.5
	Q10_channel	= 3
	ena 		= 87.39 (mV) 
	celsius (degC)
} 

STATE { 
	m 
} 

ASSIGNED { 
	ina (mA/cm2) 
	ic (mA/cm2) 
	m_inf 
	tau_m (ms) 
	g (mho/cm2) 
	alpha_m (/ms)
	beta_m (/ms)
	gbar_Q10 (mho/cm2)
} 
 
INITIAL { 
	gbar_Q10 = gbar*(Q10_diff^((celsius-30)/10))
	rate(v) 
	m = m_inf 
} 
 
BREAKPOINT { 
	SOLVE states METHOD derivimplicit 
	g = gbar_Q10*m 
	ina = g*(v - ena) 
	ic = ina
	alpha_m = alp_m(v)
	beta_m = bet_m(v)
} 
 
DERIVATIVE states { 
	rate(v) 
	m' =(m_inf - m)/tau_m 
} 

FUNCTION alp_m(v(mV))(/ms) { LOCAL Q10
	Q10 = Q10_channel^((celsius-30(degC))/10(degC))
	alp_m = Q10 * Aalpha_m*linoid(v-V0alpha_m, Kalpha_m) 
} 
 
FUNCTION bet_m(v(mV))(/ms) { LOCAL Q10
	Q10 = Q10_channel^((celsius-30(degC))/10(degC))
	bet_m = Q10 * Abeta_m*linoid(v-V0beta_m, Kbeta_m) 
} 
 
PROCEDURE rate(v (mV)) {LOCAL a_m, b_m 
	TABLE m_inf, tau_m 
	DEPEND Aalpha_m, Kalpha_m, V0alpha_m, 
	       Abeta_m, Kbeta_m, V0beta_m, celsius FROM -100 TO 30 WITH 13000
	a_m = alp_m(v)  
	b_m = bet_m(v) 
:	m_inf = a_m/(a_m + b_m) 
	m_inf = 1/(1+exp(-(v-V0_minf)/B_minf))
	tau_m = 5/(a_m + b_m) 
} 

FUNCTION linoid(x (mV),y (mV)) (mV) {
        if (fabs(x/y) < 1e-6) {
                linoid = y*(1 - x/y/2)
        }else{
                linoid = x/(exp(x/y) - 1)
        }
}


