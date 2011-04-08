TITLE Cerebellum Granule Cell Model

COMMENT
        Kir channel
   
	Author: E.D'Angelo, T.Nieus, A. Fontana
	Last revised: 8.10.2000
	Old values:
			gkbar = 0.0003 (mho/cm2) 
			
ENDCOMMENT
 
NEURON { 
	SUFFIX GRANULE_KIR
	USEION k READ ek WRITE ik 
	RANGE Q10_diff,Q10_channel,gbar_Q10
	RANGE gbar, ic, g, alpha_d, beta_d 
	RANGE Aalpha_d, Kalpha_d, V0alpha_d
	RANGE Abeta_d, Kbeta_d, V0beta_d
	RANGE d_inf, tau_d 
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
} 
 
PARAMETER { 
	Aalpha_d = 0.13289 (/ms)

	Kalpha_d = -24.3902 (mV)

	V0alpha_d = -83.94 (mV)
	Abeta_d = 0.16994 (/ms)

	Kbeta_d = 35.714 (mV)

	V0beta_d = -83.94 (mV)
	v (mV) 
	Q10_diff	= 1.5
	Q10_channel	= 3
	gbar = 0.00135 (mho/cm2) : increased by 150% for Jorntell
	ek = -84.69 (mV) 
	celsius (degC)
} 

STATE { 
	d 
} 

ASSIGNED { 
	ik (mA/cm2) 
	ic (mA/cm2) 
	d_inf 
	tau_d (ms) 
	g (mho/cm2) 
	alpha_d (/ms) 
	beta_d (/ms) 
	gbar_Q10 (mho/cm2)
} 
 
INITIAL { 
	gbar_Q10 = gbar*(Q10_diff^((celsius-30)/10))
	rate(v) 
	d = d_inf 
} 
 
BREAKPOINT { 
	SOLVE states METHOD derivimplicit
	g = gbar_Q10*d   : primo ordine!!!
	ik = g*(v - ek) 
	ic = ik
	alpha_d = alp_d(v) 
	beta_d = bet_d(v) 
} 
 
DERIVATIVE states { 
	rate(v) 
	d' =(d_inf - d)/tau_d 
} 
 
FUNCTION alp_d(v(mV))(/ms) { LOCAL Q10
	Q10 = Q10_channel^((celsius-20(degC))/10(degC))
	alp_d = Q10*Aalpha_d*exp((v-V0alpha_d)/Kalpha_d) 
} 
 
FUNCTION bet_d(v(mV))(/ms) { LOCAL Q10
	Q10 = Q10_channel^((celsius-20(degC))/10(degC))
	bet_d = Q10*Abeta_d*exp((v-V0beta_d)/Kbeta_d) 
} 
 
PROCEDURE rate(v (mV)) {LOCAL a_d, b_d 
	TABLE d_inf, tau_d  
	DEPEND Aalpha_d, Kalpha_d, V0alpha_d, 
	       Abeta_d, Kbeta_d, V0beta_d, celsius FROM -100 TO 30 WITH 13000 
	a_d = alp_d(v)  
	b_d = bet_d(v) 
	tau_d = 1/(a_d + b_d) 
	d_inf = a_d/(a_d + b_d) 
} 

