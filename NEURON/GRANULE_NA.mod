TITLE Cerebellum Granule Cell Model

COMMENT
        Na channel
	Gutfreund parametrization
   
	Author: E.D'Angelo, T.Nieus, A. Fontana
	Last revised: 8.5.2000
ENDCOMMENT
 
NEURON { 
	SUFFIX GRANULE_NA 
	USEION na READ ena WRITE ina 
	RANGE Q10_diff,Q10_channel,gbar_Q10, gbar
	RANGE gnabar, ic, g, alpha_m, beta_m, alpha_h, beta_h 
	RANGE Aalpha_m, Kalpha_m, V0alpha_m
	RANGE Abeta_m, Kbeta_m, V0beta_m

	RANGE Aalpha_h, Kalpha_h, V0alpha_h
	RANGE Abeta_h, Kbeta_h, V0beta_h

	RANGE m_inf, tau_m, h_inf, tau_h 
	RANGE alpha_m, beta_m, alpha_h, beta_h
	RANGE Q10_channel_alp_m, Q10_channel_bet_m,Q10_channel_alp_h,Q10_channel_bet_h
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
} 
 
PARAMETER { 

	Aalpha_m = -0.3 (/ms-mV)
	Kalpha_m = -10 (mV)
	V0alpha_m = -19 (mV)
	
	Abeta_m = 12 (/ms)
	Kbeta_m = -18.182 (mV)
	V0beta_m = -44 (mV)

	Aalpha_h  = 0.105 (/ms)
	Kalpha_h  = -3.333 (mV)
	V0alpha_h = -44 (mV)
 
	Abeta_h  = 0.75 (/ms)
	Kbeta_h  = -5 (mV)
	V0beta_h = -11 (mV)

	v (mV) 
	Q10_diff	= 1.5
	Q10_channel_alp_m	= 3
	Q10_channel_bet_m	= 3
	Q10_channel_alp_h	= 3
	Q10_channel_bet_h	= 3
	gbar	=  0.013 (mho/cm2)
	ena 	= 87.39 (mV) 
	celsius (degC)
} 

STATE { 
	m 
	h 
} 

ASSIGNED { 
	ina (mA/cm2) 
	ic (mA/cm2) 
	m_inf 
	h_inf 
	tau_m (ms) 
	tau_h (ms) 
	g (mho/cm2) 
	alpha_m (/ms)
	beta_m (/ms)
	alpha_h (/ms)
	beta_h (/ms)
	gbar_Q10 (mho/cm2)
} 
 
INITIAL { 
	gbar_Q10 = gbar*(Q10_diff^((celsius-30)/10))
	rate(v) 
	m = m_inf 
	h = h_inf 
} 
 
BREAKPOINT { 
	SOLVE states METHOD derivimplicit 
	g = gbar_Q10*m*m*m*h 
	ina = g*(v - ena)
	ic = ina
	alpha_m = alp_m(v)
	beta_m = bet_m(v) 
	alpha_h = alp_h(v)
	beta_h = bet_h(v) 
} 
	 
DERIVATIVE states { 
	rate(v) 
	m' =(m_inf - m)/tau_m 
	h' =(h_inf - h)/tau_h 
} 
 
FUNCTION alp_m(v(mV))(/ms) { LOCAL Q10
	Q10 = Q10_channel_alp_m^((celsius-20(degC))/10(degC)) 
	alp_m = Q10*Aalpha_m*linoid(v-V0alpha_m,Kalpha_m) 
} 
 
FUNCTION bet_m(v(mV))(/ms) { LOCAL Q10
	Q10 = Q10_channel_bet_m^((celsius-20(degC))/10(degC)) 
	bet_m = Q10*Abeta_m*exp((v-V0beta_m)/Kbeta_m) 
} 
 
FUNCTION alp_h(v(mV))(/ms) { LOCAL Q10
	Q10 = Q10_channel_alp_h^((celsius-20(degC))/10(degC)) 
	alp_h = Q10*Aalpha_h*exp((v-V0alpha_h)/Kalpha_h) 
} 
 
FUNCTION bet_h(v(mV))(/ms) { LOCAL Q10 
    Q10 = Q10_channel_bet_h^((celsius-20(degC))/10(degC)) 
    bet_h = Q10*Abeta_h/(1+exp((v-V0beta_h)/Kbeta_h))
} 
 
PROCEDURE rate(v (mV)) {LOCAL a_m, b_m, a_h, b_h 
	TABLE m_inf, tau_m, h_inf, tau_h 
	DEPEND Aalpha_m, Kalpha_m, V0alpha_m, 
	       Abeta_m, Kbeta_m, V0beta_m,
               Aalpha_h, Kalpha_h, V0alpha_h,
               Abeta_h, Kbeta_h, V0beta_h, celsius, Q10_channel_bet_h, Q10_channel_alp_h, Q10_channel_bet_m, Q10_channel_alp_m FROM -100 TO 30 WITH 13000 
	a_m = alp_m(v)  
	b_m = bet_m(v) 
	a_h = alp_h(v)  
	b_h = bet_h(v) 
	m_inf = a_m/(a_m + b_m) 
	tau_m = 1/(a_m + b_m) 
	h_inf = a_h/(a_h + b_h) 
	tau_h = 1/(a_h + b_h) 
	:if (tau_h<0.1 (ms)) {tau_h=0.1 (ms)} : riga aggiunta il 10 giugno 2003
} 

FUNCTION linoid(x (mV),y (mV)) (mV) {
        if (fabs(x/y) < 1e-6) {
                linoid = y*(1 - x/y/2)
        }else{
                linoid = x/(exp(x/y) - 1)
        }
}

