TITLE 

COMMENT
ENDCOMMENT

NEURON {
	POINT_PROCESS GRANULE_Ampa_det_Prel
	NONSPECIFIC_CURRENT i
	RANGE Q10_diff,Q10_channel
	RANGE R, g, ic
	RANGE Cdur, Erev 
	RANGE r1FIX, r2, r3,r4,r5,gmax,r1,r6,r6FIX,kB
	RANGE tau_1, tau_rec, tau_facil, U	 
	RANGE PRE,T,Tmax
		
	RANGE diffuse,Trelease,lamd
	RANGE M,Diff,Rd
	
	RANGE tspike
	RANGE nd, syntype, gmax_factor
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
	syntype
	gmax_factor = 1
	Q10_diff	= 1.5
	Q10_channel	= 2.4
	: Parametri Postsinaptici
	gmax		= 1200   (pS)	

	U 			= 0.42 (1) 	< 0, 1 >
	tau_rec 	= 35.1 (ms) 	< 1e-9, 1e9 > 	 
	tau_facil 	= 10.8 (ms) 	< 0, 1e9 > 	

	M			= 21.515 : numero di (kilo) molecole in una vescicola		
	Rd			= 1.03 (um)
	Diff		= 0.223 (um2/ms)
		 
	Cdur		= 0.3	(ms)		 
	r1FIX		= 5.4		(/ms/mM) 	 				
	r2			= 0.82	(/ms)			 
	r3			= 0		(/ms)		 
	r4			= 0		(/ms)		 
	r5			= 0.013	(/ms)			 
	r6FIX		= 1.12	(/ms/mM)		 
	Erev		= 0	(mV)
	kB			= 0.44	(mM)

	: Parametri Presinaptici
	tau_1 		= 3 (ms) 	< 1e-9, 1e9 >
	 

	
	u0 			= 0 (1) 	< 0, 1 >	: se u0=0 al primo colpo y=U
	Tmax		= 1  (mM)	
	
	: Diffusion			
	diffuse		= 1
	
	lamd		= 20 (nm)
	nd			= 1
	celsius (degC)
}


ASSIGNED {
	v		(mV)		: postsynaptic voltage
	i 		(nA)		: current = g*(v - Erev)
	ic 		(nA)		: current = g*(v - Erev)
	g 		(pS)		: conductance
	r1		(/ms)
	r6		(/ms)
	T		(mM)

	Trelease	(mM)
	tspike[5000]	(ms)
	x 
	tsyn		(ms)
	PRE[5000]
	
	Mres		(mM)	
	numpulses
	tzero
	gbar_Q10 (mho/cm2)
	Q10 (1)
}

STATE {	
	C
	O
	D
}

INITIAL {
	C=1
	O=0
	D=0
	T=0 (mM)
	numpulses=0
	Trelease=0 (mM)
	tspike[0]=1e12	(ms)
	 
	gbar_Q10 = Q10_diff^((celsius-30)/10)
	Q10 = Q10_channel^((celsius-30)/10)

	: fattore di conversione che comprende molecole -> mM
	: n molecole/(Na*V) -> M/(6.022e23*1dm^3)

	Mres = 1e3 * ( 1e3 * 1e15 / 6.022e23 * M ) : (M) to (mM) so 1e3, 1um^3=1dm^3*1e-15 so 1e15	
	numpulses=0

	FROM i=1 TO 5000 { PRE[i-1]=0 tspike[i-1]=0 }  
	tspike[0]=1e12	(ms)
	if(tau_1>=tau_rec){ 
		printf("Warning: tau_1 (%g) should never be higher neither equal to tau_rec (%g)!\n",tau_1,tau_rec)
		tau_rec=tau_1+1e-5
	} 
}

FUNCTION diffusione(){	 
	LOCAL DifWave,i	
	DifWave=0
	FROM i=1 TO numpulses{
		tzero=tspike[i-1]
		if(t>tzero){
		    :printf("%g\t",(t-tzero))
			DifWave=DifWave+Mres*exp(-Rd*Rd/(4*Diff*(t-tzero)))/((4*PI*Diff*(1e-3)*lamd)*(t-tzero))^nd
		}
	}	
	diffusione=DifWave
}

BREAKPOINT {

	if ( diffuse && (t>tspike[0]) ) { Trelease= T + diffusione() } else { Trelease=T }

	SOLVE kstates METHOD sparse
	g = gmax * gbar_Q10 * O
	i = (1e-6) * g * (v - Erev) * gmax_factor
	ic = i
}

KINETIC kstates {
	r1 = r1FIX * Trelease^2 / (Trelease + kB)^2
	r6 = r6FIX * Trelease^2 / (Trelease + kB)^2
	~ C  <-> O	(r1*Q10,r2*Q10)
	~ O  <-> D	(r3*Q10,r4*Q10)
	~ D  <-> C	(r5*Q10,r6*Q10)
	CONSERVE C+O+D = 1
}


NET_RECEIVE(weight, on, nspike, tzero (ms)) { :,y, z, u, tsyn (ms)

: *********** ATTENZIONE! ***********
:
: Qualora si vogliano utilizzare impulsi di glutammato saturanti e' 
: necessario che il pulse sia piu' corto dell'intera simulazione
: altrimenti la variabile on non torna al suo valore di default.

INITIAL {
	nspike = 1
    }
    if (flag == 0) { 
		: Qui faccio rientrare la modulazione presinaptica
		nspike = nspike + 1
		if (!on) {
			tzero = t
			on = 1				
			T=Tmax			
			tspike[numpulses] = t
			numpulses=numpulses+1
		}
		net_send(Cdur, nspike)	: !
    }
	if (flag == nspike) { 
			tzero = t
			T = 0
			on = 0
		    }
}
