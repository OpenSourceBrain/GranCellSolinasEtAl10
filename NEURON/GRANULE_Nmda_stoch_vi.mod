NEURON {
	POINT_PROCESS GRANULE_Nmda_stoch_vi
	NONSPECIFIC_CURRENT i
	RANGE Q10_diff,Q10_channel
	RANGE g , ic
	RANGE Cdur,Erev,T,Tmax
	RANGE Rb, Ru, Rd, Rr, Ro, Rc,rb,gmax,RdRate
	RANGE PRE
	RANGE Used
	RANGE MgBlock,v0_block,k_block
	RANGE diffuse,Trelease,lamd, Diff, M, Rd, nd, syntype, gmax_factor
	RANGE T_factor
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
	: Parametri Presinaptici
	gmax		= 18800  	(pS)	: 7e3 : 4e4
	Q10_diff	= 1.5
	Q10_channel	= 2.4

	M	= 21.515 : numero di (kilo) molecole in una vescicola		
	Rd	= 1.03 (um)
	Diff	= 0.223 (um2/ms)
	
	Tmax		= 1  	(mM)
	T_factor = 0.5


	: Postsinaptico, Westbrook scheme
	
	Cdur	= 0.3	(ms)
	Rb	=  5		(/ms/mM)  	: binding
	Ru	=  0.1		(/ms)		: unbinding
	RdRate	=  12e-3  	(/ms)		: desensitization
	Rr	=  9e-3		(/ms)		: resensitization
	Ro	=  3e-2 	(/ms)		: opening
	Rc	=  0.966	(/ms)		: closing
	Erev	= -3.7  (mV)	: 0 (mV)
	
	v0_block = -20 (mV)	: -16 -8.69 (mV)	: -18.69 (mV) : -32.7 (mV)
	k_block  = 13 (mV)
	nd	 = 1
	: kB	 = 0.44	(mM)

	: Diffusion			
	diffuse	= 1
	lamd	= 20 (nm)
	celsius (degC)
}


ASSIGNED {
	v		(mV)		: postsynaptic voltage
	i 		(nA)		: current = g*(v - Erev)
	ic 		(nA)		: current = g*(v - Erev)
	g 		(pS)		: actual conductance

	rb		(/ms)    : binding
	
	T		(mM)
	x 
	
	Trelease	(mM)
	tspike[100]	(ms)	: will be initialized by the pointprocess
	PRE[100]
	Mres		(mM)	
	
	MgBlock
	numpulses
	tzero
	gbar_Q10 (mho/cm2)
	Q10 (1)
}

STATE {
	: Channel states (all fractions)
	C0		: unbound
	C1		: single bound
	C2		: double bound
	D		: desensitized
	O		: open
}

INITIAL {
	rates(v)
	C0 = 1
	C1 = 0
	C2 = 0
	D  = 0
	O  = 0
	T  = 0
	numpulses=0

	gbar_Q10 = Q10_diff^((celsius-30)/10)
	Q10 = Q10_channel^((celsius-30)/10)

	Mres = 1e3 * (1e3 * 1e15 / 6.022e23 * M)     : (M) to (mM) so 1e3, 1um^3=1dm^3*1e-15 so 1e15
	FROM i=1 TO 100{ PRE[i-1]=0 tspike[i-1]=0 } :PRE_2[500]=0}
	tspike[0]=1e12	(ms)

}

FUNCTION imax(a,b) {
    if (a>b) { imax=a }
    else { imax=b }
}

FUNCTION diffusione(){	 
	LOCAL DifWave,i,cntc,fi
	DifWave=0
	cntc=imax(numpulses-100,0)
	FROM i=cntc  TO numpulses{
	    fi=fmod(i,100)
		tzero=tspike[fi]
		if(t>tzero){
			DifWave=DifWave+PRE[fi]*Mres*exp(-Rd*Rd/(4*Diff*(t-tzero)))/((4*PI*Diff*(1e-3)*lamd)*(t-tzero))^nd
		}
	}    
	diffusione=DifWave
}

BREAKPOINT {
	rates(v)
	SOLVE kstates METHOD sparse	
	g = gmax * gbar_Q10 * O * gmax_factor
	
	: E' piu' logico spostare * MgBlock * PRE sul calcolo della corrente!
	i = (1e-6) * g * (v - Erev) * MgBlock
	ic = i
    }

KINETIC kstates {	
	:if ( diffuse && (t>tspike[0]) ) { Trelease= T + diffusione() } else { Trelease=T }
	Trelease = diffusione()
	rb = Rb * Trelease	:^2 / (Trelease + kB)^2 
	~ C0 <-> C1	(rb*Q10,Ru*Q10) 	: (fattore*rb,Ru) qui 2* per descrizione part.identiche
	~ C1 <-> C2	(rb*Q10,Ru*Q10)		: (rb,fattore*Ru)	idem
	~ C2 <-> D	(RdRate*Q10,Rr*Q10)
	~ C2 <-> O	(Ro*Q10,Rc*Q10)
	CONSERVE C0+C1+C2+D+O = 1
}

PROCEDURE rates(v(mV)) {
	: E' necessario includere DEPEND v0_block,k_block per aggiornare le tabelle!
	TABLE MgBlock DEPEND v0_block,k_block FROM -120 TO 30 WITH 150
	MgBlock = 1 / ( 1 + exp ( - ( v - v0_block ) / k_block ) )
}


NET_RECEIVE(weight, on, nspike, tzero (ms)) {LOCAL fi :,y, z, u, tsyn (ms)

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
			T=Tmax * T_factor :* y			
			fi=fmod(numpulses,100)
			PRE[fi]= T_factor :y :(use with pointer)	: PRE[numpulses]=y
			tspike[fi] = t
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
