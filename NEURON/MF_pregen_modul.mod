: $Id: pregen.mod,v 1.3 2000/05/16 11:16:56 hines Exp $
: comments at end

NEURON	{ 
    POINT_PROCESS MF_SpikeGenerator_SS
    RANGE y, burst, burst_on, burst_off
	: burst stores the in_burst state, it is set to 1 when event_t == the_first_spike_of_burst and to 0 when event_t == the_firs_spike_after_the_burst
    RANGE IntraBurstInterval, InterBurstInterval, burst_len, start_background, end,delay
    RANGE BackgroundNoise, InterBurstIntervalNoise,tonic,start_burst, IntraBurstNoise, Stim
    :,sin_freq, sin_phase,sin_amp
}

UNITS {
    PI   = (pi)(1)
}
    
    PARAMETER {
    IntraBurstInterval	= 10 (ms)	: time between spikes in a burst (msec)
    IntraBurstNoise = 0                 : randomness of spike times inside a burst
    InterBurstInterval	= 1000 (ms)	: burst period (msec)
    : actually, above is interburst period in conformity with original version
    : see
    burst_len	= 2		: burst length (# spikes)
    start_background		= 50 (ms)	: start of first interburst interval
    start_burst		= 100 (ms)	: start of first burst
    end		= 1e10 (ms)	: time to stop bursting
    BackgroundNoise		= 0		: amount of randomeaness (0.0 - 1.0)
    InterBurstIntervalNoise		= 0		: InterBurstIntervalNoiseamount of randomeaness (0.0 - 1.0)
    delay		= 0
    
    tonic = 0
    sin_amp = 0
    sin_freq = 10
    sin_phase = 0
}

ASSIGNED {
    y
    last_event_t
    burst
    set_burst
    event (ms)
    burst_off (ms)
    burst_on (ms)
    next_burst_on (ms)
    toff (ms)
    on
    inter_burst_ISI (ms)
    Stim[1000]
    cnt
}

PROCEDURE seed(x) {
    set_seed(x)
}

INITIAL {
    on = 1
    toff = 1e9
    y = -90
    burst = 0
    burst_on = start_burst
    Stim[0]=burst_on
    cnt = 0
    next_burst_on = burst_on + interval(InterBurstInterval,InterBurstIntervalNoise)
    set_burst = 1
    event = start_background : - InterBurstInterval
:    if (tonic-sin_amp/2 <= 0) {
:	tonic = 1e-6
:	sin_amp = 0
:    }
    :
    event_time()
    while (on == 1 && event < 0) {
	event_time()
    }
    if (on == 1) {
	net_send(event, 1)
    }
}	

FUNCTION interval(mean (ms), noise (1)) (ms) {
    if (mean <= 0.) {
	mean = .01 (ms) : I would worry if it were 0.
	: since mean is a local variable, if the number it is set
	: to is dimensionless, mean will be dimensionless.
    }
    if (noise == 0) {
	interval = mean
    }else{
	interval = (1. - noise)*mean + noise*(mean*exprand(1)+delay)
    }
}

FUNCTION ISI_modulated(tonic (Hz), sin_amp (Hz), sin_freq (Hz), sin_phase (1), event (ms)) (ms) {
:	printf("sum %g\n",tonic + sin_amp * sin(2*PI * sin_freq * event/1000 + sin_phase))
:	if (tonic + sin_amp * sin(2*PI * sin_freq * event/1000 + sin_phase) <= 0) {
:	    if (cos(2*PI * sin_freq * event/1000 + sin_phase) <= 0) {
:		ISI_modulated = 1000*((PI/2 + asin(tonic/sin_amp)-sin_phase))/(2*PI*sin_freq)
:	    } else {
:		ISI_modulated = 1000*(asin(tonic/sin_amp)-sin_phase)/(2*PI*sin_freq)
:	    }		
:	printf("fore %g\n",ISI_modulated)
:	printf("asin %g\n",asin(tonic/sin_amp))
	
:    } else {
	ISI_modulated = 1000./(tonic + sin_amp * sin(2*PI * sin_freq * event/1000 + sin_phase))
:    }
}

PROCEDURE event_time() {
:    printf("show burst=%g burst_on=%g burst_off=%g event=%g t=%g\n",burst,burst_on,burst_off,event,t)
    if (InterBurstInterval == 0 || (burst != 0. && burst_len > 1)) {
	: We are inside a burst
	event = event + interval(IntraBurstInterval,IntraBurstNoise)
:	printf("BurstSpike event=%g\n",event)
	if (event > burst_on + burst_off) {
	    burst = 0.
	    set_burst = 0
	}
    }else{
	: We are between bursts
	: if InterBurstInterval from beginning of burst to beginning of burst
	:event = event + interval(InterBurstInterval - (burst_len-1)*IntraBurstInterval,InterBurstIntervalNoise)
	: use following if InterBurstInterval is interburst interval
	
	if (burst_len == 1 && burst == 1.) {
	    burst = 0.
	    set_burst = 0
	}
	burst_off = interval((burst_len - 1)*IntraBurstInterval,IntraBurstNoise)-1e-6
	if (set_burst == 0) {
	    burst_on = next_burst_on
	    next_burst_on = burst_on + interval(InterBurstInterval,InterBurstIntervalNoise)
	    set_burst = 1
	}
	
	inter_burst_ISI = ISI_modulated(tonic, sin_amp, sin_freq, sin_phase, event)
	
	if (start_background < start_burst || event > start_background) {
	    event = event + interval(inter_burst_ISI,BackgroundNoise)
	} else {
	    event = burst_on
	}
	
:	printf("InterBurst event=%g inter_burst_ISI=%g\n",event,inter_burst_ISI)
	
	if (event >= burst_on) {
	    :printf("begin burst")
	    burst = 1.
	    event = burst_on
            Stim[cnt]  = burst_on
            Stim[cnt+1] = burst_on + burst_off
            cnt = cnt + 2
	}
    }
    if (event > end) {
	on = 0
    }
}

NET_RECEIVE (w) {
:    printf("Pregen receive t=%g flag=%g\n", t, flag) 
    if (flag == 1 && on == 1) {
	y = 20
	net_event(t)
	event_time()
	last_event_t = t
	net_send(event - t, 1)
	net_send(.1, 2)
    }
    if (flag == 2) {
	y = -90
    }
}

COMMENT
Presynaptic spike generator
---------------------------

This mechanism has been written to be able to use synapses in a single
neuron receiving various types of presynaptic trains.  This is a "fake"
presynaptic compartment containing a fast spike generator.  The trains
of spikes can be either periodic or noisy (Poisson-distributed), and 
either tonic or bursting.

Parameters;
BackgroundNoise: 	between 0 (no noise-periodic) and 1 (fully noisy)
IntraBurstInterval: 	fast interval, mean time between spikes (ms)
InterBurstInterval:	slow interval, mean burst silent period (ms), 0=tonic train
burst_len: 	mean burst length (nb. spikes)

Written by Z. Mainen, modified by A. Destexhe, The Salk Institute

Modified by Michael Hines for use with CVode

Modified by Michael Hines to use logical event style with NET_RECEIVE
ENDCOMMENT

