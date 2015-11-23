# -*- coding: utf-8 -*-
"""
Tests the granule cell model by simulating the injection of a current pulse.
Usage: python test.py
"""
from neuron import h
import numpy as np
import sys

def testCell(run, plot, amplitude = 10.e-3):
    from GRANULE_Cell import Grc

    delay = 100.
    duration = 500.
    tstop = 700

    granule = Grc(position=np.zeros(3))
    electrode = h.IClamp(0.5, sec=granule.soma)
    electrode.delay = delay
    electrode.dur = duration
    electrode.amp = amplitude

    h.psection()
    h.celsius = 30
    print "Temperature of simulation is %f deg C "%h.celsius

    # run the simulation
    h.load_file("stdrun.hoc")
    h.tstop = tstop
    h.dt = 0.025
    h.steps_per_ms = 40
    h.v_init = -60

    if run:
        h.run()


        # convert the membrane potential and time values in numpy arrays
        time_points = np.array(granule.time)
        vm_points = np.array(granule.vm)

        spiking  = 0 
        threshold = 0
        spike_times = []

        for i in range(len(time_points)):

            if not spiking and vm_points[i] >= threshold:
                spike_times.append(float(time_points[i]))
                spiking = 1
            elif spiking and vm_points[i] < threshold:
                spiking = 0

        spike_times_file = open("spike_times_%sdeg_%snA.dat"%(h.celsius, electrode.amp), 'w')
        for st in spike_times:
            spike_times_file.write('%f\n'%st)
        spike_times_file.close()

        print "Saved spike times in file %s"%spike_times_file.name
        
        
        voltage_file = open("voltage_%sdeg_%snA.dat"%(h.celsius, electrode.amp), 'w')
        
        for i in range(len(time_points)):
            voltage_file.write('%s\t%s\n'%(time_points[i], vm_points[i]))
        voltage_file.close()

        print "Saved voltage in file %s"%voltage_file.name

        if plot:
            h.run()

            from matplotlib import pyplot as plt
            fig = plt.figure()
            trace_ax = fig.add_subplot(111)
            trace_ax.plot(time_points, vm_points)

            trace_ax.set_xlabel('Time (ms)')
            trace_ax.set_ylabel('mV')
            trace_ax.set_title('Membrane potential')

            fig.suptitle('Pulse duration: %d ms; pulse amplitude: %d pA; temperature: %d deg C' % (duration, amplitude * 1e3, h.celsius))

            plt.show()

    return granule

if __name__ == "__main__":
    
    plot = not '-nogui' in sys.argv
    testCell(True, plot)
