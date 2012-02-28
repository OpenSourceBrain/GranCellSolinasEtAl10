# -*- coding: utf-8 -*-
"""
Tests the granule cell model by simulating the injection of a current pulse.
Usage: python test.py
"""
from neuron import h
import numpy as np

def testCell(runAndPlotAlso):
    from GRANULE_Cell import Grc

    delay = 100.
    duration = 500.
    amplitude = 10.e-3
    tstop = 700

    granule = Grc(position=np.zeros(3))
    electrode = h.IClamp(0.5, sec=granule.soma)
    electrode.delay = delay
    electrode.dur = duration
    electrode.amp = amplitude

    h.psection()
    h.celsius = 37
    print "Temperature of simulation is %f deg C "%h.celsius

    # run the simulation
    h.load_file("stdrun.hoc")
    h.tstop = tstop
    h.dt = 0.025
    h.steps_per_ms = 40
    h.v_init = -60

    if runAndPlotAlso:
        h.run()

        # convert the membrane potential and time values in numpy arrays
        time_points = np.array(granule.time)
        vm_points = np.array(granule.vm)

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
    testCell(True)