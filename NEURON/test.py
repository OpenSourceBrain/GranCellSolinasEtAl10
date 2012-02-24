# -*- coding: utf-8 -*-
"""
Tests the granule cell model by simulating the injection of a current pulse.
Usage: python test.py
"""
from neuron import h
import numpy as np
from matplotlib import pyplot as plt

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

# run the simulation
h.load_file("stdrun.hoc")
h.tstop = tstop
h.dt = 0.025
h.steps_per_ms = 40

h.run()

# convert the membrane potential and time values in numpy arrays
time_points = np.array(granule.time)
vm_points = np.array(granule.vm)

# plot
fig = plt.figure()
trace_ax = fig.add_subplot(111)
trace_ax.plot(time_points, vm_points)

trace_ax.set_xlabel('Time (ms)')
trace_ax.set_ylabel('mV')
trace_ax.set_title('Membrane potential')

fig.suptitle('pulse duration: %d ms; pulse amplitude: %d pA' % (duration, amplitude * 1e3))

plt.show()
