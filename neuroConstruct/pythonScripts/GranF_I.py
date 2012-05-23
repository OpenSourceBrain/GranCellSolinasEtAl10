#
#
#   A file which generates a frequency vs current curve for the Granule cell
#
#   To execute this file, type 'nC.bat -python GranF_I.py' (Windows)
#   or 'nC.sh -python GranF_I.py' (Linux/Mac)
#
#   Author: Padraig Gleeson
#
#   This file has been developed as part of the neuroConstruct project
#   This work has been funded by the Wellcome Trust
#
#

import sys
import os

try:
	from java.io import File
except ImportError:
	print "Note: this file should be run using nC.bat -python XXX.py' or 'nC.sh -python XXX.py'"
	print "See http://www.neuroconstruct.org/docs/python.html for more details"
	quit()

from math import *

sys.path.append(os.environ["NC_HOME"]+"/pythonNeuroML/nCUtils")

import ncutils as nc
from ucl.physiol.neuroconstruct.hpc.mpi import MpiSettings

simConfig="Original mod channels"


preStimDel = 0
preStimDur = 200

stimAmpLow =  0.000
stimAmpInc =  0.002
stimAmpHigh = 0.03

stimDel = preStimDur
stimDur = 1500

simDuration = preStimDur + stimDur # ms

analyseStartTime = stimDel + 100 # So it's firing at a steady rate...
analyseStopTime = simDuration
analyseThreshold = -20 # mV

mpiConfig =            MpiSettings.LOCAL_SERIAL    # Default setting: run on one local processor
#mpiConfig =            MpiSettings.MATLEM_1PROC    # Run on one processor on UCL cluster

numConcurrentSims = 4
if mpiConfig != MpiSettings.LOCAL_SERIAL: numConcurrentSims = 30
suggestedRemoteRunTime = 9   # mins


# Load neuroConstruct project

projFile = File("../GranCellSolinasEtAl10.ncx")


simManager = nc.SimulationManager(projFile,
                                  numConcurrentSims)

simManager.generateFICurve("NEURON",
                           simConfig,
                           stimAmpLow,
                           stimAmpInc,
                           stimAmpHigh,
                           stimDel,
                           stimDur,
                           simDuration,
                           analyseStartTime,
                           analyseStopTime,
                           analyseThreshold,
                           mpiConfig =                mpiConfig,
                           suggestedRemoteRunTime =   suggestedRemoteRunTime)