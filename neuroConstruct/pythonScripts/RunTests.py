#
#
#   File to test current configuration of GranCellSolinasEtAl10 project.
#
#   To execute this type of file, type 'nC.bat -python XXX.py' (Windows)
#   or 'nC.sh -python XXX.py' (Linux/Mac). Note: you may have to update the
#   NC_HOME and NC_MAX_MEMORY variables in nC.bat/nC.sh
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

sys.path.append(os.environ["NC_HOME"]+"/pythonNeuroML/nCUtils")

import ncutils as nc # Many useful functions such as SimManager.runMultipleSims found here

projFile = File("../GranCellSolinasEtAl10.ncx")


##############  Main settings  ##################

simConfigs = []

simConfigs.append("Original mod channels")
#simConfigs.append("AllChansBothSims")

simDt =                 0.001

simulators =            ["NEURON"]

varTimestepNeuron =     False
varTimestepTolerance =  0.0001

plotSims =              True
plotVoltageOnly =       True
runInBackground =       True
analyseSims =           True
verbose =               True

#############################################


def testAll(argv=None):
    if argv is None:
        argv = sys.argv

    print "Loading project from "+ projFile.getCanonicalPath()


    simManager = nc.SimulationManager(projFile,
                                      verbose = verbose)

    simManager.runMultipleSims(simConfigs =           simConfigs,
                               simDt =                simDt,
                               simulators =           simulators,
                               runInBackground =      runInBackground,
                               varTimestepNeuron =    varTimestepNeuron,
                               varTimestepTolerance = varTimestepTolerance)

    simManager.reloadSims(plotVoltageOnly =   plotVoltageOnly,
                          analyseSims =       analyseSims)

    # These were discovered using ../../NEURON/test.py WITH DT = 0.001

    spikeTimesToCheck = {'GranCell_mod_0': [125.484, 150.146, 174.864, 199.693, 224.589, 249.517, 274.456, 299.397, 324.338, 349.276, 374.212, 399.145, 424.077, 449.008, 473.937, 498.865, 523.792, 548.719, 573.645, 598.57]}
    
    spikeTimeAccuracy = 0.01

    report = simManager.checkSims(spikeTimesToCheck = spikeTimesToCheck,
                                  spikeTimeAccuracy = spikeTimeAccuracy)

    print report

    return report


if __name__ == "__main__":
    testAll()


