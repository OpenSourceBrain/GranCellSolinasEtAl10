## Granule cell - Solinas et al. 2010

[![Continuous build using OMV](https://github.com/OpenSourceBrain/GranCellSolinasEtAl10/actions/workflows/omv-ci.yml/badge.svg)](https://github.com/OpenSourceBrain/GranCellSolinasEtAl10/actions/workflows/omv-ci.yml)

Initial version of Granule cell from: Solinas S., Nieus T, d'Angelo E. (2010) A Realistic Large-Scale Model of the Cerebellum Granular Layer Predicts Circuit Spatio-Temporal Filtering Properties. Front Cell Neurosci. 2010;4:12.

See http://www.opensourcebrain.org/projects/grancellsolinasetal10 for more details

### Installation

Details on getting a local copy of the project and running the NEURON & neuroConstruct/NeuroML versions.

**Clone the project**

To get a local copy of this project, [install Git](http://www.opensourcebrain.org/projects/gitintro/wiki/Wiki) and type:

    git clone https://github.com/OpenSourceBrain/GranCellSolinasEtAl10
    cd GranCellSolinasEtAl10

**Original NEURON version**

The PyNEURON version of this model is available [here](https://github.com/OpenSourceBrain/GranCellSolinasEtAl10/tree/master/NEURON). More details on NEURON are [here](http://www.opensourcebrain.org/projects/simulators/wiki/Wiki#NEURON).

Once the project is checked out as above, run the NEURON version with:

    cd NEURON
    nrnivmodl
    nrngui -python test.py

**neuroConstruct version**

The NeuroML conversion of this project can be executed using neuroConstruct.

See full instructions for installing neuroConstruct, and for running this and other projects from the OSB [here](http://opensourcebrain.org/docs#Using_neuroConstruct_Based_Projects).

**Granule cell properties**

See here: [Cerebellar Granule cell modelling](http://www.opensourcebrain.org/projects/cerebellar-granule-cell-modelling/wiki/Wiki)


