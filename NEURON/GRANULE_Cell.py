# -*- coding: utf-8 -*-
from neuron import h
from synapse.synapse import Synapse


class Grc:
    def __init__(self,position):
        self.soma = h.Section(name='soma')
        self.soma.nseg = 1
        self.soma.diam = 9.76
        self.soma.L = 9.76
        self.soma.cm = 1
        self.soma.Ra = 100

        #Shouldn't change a global param like this in a class...
        #h.celsius = 37

        self.whatami = "grc"

        self.soma.push()
        h.pt3dclear()
        h.pt3dadd(position.item(0), position.item(1) - self.soma.L, position.item(2), self.soma.diam)
        h.pt3dadd(position.item(0), position.item(1) + self.soma.L, position.item(2), self.soma.diam)
        h.pop_section()

        self.soma.insert('GRANULE_LKG1')
        self.soma.insert('GRANULE_LKG2')
        self.soma.insert('GRANULE_Nmda_leak')
        self.soma.insert('GRANULE_NA')
        self.soma.insert('GRANULE_NAR')
        self.soma.insert('GRANULE_PNA')
        self.soma.insert('GRANULE_KV')
        self.soma.insert('GRANULE_KA')
        self.soma.insert('GRANULE_KIR')
        self.soma.insert('GRANULE_KCA')
        self.soma.insert('GRANULE_KM')
        self.soma.insert('GRANULE_CA')
        self.soma.insert('GRANULE_CALC')

        h.usetable_GRANULE_NA = 1
        h.usetable_GRANULE_NAR = 1
        h.usetable_GRANULE_PNA = 1
        h.usetable_GRANULE_KV  = 1
        h.usetable_GRANULE_KA = 1
        h.usetable_GRANULE_KIR = 1
        h.usetable_GRANULE_KCA = 0
        h.usetable_GRANULE_KM = 1
        h.usetable_GRANULE_CA = 1

        self.soma.ena = 87.39
        self.soma.ek = -84.69
        self.soma.eca = 129.33

        self.MF_L = []
        self.GOC_L = []
        self.mfncpc = []
        self.gocncpc = []

        self.spike = h.NetStim(0.5, sec= self.soma)
        self.spike.start = -10
        self.spike.number = 1
        self.spike.interval = 1e9
        self.nc = h.NetCon(self.soma(1)._ref_v, self.spike, sec = self.soma)
        self.nc.threshold=-20
        self.nc.delay=0
        self.nc.weight[0]=1

        self.SpikeTrain_output = [h.Vector(),h.Vector()]
        self.nc.record(self.SpikeTrain_output[1], self.SpikeTrain_output[0], 1)

        self.vm = h.Vector()
        self.vm.record(self.soma(.5)._ref_v, sec = self.soma)
        self.time = h.Vector()
        self.time.record(h._ref_t)


    # To be removed in the OOP version
    ## def connect2target(self, target):
    ##     nc = h.NetCon(self.soma(1)._ref_v, target, sec = self.soma)
    ##     nc.threshold = 10
    ##     return nc

    #Synapses
    def createsyn(self,nmf,nrel = 0,ngoc = -1):
        # Use here the source target sting name
        # so the presynaptic link is not made
        # and it will have to be manged later
        # by the gid connect for parallel simulations
        #Mossy
        if ngoc <0 :
            ngoc = nmf

        for i in range(nmf):
            self.MF_L.append(Synapse('glom',self,self.soma,nrel))

            #Inibition
        for i in range(ngoc):
            self.GOC_L.append(Synapse('goc',self,self.soma))


	
    def pconnect(self,pc,source,syn_idx,type):
        if type == 'mf':
            source = int(source)
            # print len(self.MF_L), syn_idx, self.whatami
            self.mfncpc.append(pc.gid_connect(source, self.MF_L[syn_idx].input))
            self.mfncpc[-1].delay = 1
            self.mfncpc[-1].weight[0] = 1
            return self.mfncpc[-1]
        if type == 'goc':
            source = int(source)
            self.gocncpc.append(pc.gid_connect(source, self.GOC_L[syn_idx].input))
            self.gocncpc[-1].delay = 1
            self.gocncpc[-1].weight[0] = 1
            return self.gocncpc[-1]

    def getSomaSection(self):
        return self.soma

    def destroy(self):
        del self.nc
        for m in self.MF_L:
            m.destroy()
            del m
        for m in self.GOC_L:
            m.destroy()
            del m
        for m in self.mfncpc:
            del m
        for m in self.gocncpc:
            del m
        del self
