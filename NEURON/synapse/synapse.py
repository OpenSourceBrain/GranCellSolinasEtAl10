# -*- coding: utf-8 -*-
from neuron import h

class Synapse:
    def __init__(self,source,target,section,nrel = 0,syntype = 'ANK',record_all = 1, weight = 1):
		
	self.input = h.NetStim(0.5)
	self.input.start = -10
	self.input.number = 1
	self.input.interval = 1e9
        self.weight = weight

        self.nrel = nrel

        self.postsyns = {}

        if (type(source) == type('s')):
            sourcetype = source
        else:
            sourcetype = source.whatami
            
        if record_all:
            self.SpikeTrain_input = [h.Vector(),h.Vector()]
            self.netcon_in = h.NetCon(self.input,None, 0, 0, 1)
            self.netcon_in.record(self.SpikeTrain_input[1], self.SpikeTrain_input[0], 1)
           
        # Decide what kind of connection to create based on the source and target types
        if sourcetype == 'glom':
            if target.whatami == 'PF':
                ## No input to the fake grcs
                print 'Pf'
            else:
                # Make a mf synapse
                if self.nrel>0 :
                    if target.whatami == 'grc':
                        self.whatami = "syn_glom2grc_stoch"
                    elif target.whatami == 'goc':
                        self.whatami = "syn_glom2goc_stoch"
                    # Use stichastic synapses
                    self.release_sites = [h.Release_site(0.5, sec=section) for i in range(self.nrel)]
                    # Set Prob release
                    for site in self.release_sites:
                        site.U = 0.42
                    # Set netcon for spikes
                    self.netcon_out = [h.NetCon(release_site,None, 0, 0, 1) for release_site in self.release_sites]
                    # Connect input
                    self.nc_rel = [h.NetCon(self.input,release_site,0,0,1) for release_site in self.release_sites]
                    # Add poststnaptic sites
                    if 'A'  in syntype:
                        # AMPA
                        self.postsyns['AMPA'] = [h.GRANULE_Ampa_stoch_vi(.5, sec=section) for r in self.release_sites]
                        if target.whatami == 'goc':
                            for p in self.postsyns['AMPA']: # set the correct gmax for MF-> GoC
                                p.gmax = 4270

                    if 'N'  in syntype:
                        # NMDA
                        self.postsyns['NMDA'] = [h.GRANULE_Nmda_stoch_vi(.5, sec=section) for r in self.release_sites]
                        if target.whatami == 'goc':
                            for p in self.postsyns['NMDA']: # set the correct gmax for MF-> GoC
                                p.gmax = 46500

                    for rec_site in self.postsyns['AMPA']:
                        rec_site.gmax_factor = (1./0.5)/self.nrel
                    for rec_site in self.postsyns['NMDA']:
                        rec_site.gmax_factor = (1./0.5)/self.nrel
                        
                    # Connect pre to post synapse
                    self.nc_syn = [[h.NetCon(release_site, receptor[k],0,0,1) for k, release_site in enumerate(self.release_sites)] for receptor in self.postsyns.itervalues()]

                    # If required recod the timing of each vescicle release
                    if record_all:
                        self.SpikeTrain_output = [[h.Vector(),h.Vector()] for n in self.netcon_out]
                        for i,n in enumerate(self.netcon_out):
                            n.record(self.SpikeTrain_output[i][1],self.SpikeTrain_output[i][0],i+2)
                                        
                else:
                    if target.whatami == 'grc':
                        self.whatami = "syn_glom2grc_det"
                    elif target.whatami == 'goc':
                        self.whatami = "syn_glom2goc_det"
                    # Use deterministic synapses
                    if 'A'  in syntype:
                        # AMPA
                        self.postsyns['AMPA'] = [h.GRANULE_Ampa_det_vi(0.5, sec=section)]
                        if target.whatami == 'goc':
                            self.postsyns['AMPA'][0].gmax = 4270
                            
                    if 'N'  in syntype:
                        # NMDA
                        self.postsyns['NMDA'] = [h.GRANULE_Nmda_det_vi(0.5, sec=section)]
                        if target.whatami == 'goc':
                            self.postsyns['NMDA'][0].gmax = 46500

                    # Connect input to the receptors
                    self.nc_syn = [h.NetCon(self.input,receptor[0],0,0,1) for receptor in self.postsyns.itervalues()]


        elif sourcetype == 'goc':
            # Make a Golgi (GABAergic) synapse onto a granule cell
            # Use deterministic synapses
            self.whatami = "syn_goc2grc_det"
            self.postsyns['GABA'] = [h.GRANULE_Gaba_det_vi(0.5, sec=section)]
            self.nc_syn = [h.NetCon(self.input,receptor[0],0,0,1) for receptor in self.postsyns.itervalues()]

                
        elif sourcetype == 'grc':
            if target.whatami == 'goc':
                # Make a PF synapse onto a golgi cell
                # Use deterministic synapses
                self.whatami = "syn_grc2goc_det"
                if 'A'  in syntype:
                    # AMPA
                    self.postsyns['AMPA'] = [h.GRANULE_Ampa_det_vi(0.5, sec=section)]
                    self.postsyns['AMPA'][0].U = 0.1
                    self.postsyns['AMPA'][0].diffuse = 0
                    self.postsyns['AMPA'][0].gmax = 12000

                if 'N'  in syntype:
                    # NMDA
                    self.postsyns['NMDA'] = [h.GRANULE_Nmda_det_vi(0.5, sec=section)]
                    self.postsyns['NMDA'][0].gmax = 18800*8
                    self.postsyns['NMDA'][0].U = 0.1

                if 'K'  in syntype:
                    # KAINATE
                    self.postsyns['KAIN'] = [h.GRANULE_Ampa_det_vi(0.5, sec=section)]
                    self.postsyns['KAIN'][0].U = 0.1
                    self.postsyns['KAIN'][0].Cdur = 0.3/5
                    self.postsyns['KAIN'][0].r1FIX = 0.1
                    self.postsyns['KAIN'][0].r2 = 0.01
                    self.postsyns['KAIN'][0].r6FIX = 0.1
                    self.postsyns['KAIN'][0].gmax = 19000/20 
                self.nc_syn = [h.NetCon(self.input,receptor[0],0,0,1) for receptor in self.postsyns.itervalues()]

            elif target.whatami == 'stl':
                self.postsyns['AMPA'] = [h.GRANULE_Ampa_det_vi(0.5, sec=section)]
                self.postsyns['AMPA'][0].tau_facil=10.8*5
                self.postsyns['AMPA'][0].tau_rec=35.1*1
                self.postsyns['AMPA'][0].tau_1=3*5
                self.postsyns['AMPA'][0].gmax = 9500 #*area(0.5)/1.33 // 7.5   2050  3943 o 2790 valori a 25 gradia 37 gradi,1990 la media
                self.postsyns['AMPA'][0].U=0.15				
                self.nc_syn = [h.NetCon(self.input,receptor[0],0,0,1) for receptor in self.postsyns.itervalues()]
                

        elif sourcetype == 'stl':
            if target.whatami == 'stl':
                self.postsyns['GABA'] = [h.ExpSyn(0.5, sec=section)] # self.postsyns
                self.postsyns['GABA'][0].tau = 1.1
                self.postsyns['GABA'][0].e = -60
                self.nc_syn = [h.NetCon(self.input,receptor[0],0,0,1) for receptor in self.postsyns.itervalues()]
                for nc in self.nc_syn:
                    nc.weight[0] = 5e-2

#         if (type(sourcetype) != type('s')):
#             return self.nc_input

    def prel(self,prel): 
        for r in self.release_sites:
            r.U = 0.42

    def destroy(self):
        del self.netcon_in
        if self.nrel > 0:
            for r in self.netcon_out:
                del r
            for r in self.nc_rel:
                del r
        else:
            for r in self.nc_syn:
                del r
        del self
