from SimPEG import Directives
import numpy as np
class RememberXC(Directives.InversionDirective):
    xcs = []
    dpred = []
    beta = []
    objdata = []
    objmodel = []
    invfile = None

    def endIter(self):
        # Need to ask a way that I can get Initial misfit
        self.xcs.append(self.opt.xc)
        self.dpred.append(self.invProb.dpred)
        self.beta.append(self.invProb.beta)
        self.objdata.append(self.invProb.phi_d)
        self.objmodel.append(self.invProb.phi_m)
        if self.invProb.opt.iter==1:
            self.invfile = open("SimPEGinv.out", "w")
            fname = 'dobs'
            np.save(fname, self.invProb.survey.dobs)
            Header = 'iter'+'   phi_d'+'       phi_m'+'       beta  \n'
            self.invfile.writelines(Header)
            self.invfile.writelines('===================================================\n')

        self.invfile.writelines(("%3s %12.5e %10.5e %10.5e \n") % (self.invProb.opt.iter, self.invProb.phi_d, self.invProb.phi_m, self.invProb.beta))
        fname = 'dpred_'+str(self.invProb.opt.iter)
        np.save(fname, self.invProb.dpred)
        fname = 'model_'+str(self.invProb.opt.iter)
        np.save(fname, self.opt.xc)

        if self.invProb.opt.iter==self.invProb.opt.maxIter:
            self.invfile.close()
