from SimPEG import *
import simpegEM as EM
from pymatsolver import MumpsSolver

cs, ncx, ncy, ncz, npad = 20., 30, 20, 30, 15
hx = [(cs,npad,-1.4), (cs,ncx), (cs,npad,1.4)]
hy = [(cs,npad,-1.4), (cs,ncy), (cs,npad,1.4)]
hz = [(cs,npad,-1.4), (cs,ncz), (cs,npad,1.4)]
mesh = Mesh.TensorMesh([hx,hy,hz], 'CCC')
print ("Padding distance x: %10.5f m") % (np.sum(mesh.hx[:npad]))
print ("Padding distance z: %10.5f m") % (np.sum(mesh.hz[:npad]))
print ("Min dx: %10.5f m") % (mesh.hx.min())
print ("Min dz: %10.5f m") % (mesh.hz.min())

sigma = Utils.meshutils.readUBCTensorModel('sigma_realistic_late.con', mesh)

x1 = np.arange(30)*10 - 300.
y1 = np.arange(30)*10 - 150.
xyz1 = Utils.ndgrid(x1, y1, np.r_[0.])

x2 = np.arange(30)*10 + 10.
y2 = np.arange(30)*10 - 150.
xyz2 = Utils.ndgrid(x2, y2, np.r_[0.])

mapping = Maps.IdentityMap(mesh)
ntx = 2
frequency = np.r_[0.1, 1., 10.]




txList = []

for itrx in range(2):
    if itrx == 0:
        for freq in frequency:
            rxr = EM.FDEM.RxFDEM(xyz1,'bzr')
            rxi = EM.FDEM.RxFDEM(xyz1,'bzi')
            tx = EM.FDEM.TxFDEM(np.array([0., -150., 0.]), 'CircularLoop', freq, [rxr, rxi])
            tx.radius = 250.
            txList.append(tx)
    elif itrx == 1:
        for freq in frequency:
            rxr = EM.FDEM.RxFDEM(xyz2,'bzr')
            rxi = EM.FDEM.RxFDEM(xyz2,'bzi')
            tx = EM.FDEM.TxFDEM(np.array([0., 150., 0.]), 'CircularLoop', freq, [rxr, rxi])
            tx.radius = 250.
            txList.append(tx)

survey = EM.FDEM.SurveyFDEM(txList)

prb = EM.FDEM.ProblemFDEM_b(mesh, mapping=mapping, verbose=True)
prb.Solver = MumpsSolver
prb.solverOpts = {"symmetric":True}

if prb.ispaired:
    prb.unpair()
if survey.ispaired:
    survey.unpair()
prb.pair(survey)

bz = survey.dpred(sigma)
np.save('bz_FD_realistic', bz)
