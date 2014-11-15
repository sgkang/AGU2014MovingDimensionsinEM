from SimPEG import *
import simpegEM as EM
from pymatsolver import MumpsSolver

cs, ncx, ncy, ncz, npad = 20., 30, 20, 30, 12
hx = [(cs,npad,-1.4), (cs,ncx), (cs,npad,1.4)]
hy = [(cs,npad,-1.4), (cs,ncy), (cs,npad,1.4)]
hz = [(cs,npad,-1.4), (cs,ncz), (cs,npad,1.4)]
mesh = Mesh.TensorMesh([hx,hy,hz], 'CCC')
print ("Padding distance x: %10.5f m") % (np.sum(mesh.hx[:npad]))
print ("Padding distance z: %10.5f m") % (np.sum(mesh.hz[:npad]))
print ("Min dx: %10.5f m") % (mesh.hx.min())
print ("Min dz: %10.5f m") % (mesh.hz.min())
print mesh
stop
sigma = Utils.meshutils.readUBCTensorModel('sigma.con', mesh)

x1 = np.arange(30)*10 - 300.
y1 = np.arange(30)*10 - 150.
xyz1 = Utils.ndgrid(x1, y1, np.r_[0.])

x2 = np.arange(30)*10 + 10.
y2 = np.arange(30)*10 - 150.
xyz2 = Utils.ndgrid(x2, y2, np.r_[0.])

mapping = Maps.IdentityMap(mesh)
ntx = 2
time = np.logspace(-4, -2, 31)

rx1 = EM.TDEM.RxTDEM(np.array([[xyz1[:,0], xyz1[:,1], xyz1[:,2]]]), time, 'bz')
tx1 = EM.TDEM.TxTDEM(np.array([0., -150., 0.]), 'CircularLoop_MVP', [rx1])
tx1.radius = 250.
rx2 = EM.TDEM.RxTDEM(np.array([[xyz2[:,0], xyz2[:,1], xyz2[:,2]]]), time, 'bz')
tx2 = EM.TDEM.TxTDEM(np.array([0.,  150., 0.]), 'CircularLoop_MVP', [rx2])
tx2.radius = 250.

survey = EM.TDEM.SurveyTDEM([tx1, tx2])
prb = EM.TDEM.ProblemTDEM_b(mesh, mapping=mapping, verbose=True)
prb.Solver = MumpsSolver
prb.solverOpts = {"symmetric":True}
prb.timeSteps = [(1e-4/10, 15), (1e-3/10, 15), (1e-2/10, 15), (1e-1/10, 15)]
if prb.ispaired:
    prb.unpair()
if survey.ispaired:
    survey.unpair()
prb.pair(survey)

bz = survey.dpred(sigma)
np.save('bz', bz)
