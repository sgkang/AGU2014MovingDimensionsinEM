from SimPEG import *
import simpegEM as EM

from pymatsolver import MumpsSolver
from Rules import RememberXC

cs, ncx, ncy, ncz, npad = 20., 30, 20, 30, 12
hx = [(cs,npad,-1.4), (cs,ncx), (cs,npad,1.4)]
hy = [(cs,npad,-1.4), (cs,ncy), (cs,npad,1.4)]
hz = [(cs,npad,-1.4), (cs,ncz), (cs,npad,1.4)]
mesh = Mesh.TensorMesh([hx,hy,hz], 'CCC')
print ("Padding distance x: %10.5f m") % (np.sum(mesh.hx[:npad]))
print ("Padding distance z: %10.5f m") % (np.sum(mesh.hz[:npad]))
print ("Min dx: %10.5f m") % (mesh.hx.min())
print ("Min dz: %10.5f m") % (mesh.hz.min())

x1 = np.arange(30)*10 - 300.
y1 = np.arange(30)*10 - 150.
xyz1 = Utils.ndgrid(x1, y1, np.r_[0.])

x2 = np.arange(30)*10 + 10.
y2 = np.arange(30)*10 - 150.
xyz2 = Utils.ndgrid(x2, y2, np.r_[0.])

ind = xyz1[:,1] == 0.

ntx = 2
time = np.logspace(-4, -2, 31)

rx1 = EM.TDEM.RxTDEM(xyz1[ind,:], time, 'bz')
tx1 = EM.TDEM.TxTDEM(np.array([-150., 0.,  0.]), 'CircularLoop_MVP', [rx1])
tx1.radius = 250.
rx2 = EM.TDEM.RxTDEM(xyz2[ind,:], time, 'bz')
tx2 = EM.TDEM.TxTDEM(np.array([150., 0., 0.]), 'CircularLoop_MVP', [rx2])
tx2.radius = 250.


mesh2D = Mesh.TensorMesh([hx,hz], 'CC')
active = mesh2D.gridCC[:,1] < 0.
actMap = Maps.ActiveCells(mesh2D, active, np.log(1e-8), nC=mesh2D.nC)
map2to3 = Maps.Map2Dto3D(mesh, normal = 'Y')
mapping =  Maps.ExpMap(mesh) * map2to3 * actMap

survey = EM.TDEM.SurveyTDEM([tx1, tx2])
prb = EM.TDEM.ProblemTDEM_b(mesh, mapping=mapping, verbose=True)
prb.Solver = MumpsSolver
prb.solverOpts = {"symmetric":True}
prb.timeSteps = [(1e-4/10, 15), (1e-3/10, 15), (1e-2/10, 15)]
if prb.ispaired:
    prb.unpair()
if survey.ispaired:
    survey.unpair()
prb.pair(survey)

dobs = np.load('bzobs_realistic_line.npy')
std = 0.05
survey.dobs = Utils.mkvc(dobs)
survey.std = survey.dobs*0 + std
sig_half = 2e-3

sigma = np.ones(mesh2D.nC)*sig_half

m0 = np.log(sigma[active])
dmis = DataMisfit.l2_DataMisfit(survey)
dmis.Wd = 1/(abs(survey.dobs)*std)
opt = Optimization.InexactGaussNewton(maxIter = 15)
reg = Regularization.Tikhonov(mesh, mapping=mapping)
betaest = Directives.BetaEstimate_ByEig(beta0_ratio = 1e2)
beta = Directives.BetaSchedule(coolingFactor = 8, coolingRate = 3)
invprb = InvProblem.BaseInvProblem(dmis, reg, opt)
inv = Inversion.BaseInversion(invprb, directiveList = [betaest,beta, RememberXC()])
#inv = Inversion.BaseInversion(invprb, directiveList = [beta, RememberXC()])


reg.alpha_s = 1e-5
reg.alpha_x = 1.
reg.alpha_y = 1.
reg.alpha_z = 1.
reg.mref = m0
C =  Utils.Counter()
prb.counter = C
opt.counter = C
opt.LSshorten = 0.5

mopt = inv.run(m0)
opt.counter.summary()
