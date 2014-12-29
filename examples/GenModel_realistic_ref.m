proj = GIFproject;
mesh = mesh3D(proj);
mesh.readFile('mesh.msh');
%%
model = GIFmodel(proj, mesh, 1./500);
%%
xyz = mesh.fetchCenter;
active = xyz(:,3)<0.;
sigma = model.value;
sigma(~active) = 1e-8;
model.setValue(sigma);


%%

model.addBlock([-120., -400., -60.], [1000, 400., -200], 1/50.);


model.writeFile('sigma_realistic_ref.con');

