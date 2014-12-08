proj = GIFproject;
mesh = mesh3D(proj);
mesh.readFile('mesh.msh');
%%
model = GIFmodel(proj, mesh, 1./1000);
%%
xyz = mesh.fetchCenter;
active = xyz(:,3)<0.;
sigma = model.value;
sigma(~active) = 1e-8;
layerInd = xyz(:,3)<-60. & xyz(:,3)> -100.;
sigma(layerInd) = 1./50;
layerInd = xyz(:,3)<0. & xyz(:,3)> -60.;
sigma(layerInd) = 1./300;
layerInd = xyz(:,3)<-100. & xyz(:,3)> -200.;
sigma(layerInd) = 1./500;
model.setValue(sigma);
%%


%%

x1 = [0.; -40.; -300.; 0.; 0.];
y1 = [400.; 400.; -400.; -400.; 400.];
z1 = ones(5,1)*-100.;
z2 = ones(5,1)*-200.;
xBlock1 = [x1; x1];
yBlock1 = [y1; y1];
zBlock1 = [z1; z2];
model.addPolyBlock([xBlock1 yBlock1 zBlock1],  1./5);   
model.addBlock([0., -400., -100.], [1000, 400., -200], 1/5.);
model.addBlock([0., 0., -100.], [1000, 1000., -140], 1/50.);

model.writeFile('sigma.con');
% % 1
% % 2
% model.addBlock([-600-100, -600-100, 0. ], [-400-100, -400-100, -200], 1e-3*(1-0.5));
% % 3
% model.addBlock([400+100, 400+100, -50], [600+100, 600+100, -250], 1e-1*(1-0.4));
% % 4
% model.addBlock([400+100, -600-100, -150], [600+100, -400-100, -350], 1e-1*(1-0.8));
% model.writeFile('sigma.con');
