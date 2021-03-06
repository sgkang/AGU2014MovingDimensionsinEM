proj = GIFproject;
mesh = mesh3D(proj);
mesh.readFile('mesh_late.msh');
%%
model = GIFmodel(proj, mesh, 1./1000);
%%
xyz = mesh.fetchCenter;
active = xyz(:,3)<0.;
sigma = model.value;
sigma(~active) = 1e-8;
layerInd = xyz(:,3)<0. & xyz(:,3)> -60.;
sigma(layerInd) = 1./300;
layerInd = xyz(:,3)<-60. & xyz(:,3)> -160.;
sigma(layerInd) = 1./500;
model.setValue(sigma);
%%ttt


%%
x4 = [-60.; 0.; -300.; 0; -60];
x3 = [0.; 0.; -300.; 0; 0];
x2 = [60; 60; -300.+60; 60; 60];
x1 = [120.; 120.; -300.+120; 120; 120];

y1 = [400.; 400.; -400.; -400.; 400.];
z1 = ones(5,1)*-60.;
z2 = ones(5,1)*-100.;
z3 = ones(5,1)*-140.;
z4 = ones(5,1)*-160.;
xBlock1 = [x1; x2; x3; x4]+50;
yBlock1 = [y1; y1; y1; y1];
zBlock1 = [z1; z2; z3; z4];
model.addPolyBlock([xBlock1 yBlock1 zBlock1],  1./5);   
model.addBlock([+60., -400., -60.], [1000, 400., -160], 1/5.);


model.writeFile('sigma_realistic_late.con');
% % 1
% % 2
% model.addBlock([-600-100, -600-100, 0. ], [-400-100, -400-100, -200], 1e-3*(1-0.5));
% % 3
% model.addBlock([400+100, 400+100, -50], [600+100, 600+100, -250], 1e-1*(1-0.4));
% % 4
% model.addBlock([400+100, -600-100, -150], [600+100, -400-100, -350], 1e-1*(1-0.8));
% model.writeFile('sigma.con');
