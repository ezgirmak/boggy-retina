addpath([fileparts(pwd()) '\Utilities'])
addpath(fileparts(pwd()))
addpath('D:\Dropbox\PhD\MATLAB\Trees Package 1.15')

cd('D:\Dropbox\PhD\MATLAB\Data\RGC Cells\badea\CNG version');
% cd('C:\Users\tesler\Dropbox\PhD\MATLAB\Data\RGC Cells\chalupa\CNG version');
files = dir;
start_trees;
d_z = 2e-6;

ColorSet = varycolor(length(files)-2);
P4 = figure('units','normalized','outerposition',[1.5 -0.03 0.64 1.03],'Color',[1 1 1]);
ax = axes;
set(ax, 'ColorOrder', ColorSet);
hold on

for i = 3:length(files)
    load_tree(files(i).name);
    [Xa, Ya, Za, AIS, PL] = getAxonXYZ(trees{i-2},d_z*1e6,'extend');
    
    XYZ = {Xa, Ya, Za};
    XYZ = smoothn(XYZ,50);
    Xa = XYZ{1};
    Ya = XYZ{2};
    Za = XYZ{3};
    Xa(1) = 0;
    Ya(1) = 0;
    Za(1) = 0;
    
    zp_max = PL*1e-6/2;
    kzp_max = pi/d_z;
    nzp = (fix(zp_max/d_z));
    
    %% Plot geometry of stimulation
    
    h1 = plot3(Za(nzp+1:end)*1e6,Xa(nzp+1:end)*1e6,Ya(nzp+1:end)*1e6,'-d','LineWidth',2,'MarkerSize',1);
    title('Axon Geometries')
end

hold off
xlabel('Z')
ylabel('X')
zlabel('Y')
axis equal
grid minor
% zlim([-40 75])
% xlim([-420 600])
view(45,20)