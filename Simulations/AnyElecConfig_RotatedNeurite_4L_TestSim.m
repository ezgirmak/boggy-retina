addpath([fileparts(pwd()) '\Utilities'])
addpath([fileparts(pwd()) '\Utilities\altmany-export_fig-04ca93c'])
addpath(fileparts(pwd()))

%% Parameters
    
% Define simulation size and step sizes
x_max = 2000e-6;
z_max = 2000e-6;
t_max = 600e-6;
d_x = 10e-6;
d_z = 10e-6;
d_t = 4e-6;

Z = -z_max:d_z:z_max;
T = -t_max:d_t:t_max;
X = -x_max:d_x:x_max;
X_Centre = find(X == 0);

%% One Electrode Simulation
    
% Define neurite rotation in xz plane
rot = 0;

% Define source locations
Xi = 0e-6;
Yi = 200e-6;
Zi = 0e-6;

% Define electrode radius
Ri = 50e-6;

% Define axon locations
h_F = 100e-6;
Ya = -10e-6;

% Define source amplitudes and pulse durations
I_M = -50e-6;
I_D = 200e-6;

% Define location on neurite to calculate the membrane potential
phi = 'magnitude';

%% Calculate the longitudinal and transverse components separately
    
[Ve_L_f, Ve_T_X_f, Ve_T_Y_f, ~,n_e] = CellComp4Layer_Ve_f_Plane(...
    Xi, Yi, Zi, I_M, I_D, ...                   % Electrode parameters
    x_max, z_max, t_max, d_x, d_z, d_t, ...     % Spatial sampling
    h_F, Ya, rot, Ri ...
    );

% Get Vm_L
Vm_L_f = n_e.*Ve_L_f;
Vm_L = real(fftshift(ifftn(ifftshift((2*pi)^(3/2)/d_z/d_x/d_t*Vm_L_f))));
% Get Vm_T
Ve_T_X = real(fftshift(ifftn(ifftshift((2*pi)^(3/2)/d_z/d_x/d_t*Ve_T_X_f))));
Ve_T_Y = real(fftshift(ifftn(ifftshift((2*pi)^(3/2)/d_z/d_x/d_t*Ve_T_Y_f))));

Vm_L_rot = squeeze(Vm_L(X_Centre,:,:));
Vm_T_rot = 2*sqrt(squeeze(Ve_T_X(X_Centre,:,:)).^2 + ...
    squeeze(Ve_T_Y(X_Centre,:,:)).^2);

%% Sum the longitudinal and tranverse membrane potentials

Vm_rot = Vm_L_rot + Vm_T_rot;

%% Plot

figSize = [5 2 8 6]*3;
fontSize = 6*3;

fig = figure('Units','centimeters','Position',figSize);
set(fig, 'PaperPosition',figSize)

ax = subplot(1,1,1);
hold on
plot(Z*1e6,squeeze(Vm_rot(:,floor((t_max+I_D(1))/d_t)))*1e3,'LineWidth',1.5)
xlabel('Neurite axis (\mum)')
ylabel('V_{m} (mV)')
xlim([-1000 1000])
ylim([-6 14])
set(ax,'XTick',[-1000 0 1000])
set(ax,'YTick',[-4 0 4 8 12])
set(gca,'fontsize',fontSize)
ax.Box = 'off';

outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

%% Show SketchUp exported image

% winopen(['Figures\' filename])