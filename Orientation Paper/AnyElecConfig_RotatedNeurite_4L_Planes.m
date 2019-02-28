addpath([fileparts(pwd()) '\Utilities'])
addpath([fileparts(pwd()) '\Utilities\altmany-export_fig-04ca93c'])
addpath(fileparts(pwd()))

%% Parameters

% Define simulation size and step sizes
x_max = 4000e-6;
z_max = 4000e-6;
t_max = 2000e-6;
d_x = 20e-6;
d_z = 20e-6;
d_t = 10e-6;

Z = -z_max:d_z:z_max;
T = -t_max:d_t:t_max;
X = -x_max:d_x:x_max;

% simulation = 7;

SimulationParameters

% Define membrane thresholds
Th = [12.09e-3 6.30e-3*ones(1,length(rot)-1)];

%% Error checking

if max(rot) ~= pi
    error('"rot" should vary between 0 and pi')
end

%% Calculate the longitudinal and transverse components separately

h = waitbar(0,'Calculating Rotated Membrane Potential');

Vm_max = zeros(length(rot),length(X),length(Z));
iterations = length(rot);
disp(['Iterations: ' num2str(iterations)])

for k = 1:length(rot)
    
    [Ve_L_f, Ve_T_X_f, Ve_T_Y_f, ~,n_e] = CellComp4Layer_Ve_f_Plane(...
        Xi, Yi, Zi, I_M, I_D, ...                   % Electrode parameters
        x_max, z_max, t_max, d_x, d_z, d_t, ...     % Spatial sampling
        h_F, Ya(k), rot(k), Ri ...
        );
    
    % Get Vm_L
    Vm_L_f = n_e.*Ve_L_f;
    clear Ve_L_f
    Vm_L = real(fftshift(ifftn(ifftshift((2*pi)^(3/2)/d_z/d_x/d_t*Vm_L_f))));
    clear Vm_L_f
    % Get Vm_T
    Ve_T_X = real(fftshift(ifftn(ifftshift((2*pi)^(3/2)/d_z/d_x/d_t*Ve_T_X_f))));
    clear Ve_T_X_f
    Ve_T_Y = real(fftshift(ifftn(ifftshift((2*pi)^(3/2)/d_z/d_x/d_t*Ve_T_Y_f))));
    clear Ve_T_Y_f
    
    Vm_L_rot = squeeze(Vm_L);
    Vm_T_rot = 2*squeeze(sqrt(Ve_T_X.^2 + Ve_T_Y.^2));
    clear Vm_L Ve_T_X Ve_T_Y
    
    %% Sum the longitudinal and tranverse membrane potentials
    
    Vm_rot = Vm_L_rot + Vm_T_rot;
    
    Vm_max(k,:,:) = squeeze(max(abs(Vm_rot),[],3));
    
    waitbar(k/length(rot),h)
    
end
close(h)

%% Pull out AOP and AIS components and scale

Vm_max_AOP = Vm_max(1,:,:);

Vm_max = Vm_max(2:end,:,:);

rot_AOP = rot(1);
Ya_AOP = Ya(1);
rot = rot(2:end);
Ya = Ya(2:end);

%% Upsample w.r.t. orientation values

% Interestingly, Vm_max vs rot is predicted perfectly by a cubic
% polynomial, so we use this to upsample
% TO-DO: it is actually fit perfectly by a sum of sines and cosines, so use
% that
rot_upsamp = 0:d_theta_upsamp:max(rot);

Vm_max_upsamp = interp1(rot,Vm_max,rot_upsamp,'pchip');

clear Vm_max

% Scale the results by the amount of extra current required to cause
% activation in the NFL
scale = 1;
if scaleOut && islogical(scaleOut)
    scale = Th(1)/max(abs(Vm_max_AOP(:)));
elseif ~islogical(scaleOut)
    Vm_maxPerRot = sort(squeeze(max(max(abs(Vm_max_upsamp),[],3),[],2)));
	indScale = round(length(rot_upsamp)*(1-scaleOut));
	scale = Th(2)/Vm_maxPerRot(indScale);
end
Vm_max_scaled = Vm_max_upsamp*scale;

%% Get mean preferential activation

meanPrefAct = squeeze(mean(Vm_max_scaled > Th(2),1));

%% Plot Mean prefential activation

fig = figure('Units','centimeters','Position',[5 3 24 16]);
set(fig, 'PaperPosition',[5 3 24 16])

ax = subplot(1,1,1);
imagesc(X*1e6, Z*1e6, meanPrefAct,[0 1]);
axis xy
colormap(hot)
cb = colorbar;
xlabel('X (\mum)')
ylabel('Y (\mum)')
ylabel(cb,'Mean prefential activation (mV)')
xlim([-1500 1500])
ylim([-1500 1500])

hold on
contour(X*1e6, Z*1e6, meanPrefAct, [0.001 0.001], ...
    'color', 'w', 'LineWidth', 1.2,'LineStyle','--');

set(ax,'fontsize',16)
set(cb,'fontsize',16)

saveas(fig, ['Figures/' folderName '/' figName],'fig')
saveas(fig, ['Figures/' folderName '/' figName],'jpeg')
saveas(fig, ['Figures/' folderName '/' figName],'eps')

