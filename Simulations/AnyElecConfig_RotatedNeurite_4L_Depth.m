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
Z_Centre = find(Z == 0);

simulation = 17;

SimulationParameters

% Define membrane thresholds
Th = [12.09e-3 6.30e-3*ones(1,length(rot)-1)];

%% Calculate the longitudinal and transverse components separately

Ve_z = zeros(length(rot),length(Z));
Ve_x = zeros(length(rot),length(X));

h = waitbar(0,'Calculating Rotated Membrane Potential');


for i = 1:length(rot)
    
    [Ve_f, ~, ~, ~,~] = CellComp4Layer_Ve_f_Plane(...
        Xi, Yi, Zi, I_M, I_D, ...                   % Electrode parameters
        x_max, z_max, t_max, d_x, d_z, d_t, ...     % Spatial sampling
        h_F, Ya(i), rot(i), Ri ...
        );
    
    % Get Ve
    Ve_all = real(fftshift(ifftn(ifftshift((2*pi)^(3/2)/d_z/d_x/d_t*Ve_f))));
    
    Ve_z(i,:) = squeeze(max(squeeze(Ve_all(X_Centre,:,:)),[],2));
    Ve_x(i,:) = squeeze(max(squeeze(Ve_all(:,Z_Centre,:)),[],2));
    
    waitbar(i/length(rot),h)
    
end

close(h)

%% Normalise layers so we can clearly see the spread

Ve_z_tmp = Ve_z - repmat(min(Ve_z,[],2),1,length(Z));
Ve_x_tmp = Ve_x - repmat(min(Ve_x,[],2),1,length(X));
Ve_z_norm = Ve_z_tmp./repmat(max(Ve_z_tmp,[],2),1,length(Z));
Ve_x_norm = Ve_x_tmp./repmat(max(Ve_x_tmp,[],2),1,length(X));

%% Plot normalised spread with normalised contours

figz = figure;
imagesc(Z*1e6,Ya*1e6,Ve_z_norm), axis xy, xlim([-1000 1000]), ylim([Ya(end) Ya(1)]*1e6)
hold on
contour(Z*1e6,Ya*1e6,Ve_z_norm,[0.5 0.5],'k')

figx = figure;
imagesc(X*1e6,Ya*1e6,Ve_x_norm), axis xy, xlim([-1000 1000]), ylim([Ya(end) Ya(1)]*1e6)
hold on
contour(X*1e6,Ya*1e6,Ve_x_norm,[0.5 0.5],'k')

saveas(figz, ['Figures/' folderName '/' figName 'z'],'fig')
saveas(figx, ['Figures/' folderName '/' figName 'x'],'fig')

%% Plot spread with normalised contours

Ve_z_shift = Ve_z - min(Ve_z(:));
Ve_x_shift = Ve_x - min(Ve_x(:));

figz = figure;
imagesc(Z*1e6,Ya*1e6,Ve_z_shift), axis xy, xlim([-1000 1000]), ylim([Ya(end) Ya(1)]*1e6);
hold on
clim = figz.CurrentAxes.CLim;
contour(Z*1e6,Ya*1e6,Ve_z_norm,[0.5 0.5],'k')
figz.CurrentAxes.CLim = clim;

figx = figure;
imagesc(X*1e6,Ya*1e6,Ve_x_shift), axis xy, xlim([-1000 1000]), ylim([Ya(end) Ya(1)]*1e6)
hold on
clim = figz.CurrentAxes.CLim;
contour(X*1e6,Ya*1e6,Ve_x_norm,[0.5 0.5],'k')
figz.CurrentAxes.CLim = clim;

saveas(figz, ['Figures/' folderName '/' figName 'z_unnorm'],'fig')
saveas(figx, ['Figures/' folderName '/' figName 'x_unnorm'],'fig')
