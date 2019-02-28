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

% simulation = 4;

SimulationParameters

%% Calculate the longitudinal and transverse components separately

h = waitbar(0,'Calculating Rotated Membrane Potential');
count = 0;
iterations = length(Yi_V)*length(I_D_V)*length(rot);
disp(iterations)

Vm_max = zeros(length(rot),length(Yi_V),length(I_D_V));

for i = 1:length(Yi_V)
    
    Yi = ones(1,length(I_M))*Yi_V(i);
    
    for j = 1:length(I_D_V)
        
        I_D = ones(1,length(I_M))*I_D_V(j);
        
        for k = 1:length(rot)
            tic
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
            
            Vm_max(k,i,j) = max(abs(Vm_rot(:)));
            
            count = count + 1;
            waitbar(count/iterations,h)
            
            toc
        end
    end
end
close(h)

% Scale the results by the amount of extra current required to cause
% activation in the NFL
Vm_max_scaled = Vm_max(2:end,:,:).*repmat(Th(1)./Vm_max(1,:,:),length(rot)-1,1);

% Interestingly, Vm_max vs rot is predicted perfectly by a cubic
% polynomial, so we use this to upsample
% TO-DO: it is implemented perfectly by a sum of sines and cosines, so use
% that
rot_upsamp = 0:d_theta_upsamp:pi/2;
Vm_max_scaled_upsamp = zeros(length(rot_upsamp),length(Yi_V),length(I_D_V));
for i = 1:length(Yi_V)
    for j = 1:length(I_D_V)
        cubicFit = polyfit(rot(2:end),Vm_max_scaled(:,i,j)',3);
        Vm_max_scaled_upsamp(:,i,j) = polyval(cubicFit, rot_upsamp);
    end
end

Vm_threshold = double(Vm_max_scaled_upsamp >= Th(2));

prefAct = squeeze(sum(Vm_threshold, 1))./length(rot_upsamp);

% Interpolate points for smoothness (fine since everything should vary
% smoothly)
Yi_upsamp = (Yi_V(1):d_Yi_upsamp:Yi_V(end))';
I_D_upsamp = I_D_V(1):d_I_D_upsamp:I_D_V(end);
prefAct_smooth = interp2(I_D_V, Yi_V, prefAct, I_D_upsamp, Yi_upsamp);

%% Plot for paper
fig = figure('Units','centimeters','Position',[5 -5 20 17]);
set(fig, 'PaperPosition',[5 -5 20 17])

ax = subplot(1,1,1);
imagesc(Yi_upsamp*1e6, I_D_upsamp*1e6, prefAct_smooth, [0 1]);
axis xy
colormap(hot)
cb = colorbar;
xlabel('Pulse Duration (\mus)')
ylabel('Electrode Height (\mum)')
ylabel(cb,'AIS proportion preferentially activated')

set(ax,'YTick',0:100:max(Yi_V)*1e6, ...
    'Color','w', ...
    'XTick',0:100:max(I_D_V)*1e6)
set(cb,'YTick',0:0.1:1)

hold on
[C, H] = contour(I_D_upsamp*1e6,Yi_upsamp*1e6,prefAct_smooth,contLevels, ...
    'color', 'w', 'LineWidth', 1.2,'LineStyle','--');
% text(I_D_upsamp(end)*1e6, ...
%     Yi_upsamp(find(prefAct_smooth(:,end)>contLevels(1),1)-1)*1e6,'Low', ...
%     'HorizontalAlignment','right','VerticalAlignment','top', ...
%     'FontSize',16,'Color','w')
% text(I_D_upsamp(end)*1e6, ...
%     Yi_upsamp(find(prefAct_smooth(:,end)>contLevels(2),1)-1)*1e6,'Medium', ...
%     'HorizontalAlignment','right','VerticalAlignment','top', ...
%     'FontSize',16,'Color','w')
% text(I_D_upsamp(end)*1e6, ...
%     Yi_upsamp(end)*1e6,'High', ...
%     'HorizontalAlignment','right','VerticalAlignment','top', ...
%     'FontSize',16,'Color','w')

set(ax,'fontsize',16)
set(cb,'fontsize',16)

saveas(fig, ['Figures/' folderName '/' figName],'fig')
saveas(fig, ['Figures/' folderName '/' figName],'jpeg')
saveas(fig, ['Figures/' folderName '/' figName],'eps')

