function PlotCurrentSpread(Ya, rot)
% Basic setup to plot current spread at  a given axon location and rotation

% if rotation is not given
if nargin<2
    rot=0;
end
% if axon location is not given
if nargin<1
    Ya=-0.015;
end


%%

% Define simulation size and step sizes
x_max = 3000e-6;
z_max = 3000e-6;
t_max = 600e-6;
d_x = 10e-5;
d_z = 10e-5;
d_t = 5e-5;

Z = -z_max:d_z:z_max;
T = -t_max:d_t:t_max;
X = -x_max:d_x:x_max;
X_Centre = find(X == 0);
Z_Centre = find(Z == 0);



% Define source locations
Xi = [0e-6];
Yi = [100e-6];
Zi = [0e-6];

% Define electrode radius
Ri = 50e-6;


% Define source amplitudes and pulse durations
I_M = [-50e-6];
I_D = [200e-6];

% Define axon locations
h_F = 100e-6;
% Define membrane thresholds
%these thresholds were simulated
Th = [12.09e-3 6.30e-3*ones(1,length(rot)-1)];

[Ve_f, ~, ~, ~,~] = CellComp4Layer_Ve_f_Plane(...
    Xi, Yi, Zi, I_M, I_D, ...                   % Electrode parameters
    x_max, z_max, t_max, d_x, d_z, d_t, ...     % Spatial sampling
    h_F, Ya, rot, Ri ...
    );

Ve_all = real(fftshift(ifftn(ifftshift((2*pi)^(3/2)/d_z/d_x/d_t*Ve_f))));

% cathodic is negative current, gives negative potential
% normalize to [0, 1]:
minVe = min(Ve_all(:));
maxVe = max(Ve_all(:));
Ve_norm = (Ve_all - minVe) / (maxVe - minVe);

% however, we think of cathodic (negative) as bright, so flip the color:
Ve_norm = 1 - Ve_norm;

%%
%looking at x-z across time
for i =1:size(Ve_norm,3)
    imagesc(Ve_norm(:,:,i));
    title(['Frame t=' num2str(T(i))]);
    drawnow; pause(1);
end

end