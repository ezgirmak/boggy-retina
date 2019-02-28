function [Ve_T,phi_max] = CellComp_VeT_Plane(Xi, Yi, Zi, I_M, I_D, x_max, z_max, t_max, d_x, d_z, d_t, phi)
%%CELLCOMP_VET_PLANE
% CELLCOMP_VET_PLANE returns the transverse component of the
% electric field for a flat plane in tissue comprised of parallel
% fibers given a set of stimulating electrodes, modelled as point sources.
%
% The membrane voltage is calculated using a modified version of the
% self-consistent, linear, sub-threshold model presented in:
%
%   B. Tahayori, H. Meffin, E.N. Sergeev, I.M.Y. Mareels, A.N. Burkitt, and
%   D.N. Grayden (2014), "Modelling extracellular electrical stimulation:
%   IV. Effect of the cellular composition of neural tissue on its
%   spatio-temporal filtering properties", J. Neural Eng. 11.
%
% INPUTS:
%
% Xi, Yi and Zi         the x-, y-, and z- coordinates of the point source
%                       electrodes.
% I_M and I_D           the stimulation amplitude and duration for each
%                       electrode.
% z_max, x_max, t_max   z-, x- and time- extent of the simulation.
% d_z, d_x, d_t         z- and x-direction and time step sizes.
% rot                   coordinate rotation (use if the output will be used
%                       to calculate the membrane potential for an
%                       anomalous rotated fibre.
%
% OUTPUTS:
%
% Ve_T                  The calculated extracellular potential along the
%                       z-direction of the neurite for a full plane.
%
%
% EXAMPLE USAGE:
%
% Xi=0; Yi=30e-6; Zi=0;
% I_M=-1e-6; I_D=100e-6;
% z_max=300e-6; x_max=300e-6; t_max=600e-6;
% d_z=4e-6; d_x=4e-6; d_t=2e-6;
% 
% [Ve_L] = CellComp_VeT_Plane(Xi,Yi,Zi,I_M,I_D,z_max,x_max,t_max,d_z,d_x,d_t,phi);
%
% Created by: Tim Esler, 2015

% addpath('Utilities')

%% Define parameters

p = NTESparams;         % Call parameter function

% Unpack parameters
b = p.b;                % NTES radius (m)
d = p.d;                % Width of extracellular sheath (m)
a = p.a;                % Neurite radius (m)

C_m = p.C_m;            % Membrane capacitance (F/m^2)
R_m = p.R_m;            % Membrane unit area resistance (ohm.m^2)

rho_i = p.rho_i;        % Intracellular resistivity (ohm.m)
rho_e = p.rho_e;        % Extracellular resistivity (ohm.m)
r_m = p.r_m;            % Membrane unit length resistance (ohm.m)
r_i = p.r_i;            % Intracellular resistance (ohm/m)
r_e = p.r_e;            % Extracellular resistance (ohm/m)

%% Define sampling in time-space and Fourier domains

% Sampling space Fourier domains
kz_max = pi/d_z;
nz = (fix(z_max/d_z));
d_kz = kz_max/nz;
kz = single(-kz_max:d_kz:kz_max);
kz(fix(length(kz)/2+1)) = 1e-20;

kx_max = pi/d_x;
nx = (fix(x_max/d_x));
d_kx = kx_max/nx;
kx = single(-kx_max:d_kx:kx_max);
kx(fix(length(kx)/2+1)) = 1e-20;

% Sampling time Fourier domain
w_max = pi/d_t;
nt = (fix(t_max/d_t));
d_w = w_max/nt;
w = single(-w_max:d_w:w_max);
w(fix(length(w)/2+1)) = 1e-20;

% Sampling space domain
Z = single(-z_max:d_z:z_max);
Z(fix(length(kz)/2+1)) = 1e-20;

X = single(-x_max:d_x:x_max);
X(fix(length(kx)/2+1)) = 1e-20;

% Sampling time domain
T = single(-t_max:d_t:t_max);
T(fix(length(w)/2+1)) = 1e-20;

% Create sample mesh
[kx_m,kz_m,w_m] = ndgrid(kx,kz,w);

%% Define electrotonic length constants, time constants and admittivities
% in the Fourier domain

tau_m = R_m*C_m;                % Membrane time constant (s)

% Electrotonic length constants (static and frequency-dependent, for both
% current density (J) and voltage (V) boundary conditions)
L_0J = sqrt(r_m/(r_e+r_i));
L_0V = sqrt(r_m/r_i);
L_J_m = L_0J./sqrt(1+1i*w_m*tau_m);
L_V_m = L_0V./sqrt(1+1i*w_m*tau_m);

% Longitudinal and transverse admittivities in the Fourier domain
xi_L_f_m = 1/rho_i * (1+(kz_m.^2).*L_J_m.^2)./(1+(kz_m.^2).*L_V_m.^2);
xi_T_f = d/b/rho_e;

% Anisotropy ratio
chi_f_m = sqrt(xi_L_f_m/xi_T_f);

%% Iterate through the point sources

Ve_dx = zeros(length(X),length(Z),length(T));
Ve_dy = zeros(size(Ve_dx));

for i = 1:length(I_M)
    %% Define point source stimulation in the time Fourier domain
    
    % Monophasic
    % I_hat = I_M(i)/sqrt(2*pi)*I_D(i)*exp(-1i*I_D(i)*w_m/2) ...
    %     .*sinc(I_D(i)*w_m/2/pi);
    
    % Biphasic
    I_hat = 1i*2*I_D(i)*I_M(i)/sqrt(2*pi)*exp(-1i*I_D(i)*w_m) ...
        .*sinc(I_D(i)*w_m/2/pi).*sin(I_D(i)*w_m/2);
    
    %% Calculate extracellular voltage contributed by this source
    % in the Fourier domain
    
    % Define alpha
    aa = sqrt(kx_m.^2 + chi_f_m.^2.*kz_m.^2);
    
    % Extracellular voltage in Fourier domain. NOTE: since the rotated
    % fibre lies in the plane defined by y = 0, we have y-Yi(i) = -Yi(i).
    Ve_dx_f = 1i*kx_m.*I_hat.*exp(-1i*(Xi(i)*kx_m+Zi(i)*kz_m)).* ...
        exp(-aa*abs(Yi(i)))./(4*pi*xi_T_f*aa);
    Ve_dy_f = -sign(-Yi(i))*I_hat.*exp(-1i*(Xi(i)*kx_m+Zi(i)*kz_m)).* ...
        exp(-aa*abs(Yi(i)))./(4*pi*xi_T_f);
    
%     % Radial unit vector in terms of rotated coordinates
%     vec_rp = [-Xip(i)+X; ...
%             -ones(1,length(X))*Yi(i); ...
%             zeros(1,length(X))];
%     u_rp = vec_rp./repmat(sqrt(vec_rp(1,:).^2 + vec_rp(2,:).^2 + vec_rp(3,:).^2),3,1);
%     u_rp_xp = sign(I_M(i))*repmat(u_rp(1,:),length(Z),1,length(T));
%     u_rp_y = repmat(u_rp(2,:),length(Z),1,length(T));
    
    Ve_dx = Ve_dx + b/2*real(fftshift(ifftn(ifftshift((2*pi)^(3/2)/d_z/d_x/d_t*Ve_dx_f))));
    Ve_dy = Ve_dy + b/2*real(fftshift(ifftn(ifftshift((2*pi)^(3/2)/d_z/d_x/d_t*Ve_dy_f))));
end

%% Calculate total extracellular voltage in the Fourier domain

Ve_T_mag = sqrt(Ve_dx.^2 + Ve_dy.^2);

Ve_T_Centre = squeeze(Ve_T_mag(X == single(1e-20),:,:));
[~,Ve_T_maxInd] = max(Ve_T_Centre(:));
[Ve_T_Zind,Ve_T_Tind] = ind2sub(size(Ve_T_Centre),Ve_T_maxInd);

phi_max = atan2(Ve_dy(Ve_T_Zind,Ve_T_Tind),-Ve_dx(Ve_T_Zind,Ve_T_Tind));

if strcmpi(phi, 'phi_max')
    Ve_T = (Ve_dx*cos(phi_max) + Ve_dy*sin(phi_max));
elseif isnumeric(phi)
    Ve_T = (Ve_dx*cos(phi) + Ve_dy*sin(phi));
elseif strcmpi(phi, 'magnitude')
    Ve_T = Ve_T_mag;
else
    error(['Argument ''phi'' must be either ',...
        '''phi_max'', ',...
        '''magnitude'' or a numeric constant'])
end

end
