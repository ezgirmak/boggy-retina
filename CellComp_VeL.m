function [Ve_L] = CellComp_VeL(Xi, Yi, Zi, I_M, I_D, z_max, t_max, d_z, d_t)
%%CELLCOMP_VEL
% CELLCOMP_VEL returns the longitudinal component of the
% extracellular voltage for a neurite in tissue comprised of parallel
% fibers given a set of stimulating electrodes, modelled as point sources.
% The extracellular voltage is calculated using a modified version of the
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
% z_max and t_max       the z- and time- extent of the simulation.
% d_z and d_t           the z-direction and time step sizes.
%
% OUTPUTS:
%
% Ve_L                  The calculated extracellular potential along the
%                       z-direction of the neurite.
%
% Created by: Tim Esler, 2015

% addpath('Utilities')

%% Define parameters

p = NTESparams('single');

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

% Convert Xi, Yi source coordinates to Ri and Thetai values
Ri = sqrt(Xi.^2 + Yi.^2);
Thetai = atan(Yi./Xi);

%% Define sampling in time-space and Fourier domains

% Sampling space Fourier domain
kz_max = pi/d_z;
nz = (fix(z_max/d_z));
d_kz = kz_max/nz;
kz = single(-kz_max:d_kz:kz_max);
kz(fix(length(kz)/2+1)) = 1e-20;

% Sampling time Fourier domain
w_max = pi/d_t;
nt = (fix(t_max/d_t));
d_w = w_max/nt;
w = single(-w_max:d_w:w_max);
w(fix(length(w)/2+1)) = 1e-20;

% Sampling space domain
Z = single(-z_max:d_z:z_max);
Z(fix(length(kz)/2+1)) = 1e-20;

% Sampling time domain
T = single(-t_max:d_t:t_max);
T(fix(length(w)/2+1)) = 1e-20;

% Create sample mesh
[w_m,kz_m] = meshgrid(w,kz);

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

Ve_L_f_sources = zeros(length(I_M),length(Z),length(T),'single');

for i = 1:length(I_M)
    %% Define point source stimulation in the time Fourier domain
    
    % Monophasic
    % I_hat = I_M/sqrt(2*pi)*I_D*exp(-1i*I_D*w_m/2) ...
    %     .*sinc(I_D*w_m/2/pi);
    
    % Biphasic
    I_hat = 1i*2*I_D(i)*I_M(i)/sqrt(2*pi)*exp(-1i*I_D(i)*w_m) ...
        .*sinc(I_D(i)*w_m/2/pi).*sin(I_D(i)*w_m/2);
    
    %% Calculate extracellular voltage contributed by this source
    
    % Extracellular voltage in Fourier domain
    Ve_L_f_sources(i,:,:) = I_hat.*exp(-1i*kz_m*Zi(i))./(xi_T_f*(2*pi)^(3/2)) ...
        .*besselk(0,Ri(i)*chi_f_m.*abs(kz_m));
end

%% Calculate total extracellular voltage in the Fourier domain

if length(I_M) > 1
    Ve_L_f = squeeze(sum(Ve_L_f_sources));
else
    Ve_L_f = squeeze(Ve_L_f_sources);
end

%% Calculate extracellular voltage in the time-space domain

Ve_L = (2*pi)/d_z/d_t*real(fftshift(ifftn(ifftshift(Ve_L_f))));

end