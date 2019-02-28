function [Ve_T,phi_max] = CellComp_VeT(Xi, Yi, Zi, I_M, I_D, phi, z_max, t_max, d_z, d_t)
%%CELLCOMP_VET
% CELLCOMP_VET returns the transverse component of the
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
% Ve_T                  The calculated extracellular potential along the
%                       z-direction of the neurite.
%
% Created by: Tim Esler, 2015

% addpath('Utilities')

if range(I_D) ~= 0
    error('Code cannot handle multiple pulse durations yet...')
end

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
Thetai = atan2(Yi,Xi);

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
% [w_m,kz_m] = meshgrid(w,kz);
kz_m = kz;

Ve_T_sources = zeros(length(I_M),length(Z),length(T),'single');
Ve_T_sourcesX = zeros(size(Ve_T_sources),'single');
Ve_T_sourcesY = zeros(size(Ve_T_sources),'single');

for j = 1:length(w)
    w_m = w(j);

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



    for i = 1:length(I_M)
        %% Define point source stimulation in the time Fourier domain

        % Monophasic
        % I_hat = I_M/sqrt(2*pi)*I_D*exp(-1i*I_D*w_m/2) ...
        %     .*sinc(I_D*w_m/2/pi);

        % Biphasic
        I_hat = 1i*2*I_D(i)*I_M(i)/sqrt(2*pi)*exp(-1i*I_D(i)*w_m) ...
            .*sinc(I_D(i)*w_m/2/pi).*sin(I_D(i)*w_m/2);

        %% Calculate extracellular voltage contributed by this source
        % Contributions for each source are calculated as magnitudes only,
        % neglecting the rotational variation around the surface of the
        % neurite. The corresponding vectors from each source are then summed
        % and the angle of the resultant vector determines the variation around
        % the surface of the neurite via cos(phi - ThetaT).

        % Extracellular voltage in Fourier domain
        Ve_T_f = I_hat.*exp(-1i*kz_m*Zi(i))*b.*chi_f_m.*abs(kz_m) ...
            ./(xi_T_f*2*(2*pi)^(3/2)).*besselk(1,Ri(i)*chi_f_m.*abs(kz_m));

        %% Calculate extracellular voltage in the time-space domain
        Ve_T_sources(i,:,j) = (2*pi)/d_z/d_t*real(fftshift(ifft(ifftshift(Ve_T_f))));

        % Get X and Y components of the corresponding transverse vector
        Ve_T_sourcesX(i,:,j) = Ve_T_sources(i,:,j)*cos(Thetai(i));
        Ve_T_sourcesY(i,:,j) = Ve_T_sources(i,:,j)*sin(Thetai(i));
    end
end
%% Calculate total extracellular potential

Ve_T_X = squeeze(sum(Ve_T_sourcesX,1))+1e-10;
Ve_T_Y = squeeze(sum(Ve_T_sourcesY,1))+1e-10;
[~,Ve_T_maxInd] = max(Ve_T_X(:).^2 + Ve_T_Y(:).^2);
[Ve_T_Zind,Ve_T_Tind] = ind2sub(size(Ve_T_X),Ve_T_maxInd);

phi_max = atan2(sin(-Thetai+pi)*Ve_T_sources(:,Ve_T_Zind,Ve_T_Tind), ...
    -cos(-Thetai+pi)*Ve_T_sources(:,Ve_T_Zind,Ve_T_Tind));

if strcmpi(phi, 'phi_max')
    phi_scale = repmat(cos(phi_max-Thetai)',1,length(Z),length(T));
    
    % Output extracellular voltage at phi
    Ve_T = squeeze(sum(Ve_T_sources.*-phi_scale,1));
elseif isnumeric(phi)
    phi_scale = repmat(cos(phi-Thetai)',1,length(Z),length(T));
    
    % Output extracellular voltage at phi
    Ve_T = squeeze(sum(Ve_T_sources.*-phi_scale,1));
elseif strcmpi(phi, 'magnitude')
    % Output extracellular voltage at phi
    Ve_T = sqrt(Ve_T_X.^2 + Ve_T_Y.^2);
else
    error(['Argument ''phi'' must be one of',...
        '''phi_max''',...
        '''max_phi'' or a numeric constant'])
end

end