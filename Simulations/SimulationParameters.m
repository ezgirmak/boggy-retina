if simulation == -1
    
    % Calculate distribution of orientation from NeuroMorpho
    
    % Define AIS and AOP location
    AISLen = 40;
    AISPos = 40;
    AOPLen = 20;
    
    % Define distance from soma at which AOP starts
    AOPPos = [100 300 500];
    
    % Histogram bins
    bins = 6;
    
    % Folder name
    folderName = 'Figure - Orientation Distribution';
    
    % Plot filename
    figName = 'OrientationDistribution';
    
elseif simulation == 0
    
    %% Generate SketchUp Model
    
    % Define neurite rotations in xz plane
    d_rot = pi/18;
    rot = [0 0:d_rot:pi/2+pi/10000];
    
    % Define source locations
    Xi = [-200e-6 -200e-6 -200e-6 -200e-6 ...
        0e-6 0e-6 0e-6 0e-6 ...
        200e-6 200e-6 200e-6 200e-6];
    Yi = [100e-6 100e-6 100e-6 100e-6 ...
        100e-6 100e-6 100e-6 100e-6 ...
        100e-6 100e-6 100e-6 100e-6];
    Zi = [-300e-6 -100e-6 100e-6 300e-6 ...
        -300e-6 -100e-6 100e-6 300e-6 ...
        -300e-6 -100e-6 100e-6 300e-6];
    
    % Define electrode radius
    Ri = [50e-6 50e-6 50e-6 50e-6 ...
        50e-6 50e-6 50e-6 50e-6 ...
        50e-6 50e-6 50e-6 50e-6];
    
    % Define axon locations
    h_F = 100e-6;
    Ya = [-10e-6 -110e-6*ones(1,length(rot)-1)];
    
    % Define source amplitudes and pulse durations
    I_M = [-50e-6 -50e-6 -50e-6 -50e-6 ...
        -50e-6 -50e-6 -50e-6 -50e-6 ...
        -50e-6 -50e-6 -50e-6 -50e-6];
    I_D = [100e-6 100e-6 100e-6 100e-6 ...
        100e-6 100e-6 100e-6 100e-6 ...
        100e-6 100e-6 100e-6 100e-6];
    
    % Define location on neurite to calculate the membrane potential
    phi = 'magnitude';
    
    % Should we scale the output so the AOP just hits threshold?
    scaleOut = true;
    
    % Folder name
    folderName = 'Figure - Single Voltage Traces';
    
    % Sketchup filename
    filename = 'SimGeom_AllElectrodes.jpg';
    
elseif simulation == 1
    
    %% One Electrode Simulation
    
    % Define neurite rotations in xz plane
    d_rot = pi/18;
    rot = [0 0:d_rot:pi/2+pi/10000];
    
    % Define source locations
    Xi = [0e-6];
    Yi = [200e-6];
    Zi = [0e-6];
    
    % Define electrode radius
    Ri = 50e-6;
    
    % Define axon locations
    h_F = 100e-6;
    Ya = [-10e-6 -110e-6*ones(1,length(rot)-1)];
    
    % Define source amplitudes and pulse durations
    I_M = [-50e-6];
    I_D = [200e-6];
    
    % Define location on neurite to calculate the membrane potential
    phi = 'magnitude';
    
    % Should we scale the output so the AOP just hits threshold?
    scaleOut = true;
    
    % Folder name
    folderName = 'Figure - Single Voltage Traces';
    
    % Sketchup filename
    filename = 'SimGeom_OneElectrode.jpg';
    
    % Plot filename
    figName = 'RotNeurite_OneElectrode';
    
elseif simulation == 2
    
    %% Two Electrode Simulation
    
    % Define neurite rotations in xz plane
    d_rot = pi/18;
    rot = [0 0:d_rot:pi/2+pi/10000];
    
    % Define source locations
    Xi = [0e-6 0e-6];
    Yi = [200e-6 200e-6];
    Zi = [-100e-6 100e-6];
    
    % Define electrode radius
    Ri = [50e-6 50e-6];
    
    % Define axon locations
    h_F = 100e-6;
    Ya = [-10e-6 -110e-6*ones(1,length(rot)-1)];
    
    % Define source amplitudes and pulse durations
    I_M = [-50e-6 -50e-6];
    I_D = [200e-6 200e-6];
    
    % Define location on neurite to calculate the membrane potential
    phi = 'magnitude';
    
    % Should we scale the output so the AOP just hits threshold?
    scaleOut = true;
    
    % Folder name
    folderName = 'Figure - Single Voltage Traces';
    
    % Sketchup filename
    filename = 'SimGeom_TwoElectrode.jpg';
    
    % Plot filename
    figName = 'RotNeurite_TwoElectrode';
    
elseif simulation == 3
    
    %% Four Electrode Simulation
    
    % Define neurite rotations in xz plane
    d_rot = pi/18;
    rot = [0 0:d_rot:pi/2+pi/10000];
    
    % Define source locations
    Xi = [0e-6 0e-6 0e-6 0e-6];
    Yi = [200e-6 200e-6 200e-6 200e-6];
    Zi = [-300e-6 -100e-6 100e-6 300e-6];
    
    % Define electrode radius
    Ri = [50e-6 50e-6 50e-6 50e-6];
    
    % Define axon locations
    h_F = 100e-6;
    Ya = [-10e-6 -110e-6*ones(1,length(rot)-1)];
    
    % Define source amplitudes and pulse durations
    I_M = [-50e-6 -50e-6 -50e-6 -50e-6];
    I_D = [200e-6 200e-6 200e-6 200e-6];
    
    % Define location on neurite to calculate the membrane potential
    phi = 'magnitude';
    
    % Should we scale the output so the AOP just hits threshold?
    scaleOut = true;
    
    % Folder name
    folderName = 'Figure - Single Voltage Traces';
    
    % Sketchup filename
    filename = 'SimGeom_FourElectrode.jpg';
    
    % Plot filename
    figName = 'RotNeurite_FourElectrode';
    
elseif simulation == 4
    
    %% One electrode parameter search
    
    % Define neurite rotations in xz plane
    d_rot = pi/6;
    rot = [0 0:d_rot:pi/2+pi/10000];
    
    % Define source locations
    Xi = [0e-6];
    Yi_V = [10 (20:40:500)]*1e-6;
    Zi = [0e-6];
    
    % Define electrode radius
    Ri = [50e-6];
    
    % Define axon locations
    h_F = 100e-6;
    Ya = [-10e-6 -110e-6*ones(1,length(rot)-1)];
    
    % Define source amplitudes and pulse durations
    I_M = [-1e-6];
    I_D_V = [10 (20:40:500)]*1e-6;
    
    % Define location on neurite to calculate the membrane potential
    phi = 'magnitude';
    
    % Define post-simulation resampling
    d_theta_upsamp = pi/90;
    d_Yi_upsamp = 10e-6;
    d_I_D_upsamp = 10e-6;
    
    % Indicator contours
    contLevels = [0.1 0.4];
    
    % Folder name
    folderName = 'Figure - Param Search';
    
    figName = 'ParamSearch_OneElectrode';
    
    % Define membrane thresholds
    Th = [12.09e-3 6.30e-3*ones(1,length(rot)-1)];
    
elseif simulation == 5
    
    %% Two electrode parameter search
    
    % Define neurite rotations in xz plane
    d_rot = pi/6;
    rot = [0 0:d_rot:pi/2+pi/10000];
    
    % Define source locations
    Xi = [0e-6 0e-6];
    Yi_V = [10 (20:40:500)]*1e-6;
    Zi = [-100e-6 100e-6];
    
    % Define electrode radius
    Ri = [50e-6 50e-6];
    
    % Define axon locations
    h_F = 100e-6;
    Ya = [-10e-6 -110e-6*ones(1,length(rot)-1)];
    
    % Define source amplitudes and pulse durations
    I_M = [-1e-6 -1e-6];
    I_D_V = [10 (20:40:500)]*1e-6;
    
    % Define location on neurite to calculate the membrane potential
    phi = 'magnitude';
    
    % Define post-simulation resampling
    d_theta_upsamp = pi/90;
    d_Yi_upsamp = 10e-6;
    d_I_D_upsamp = 10e-6;
    
    % Indicator contours
    contLevels = [0.1 0.4];
    
    % Folder name
    folderName = 'Figure - Param Search';
    
    figName = 'ParamSearch_TwoElectrode';
    
    % Define membrane thresholds
    Th = [12.09e-3 6.30e-3*ones(1,length(rot)-1)];
    
elseif simulation == 6
    
    %% Four electrode parameter search
    
    % Define neurite rotations in xz plane
    d_rot = pi/6;
    rot = [0 0:d_rot:pi/2+pi/10000];
    
    % Define source locations
    Xi = [0e-6 0e-6 0e-6 0e-6];
    Yi_V = [10 (20:40:500)]*1e-6;
    Zi = [-300e-6 -100e-6 100e-6 300e-6];
    
    % Define electrode radius
    Ri = [50e-6 50e-6 50e-6 50e-6];
    
    % Define axon locations
    h_F = 100e-6;
    Ya = [-10e-6 -110e-6*ones(1,length(rot)-1)];
    
    % Define source amplitudes and pulse durations
    I_M = [-1e-6 -1e-6 -1e-6 -1e-6];
    I_D_V = [10 (20:40:500)]*1e-6;
    
    % Define location on neurite to calculate the membrane potential
    phi = 'magnitude';
    
    % Define post-simulation resampling
    d_theta_upsamp = pi/90;
    d_Yi_upsamp = 10e-6;
    d_I_D_upsamp = 10e-6;
    
    % Indicator contours
    contLevels = [0.1 0.4];
    
    % Folder name
    folderName = 'Figure - Param Search';
    
    figName = 'ParamSearch_FourElectrode';
    
    % Define membrane thresholds
    Th = [12.09e-3 6.30e-3*ones(1,length(rot)-1)];
    
elseif simulation == 6.1
    
    %% Four Electrode Simulation
    
    % Define neurite rotations in xz plane
    d_rot = pi/18;
    rot = [0 0:d_rot:pi/2+pi/10000];
    
    % Define source locations
    Xi = [0e-6 0e-6 0e-6 0e-6];
    Yi = [100e-6 100e-6 100e-6 100e-6];
    Zi = [-300e-6 -100e-6 100e-6 300e-6];
    
    % Define electrode radius
    Ri = [50e-6 50e-6 50e-6 50e-6];
    
    % Define axon locations
    h_F = 100e-6;
    Ya = [-10e-6 -110e-6*ones(1,length(rot)-1)];
    
    % Define source amplitudes and pulse durations
    I_M = [-50e-6 -50e-6 -50e-6 -50e-6];
    I_D = [100e-6 100e-6 100e-6 100e-6];
    
    % Define location on neurite to calculate the membrane potential
    phi = 'magnitude';
    
    % Should we scale the output so the AOP just hits threshold?
    scaleOut = true;
    
    % Folder name
    folderName = 'Figure - Param Search';
    
    % Sketchup filename
    filename = 'SimGeom_FourElectrode.jpg';
    
    % Plot filename
    figName = 'RotNeurite_FourElectrode_Y100D100';
    
elseif simulation == 6.2
    
    %% Four Electrode Simulation
    
    % Define neurite rotations in xz plane
    d_rot = pi/18;
    rot = [0 0:d_rot:pi/2+pi/10000];
    
    % Define source locations
    Xi = [0e-6 0e-6 0e-6 0e-6];
    Yi = [300e-6 300e-6 300e-6 300e-6];
    Zi = [-300e-6 -100e-6 100e-6 300e-6];
    
    % Define electrode radius
    Ri = [50e-6 50e-6 50e-6 50e-6];
    
    % Define axon locations
    h_F = 100e-6;
    Ya = [-10e-6 -110e-6*ones(1,length(rot)-1)];
    
    % Define source amplitudes and pulse durations
    I_M = [-50e-6 -50e-6 -50e-6 -50e-6];
    I_D = [100e-6 100e-6 100e-6 100e-6];
    
    % Define location on neurite to calculate the membrane potential
    phi = 'magnitude';
    
    % Should we scale the output so the AOP just hits threshold?
    scaleOut = true;
    
    % Folder name
    folderName = 'Figure - Param Search';
    
    % Sketchup filename
    filename = 'SimGeom_FourElectrode.jpg';
    
    % Plot filename
    figName = 'RotNeurite_FourElectrode_Y300D100';
    
elseif simulation == 6.3
    
    %% Four Electrode Simulation
    
    % Define neurite rotations in xz plane
    d_rot = pi/18;
    rot = [0 0:d_rot:pi/2+pi/10000];
    
    % Define source locations
    Xi = [0e-6 0e-6 0e-6 0e-6];
    Yi = [100e-6 100e-6 100e-6 100e-6];
    Zi = [-300e-6 -100e-6 100e-6 300e-6];
    
    % Define electrode radius
    Ri = [50e-6 50e-6 50e-6 50e-6];
    
    % Define axon locations
    h_F = 100e-6;
    Ya = [-10e-6 -110e-6*ones(1,length(rot)-1)];
    
    % Define source amplitudes and pulse durations
    I_M = [-50e-6 -50e-6 -50e-6 -50e-6];
    I_D = [300e-6 300e-6 300e-6 300e-6];
    
    % Define location on neurite to calculate the membrane potential
    phi = 'magnitude';
    
    % Should we scale the output so the AOP just hits threshold?
    scaleOut = true;
    
    % Folder name
    folderName = 'Figure - Param Search';
    
    % Sketchup filename
    filename = 'SimGeom_FourElectrode.jpg';
    
    % Plot filename
    figName = 'RotNeurite_FourElectrode_Y100D300';
    
elseif simulation == 6.4
    
    %% Four electrode parameter search
    
    % Define neurite rotations in xz plane
    d_rot = pi/6;
    rot = [0 0:d_rot:pi/2+pi/10000];
    
    % Define source locations
    Xi = [0e-6 0e-6 0e-6 0e-6];
    Yi_V = [10 (20:40:500)]*1e-6;
    Zi = [-300e-6 -100e-6 100e-6 300e-6];
    
    % Define electrode radius
    Ri = [50e-6 50e-6 50e-6 50e-6];
    
    % Define axon locations
    h_F = 100e-6;
    Ya = [-10e-6 -110e-6*ones(1,length(rot)-1)];
    
    % Define source amplitudes and pulse durations
    I_M = [-1e-6 -1e-6 -1e-6 -1e-6];
    I_D_V = [10 (20:40:500)]*1e-6;
    
    % Define location on neurite to calculate the membrane potential
    phi = 'magnitude';
    
    % Define post-simulation resampling
    d_theta_upsamp = pi/90;
    d_Yi_upsamp = 10e-6;
    d_I_D_upsamp = 10e-6;
    
    % Indicator contours
    contLevels = [0.1 0.4];
    
    % Folder name
    folderName = 'Figure - Param Search';
    
    figName = 'ParamSearch_FourElectrode_ratio1p5';
    
    % Define membrane thresholds
    Th = [12.09e-3 (12.09+6.30)/2*1e-3*ones(1,length(rot)-1)];
    
elseif simulation == 6.5
    
    %% Four electrode parameter search
    
    % Define neurite rotations in xz plane
    d_rot = pi/6;
    rot = [0 0:d_rot:pi/2+pi/10000];
    
    % Define source locations
    Xi = [0e-6 0e-6 0e-6 0e-6];
    Yi_V = [10 (20:40:500)]*1e-6;
    Zi = [-300e-6 -100e-6 100e-6 300e-6];
    
    % Define electrode radius
    Ri = [50e-6 50e-6 50e-6 50e-6];
    
    % Define axon locations
    h_F = 100e-6;
    Ya = [-10e-6 -110e-6*ones(1,length(rot)-1)];
    
    % Define source amplitudes and pulse durations
    I_M = [-1e-6 -1e-6 -1e-6 -1e-6];
    I_D_V = [10 (20:40:500)]*1e-6;
    
    % Define location on neurite to calculate the membrane potential
    phi = 'magnitude';
    
    % Define post-simulation resampling
    d_theta_upsamp = pi/90;
    d_Yi_upsamp = 10e-6;
    d_I_D_upsamp = 10e-6;
    
    % Indicator contours
    contLevels = [0.1 0.4];
    
    % Folder name
    folderName = 'Figure - Param Search';
    
    figName = 'ParamSearch_FourElectrode_ratio1';
    
    % Define membrane thresholds
    Th = [12.09e-3 12.09e-3*ones(1,length(rot)-1)];
    
elseif simulation == 7
    
    %% One Electrode Plane Simulations
    
    % Define neurite rotations in xz plane
    d_rot = pi/6;
    rot = [0 0:d_rot:pi+pi/10000];
    
    % Define source locations
    Xi = [0e-6];
    Yi = [100]*1e-6;
    Zi = [0e-6];
    
    % Define electrode radius
    Ri = [50e-6];
    
    % Define axon locations
    h_F = 100e-6;
    Ya = [-10e-6 -110e-6*ones(1,length(rot)-1)];
    
    % Define source amplitudes and pulse durations
    I_M = [-1e-6];
    I_D = [200e-6];
    
    % Define location on neurite to calculate the membrane potential
    phi = 'magnitude';
    
    % Define post-simulation resampling
    d_theta_upsamp = pi/90;
    
    % Should we scale the output so the AOP just hits threshold?
    scaleOut = 0.5;
    
    % Folder name
    folderName = 'Figure - Current vs Activation';
    
    figName = 'RotNeurite_OneElectrode_Plane100um';
    
elseif simulation == 8
    
    %% Two Electrode Plane Simulations
    
    % Define neurite rotations in xz plane
    d_rot = pi/6;
    rot = [0 0:d_rot:pi+pi/10000];
    
    % Define source locations
    Xi = [0e-6 0e-6];
    Yi = [100 100]*1e-6;
    Zi = [-100e-6 100e-6];
    
    % Define electrode radius
    Ri = [50e-6 50e-6];
    
    % Define axon locations
    h_F = 100e-6;
    Ya = [-10e-6 -110e-6*ones(1,length(rot)-1)];
    
    % Define source amplitudes and pulse durations
    I_M = [-1e-6 -1e-6];
    I_D = [200e-6 200e-6];
    
    % Define location on neurite to calculate the membrane potential
    phi = 'magnitude';
    
    % Define post-simulation resampling
    d_theta_upsamp = pi/90;
    
    % Should we scale the output so the AOP just hits threshold?
    scaleOut = 0.5;
    
    % Folder name
    folderName = 'Figure - Current vs Activation';
    
    figName = 'RotNeurite_TwoElectrode_Plane100um';
    
elseif simulation == 9
    
    %% Four Electrode Plane Simulations
    
    % Define neurite rotations in xz plane
    d_rot = pi/6;
    rot = [0 0:d_rot:pi+pi/10000];
    
    % Define source locations
    Xi = [0e-6 0e-6 0e-6 0e-6];
    Yi = [100 100 100 100]*1e-6;
    Zi = [-300e-6 -100e-6 100e-6 300e-6];
    
    % Define electrode radius
    Ri = [50e-6 50e-6 50e-6 50e-6];
    
    % Define axon locations
    h_F = 100e-6;
    Ya = [-10e-6 -110e-6*ones(1,length(rot)-1)];
    
    % Define source amplitudes and pulse durations
    I_M = [-1e-6 -1e-6 -1e-6 -1e-6];
    I_D = [200e-6 200e-6 200e-6 200e-6];
    
    % Define location on neurite to calculate the membrane potential
    phi = 'magnitude';
    
    % Define post-simulation resampling
    d_theta_upsamp = pi/90;
    
    % Should we scale the output so the AOP just hits threshold?
    scaleOut = 0.5;
    
    % Folder name
    folderName = 'Figure - Current vs Activation';
    
    figName = 'RotNeurite_FourElectrode_Plane100um';
    
elseif simulation == 10
    
    %% One, two and four Electrode Simulation
    
    % Define neurite rotations in xz plane
    d_rot = pi/6;
    rot = [0 0:d_rot:pi+pi/10000];
    
    % Define source locations
    Xi_V = {[0e-6], [0e-6], ...
        [0e-6 0e-6], [0e-6 0e-6], ...
        [0e-6 0e-6 0e-6 0e-6], [0e-6 0e-6 0e-6 0e-6]};
    Yi_V = {[100e-6], [300e-6], ...
        [100e-6 100e-6], [300e-6 300e-6], ...
        [100e-6 100e-6 100e-6 100e-6], [300e-6 300e-6 300e-6 300e-6]};
    Zi_V = {[0e-6], [0e-6], ...
        [-100e-6 100e-6], [-100e-6 100e-6], ...
        [-300e-6 -100e-6 100e-6 300e-6], [-300e-6 -100e-6 100e-6 300e-6]};
    
    % Define electrode radius
    Ri_V = {[50e-6], [50e-6], ...
        [50e-6 50e-6], [50e-6 50e-6], ...
        [50e-6 50e-6 50e-6 50e-6], [50e-6 50e-6 50e-6 50e-6]};
    
    % Define axon locations
    h_F = 100e-6;
    Ya = [-10e-6 -110e-6*ones(1,length(rot)-1)];
    
    % Define source amplitudes and pulse durations
    I_M_V = {[-1e-6], [-1e-6], ...
        [-1e-6 -1e-6], [-1e-6 -1e-6], ...
        [-1e-6 -1e-6 -1e-6 -1e-6], [-1e-6 -1e-6 -1e-6 -1e-6]};
    I_D_V = {[200e-6], [200e-6], ...
        [200e-6 200e-6], [200e-6 200e-6], ...
        [200e-6 200e-6 200e-6 200e-6], [200e-6 200e-6 200e-6 200e-6]};
    
    % Define location on neurite to calculate the membrane potential
    phi = 'magnitude';
    
    % Define post-simulation resampling
    d_theta_upsamp = pi/720;
    
    %Plot scale
    pScale = 'linear';
    crossing = 0.5;
    
    % Folder name
    folderName = 'Figure - Current vs Activation';
    
    % Plot filename
    figName = 'CurrentVsActAndRadius';
    
elseif simulation == 11
    
    % Eight electrode simulation
    
    % Define neurite rotations in xz plane
    d_rot = pi/18;
    rot = [0 0:d_rot:pi/2+pi/10000];
    
    % Define source locations
    Xi = [-100e-6 -100e-6 -100e-6 -100e-6 100e-6 100e-6 100e-6 100e-6];
    Yi = [100e-6 100e-6 100e-6 100e-6 100e-6 100e-6 100e-6 100e-6];
    Zi = [-300e-6 -100e-6 100e-6 300e-6 -300e-6 -100e-6 100e-6 300e-6];
    
    % Define electrode radius
    Ri = [50e-6 50e-6 50e-6 50e-6 50e-6 50e-6 50e-6 50e-6];
    
    % Define axon locations
    h_F = 100e-6;
    Ya = [-10e-6 -110e-6*ones(1,length(rot)-1)];
    
    % Define source amplitudes and pulse durations
    I_M = [-50e-6 -50e-6 -50e-6 -50e-6 -50e-6 -50e-6 -50e-6 -50e-6]/2;
    I_D = [200e-6 200e-6 200e-6 200e-6 200e-6 200e-6 200e-6 200e-6];
    
    % Define location on neurite to calculate the membrane potential
    phi = 'magnitude';
    
    % Should we scale the output so the AOP just hits threshold?
    scaleOut = true;
    
    % Folder name
    folderName = 'Figure - Non-ideal Conditions';
    
    % Plot filename
    figName = 'RotNeurite_EightElectrode';
    
elseif simulation == 11.1
    
    % Eight electrode simulation
    
    % Define neurite rotations in xz plane
    d_rot = pi/6;
    rot = [0 0:d_rot:pi+pi/10000];
    
    % Define source locations
    Xi = [-100e-6 -100e-6 -100e-6 -100e-6 100e-6 100e-6 100e-6 100e-6];
    Yi = [100e-6 100e-6 100e-6 100e-6 100e-6 1 00e-6 100e-6 100e-6];
    Zi = [-300e-6 -100e-6 100e-6 300e-6 -300e-6 -100e-6 100e-6 300e-6];
    
    % Define electrode radius
    Ri = [50e-6 50e-6 50e-6 50e-6 50e-6 50e-6 50e-6 50e-6];
    
    % Define axon locations
    h_F = 100e-6;
    Ya = [-10e-6 -110e-6*ones(1,length(rot)-1)];
    
    % Define source amplitudes and pulse durations
    I_M = [-50e-6 -50e-6 -50e-6 -50e-6 -50e-6 -50e-6 -50e-6 -50e-6]/2;
    I_D = [200e-6 200e-6 200e-6 200e-6 200e-6 200e-6 200e-6 200e-6];
    
    % Define location on neurite to calculate the membrane potential
    phi = 'magnitude';
    
    % Define post-simulation resampling
    d_theta_upsamp = pi/90;
    
    % Should we scale the output so the AOP just hits threshold?
    scaleOut = 0.5;
    
    % Folder name
    folderName = 'Figure - Non-ideal Conditions';
    
    % Plot filename
    figName = 'RotNeurite_EightElectrode_Plane';
    
elseif simulation == 11.2
    
    % Eight electrode simulation
    
    % Define neurite rotations in xz plane
    d_rot = pi/6;
    rot = [0 0:d_rot:pi+pi/10000];
    
    % Define source locations
    Xi_V = {[0e-6 0e-6 0e-6 0e-6], ...
        [-100e-6 -100e-6 -100e-6 -100e-6 100e-6 100e-6 100e-6 100e-6]};
    Yi_V = {[100e-6 100e-6 100e-6 100e-6], ...
        [100e-6 100e-6 100e-6 100e-6 100e-6 100e-6 100e-6 100e-6]};
    Zi_V = {[-300e-6 -100e-6 100e-6 300e-6], ...
        [-300e-6 -100e-6 100e-6 300e-6 -300e-6 -100e-6 100e-6 300e-6]};
    
    % Define electrode radius
    Ri_V = {[50e-6 50e-6 50e-6 50e-6], ...
        [50e-6 50e-6 50e-6 50e-6 50e-6 50e-6 50e-6 50e-6]};
    
    % Define axon locations
    h_F = 100e-6;
    Ya = [-10e-6 -110e-6*ones(1,length(rot)-1)];
    
    % Define source amplitudes and pulse durations
    I_M_V = {[-1e-6 -1e-6 -1e-6 -1e-6], ...
        [-1e-6 -1e-6 -1e-6 -1e-6 -1e-6 -1e-6 -1e-6 -1e-6]};
    I_D_V = {[200e-6 200e-6 200e-6 200e-6], ...
        [200e-6 200e-6 200e-6 200e-6 200e-6 200e-6 200e-6 200e-6]};
    
    % Define location on neurite to calculate the membrane potential
    phi = 'magnitude';
    
    % Define post-simulation resampling
    d_theta_upsamp = pi/720;
    
    %Plot scale
    pScale = 'linear';
    crossing = 0.5;
    
    % Folder name
    folderName = 'Figure - Non-ideal Conditions';
    
    % Plot filename
    figName = 'CurrentVsActAndRadius_EightElectrode';
    
elseif simulation == 12
    
    % 22.5 Deg electrode simulation
    
    % Define neurite rotations in xz plane
    d_rot = pi/18;
    rot = [0 0:d_rot:pi/2+pi/10000];
    
    % Define rotation of NFL
    rotNFL = 22.5/180*pi;
    rotNFL_mat = [cos(rotNFL) -sin(rotNFL); ...
        sin(rotNFL)  cos(rotNFL)];
    
    % Define source locations
    Xip = [-200e-6 0e-6 0e-6 200e-6];
    Yi = [100e-6 100e-6 100e-6 100e-6];
    Zip = [-300e-6 -100e-6 100e-6 300e-6];
    
    % Apply NFL rotation to electrodes
    XiZi = rotNFL_mat*[Xip; Zip];
    Xi = XiZi(1,:);
    Zi = XiZi(2,:);
    
    % Define electrode radius
    Ri = [50e-6 50e-6 50e-6 50e-6];
    
    % Define axon locations
    h_F = 100e-6;
    Ya = [-10e-6 -110e-6*ones(1,length(rot)-1)];
    
    % Define source amplitudes and pulse durations
    I_M = [-50e-6 -50e-6 -50e-6 -50e-6];
    I_D = [200e-6 200e-6 200e-6 200e-6];
    
    % Define location on neurite to calculate the membrane potential
    phi = 'magnitude';
    
    % Should we scale the output so the AOP just hits threshold?
    scaleOut = true;
    
    % Folder name
    folderName = 'Figure - Non-ideal Conditions';
    
    % Plot filename
    figName = 'RotNeurite_225DegAOP';
    
elseif simulation == 12.1
    
    % 22.5 Deg electrode simulation
    
    % Define neurite rotations in xz plane
    d_rot = pi/6;
    rot = [0 0:d_rot:pi+pi/10000];
    
    % Define rotation of NFL
    rotNFL = 22.5/180*pi;
    rotNFL_mat = [cos(rotNFL) -sin(rotNFL); ...
        sin(rotNFL)  cos(rotNFL)];
    
    % Define source locations
    Xip = [-200e-6 0e-6 0e-6 200e-6];
    Yi = [100e-6 100e-6 100e-6 100e-6];
    Zip = [-300e-6 -100e-6 100e-6 300e-6];
    
    % Apply NFL rotation to electrodes
    XiZi = rotNFL_mat*[Xip; Zip];
    Xi = XiZi(1,:);
    Zi = XiZi(2,:);
    
    % Define electrode radius
    Ri = [50e-6 50e-6 50e-6 50e-6];
    
    % Define axon locations
    h_F = 100e-6;
    Ya = [-10e-6 -110e-6*ones(1,length(rot)-1)];
    
    % Define source amplitudes and pulse durations
    I_M = [-50e-6 -50e-6 -50e-6 -50e-6];
    I_D = [200e-6 200e-6 200e-6 200e-6];
    
    % Define location on neurite to calculate the membrane potential
    phi = 'magnitude';
    
    % Define post-simulation resampling
    d_theta_upsamp = pi/90;
    
    % Should we scale the output so the AOP just hits threshold?
    scaleOut = 0.5;
    
    % Folder name
    folderName = 'Figure - Non-ideal Conditions';
    
    % Plot filename
    figName = 'RotNeurite_225DegAOP_Plane';
    
elseif simulation == 12.2
    
    % 22.5 degree electrode simulation
    
    % Define neurite rotations in xz plane
    d_rot = pi/6;
    rot = [0 0:d_rot:pi+pi/10000];
    
    % Define rotation of NFL
    rotNFL = 22.5/180*pi;
    rotNFL_mat = [cos(rotNFL) -sin(rotNFL); ...
        sin(rotNFL)  cos(rotNFL)];
    
    % Define source locations
    Xip = [-200e-6 0e-6 0e-6 200e-6];
    Yi = [100e-6 100e-6 100e-6 100e-6];
    Zip = [-300e-6 -100e-6 100e-6 300e-6];
    
    % Apply NFL rotation to electrodes
    XiZi = rotNFL_mat*[Xip; Zip];
    
    Xi_V = {[0e-6 0e-6 0e-6 0e-6],XiZi(1,:)};
    Yi_V = {[100e-6 100e-6 100e-6 100e-6],Yi};
    Zi_V = {[-300e-6 -100e-6 100e-6 300e-6],XiZi(2,:)};
    
    % Define electrode radius
    Ri_V = {[50e-6 50e-6 50e-6 50e-6], ...
        [50e-6 50e-6 50e-6 50e-6]};
    
    % Define axon locations
    h_F = 100e-6;
    Ya = [-10e-6 -110e-6*ones(1,length(rot)-1)];
    
    % Define source amplitudes and pulse durations
    I_M_V = {[-1e-6 -1e-6 -1e-6 -1e-6], ...
        [-1e-6 -1e-6 -1e-6 -1e-6]};
    I_D_V = {[200e-6 200e-6 200e-6 200e-6], ...
        [200e-6 200e-6 200e-6 200e-6]};
    
    % Define location on neurite to calculate the membrane potential
    phi = 'magnitude';
    
    % Define post-simulation resampling
    d_theta_upsamp = pi/720;
    
    %Plot scale
    pScale = 'linear';
    crossing = 0.5;
    
    % Folder name
    folderName = 'Figure - Non-ideal Conditions';
    
    % Plot filename
    figName = 'CurrentVsActAndRadius_225DegAOP';
    
elseif simulation == 13
    
    % 45 Deg offset electrode simulation
    
    % Define neurite rotations in xz plane
    d_rot = pi/18;
    rot = [0 0:d_rot:pi/2+pi/10000];
    
    % Define rotation of NFL
    rotNFL = 45/180*pi;
    rotNFL_mat = [cos(rotNFL) -sin(rotNFL); ...
        sin(rotNFL)  cos(rotNFL)];
    
    % Define source locations
    Xip = [-200e-6 -200e-6 0e-6 0e-6 200e-6 200e-6];
    Yi = [100e-6 100e-6 100e-6 100e-6 100e-6 100e-6];
    Zip = [-300e-6 -100e-6 -100e-6 100e-6 100e-6 300e-6];
    
    % Apply NFL rotation to electrodes
    XiZi = rotNFL_mat*[Xip; Zip];
    Xi = XiZi(1,:);
    Zi = XiZi(2,:);
    
    % Define electrode radius
    Ri = [50e-6 50e-6 50e-6 50e-6 50e-6 50e-6];
    
    % Define axon locations
    h_F = 100e-6;
    Ya = [-10e-6 -110e-6*ones(1,length(rot)-1)];
    
    % Define source amplitudes and pulse durations
    I_M = [-50e-6 -50e-6 -50e-6 -50e-6 -50e-6 -50e-6];
    I_D = [200e-6 200e-6 200e-6 200e-6 200e-6 200e-6];
    
    % Define location on neurite to calculate the membrane potential
    phi = 'magnitude';
    
    % Should we scale the output so the AOP just hits threshold?
    scaleOut = true;
    
    % Folder name
    folderName = 'Figure - Non-ideal Conditions';
    
    % Plot filename
    figName = 'RotNeurite_45DegOffsetAOP';
    
elseif simulation == 13.1
    
    % 45 Deg offset electrode simulation
    
    % Define neurite rotations in xz plane
    d_rot = pi/6;
    rot = [0 0:d_rot:pi+pi/10000];
    
    % Define rotation of NFL
    rotNFL = 45/180*pi;
    rotNFL_mat = [cos(rotNFL) -sin(rotNFL); ...
        sin(rotNFL)  cos(rotNFL)];
    
    % Define source locations
    Xip = [-200e-6 -200e-6 0e-6 0e-6 200e-6 200e-6];
    Yi = [100e-6 100e-6 100e-6 100e-6 100e-6 100e-6];
    Zip = [-300e-6 -100e-6 -100e-6 100e-6 100e-6 300e-6];
    
    % Apply NFL rotation to electrodes
    XiZi = rotNFL_mat*[Xip; Zip];
    Xi = XiZi(1,:);
    Zi = XiZi(2,:);
    
    % Define electrode radius
    Ri = [50e-6 50e-6 50e-6 50e-6 50e-6 50e-6];
    
    % Define axon locations
    h_F = 100e-6;
    Ya = [-10e-6 -110e-6*ones(1,length(rot)-1)];
    
    % Define source amplitudes and pulse durations
    I_M = [-50e-6 -50e-6 -50e-6 -50e-6 -50e-6 -50e-6];
    I_D = [200e-6 200e-6 200e-6 200e-6 200e-6 200e-6];
    
    % Define location on neurite to calculate the membrane potential
    phi = 'magnitude';
    
    % Define post-simulation resampling
    d_theta_upsamp = pi/90;
    
    % Should we scale the output so the AOP just hits threshold?
    scaleOut = 0.5;
    
    % Folder name
    folderName = 'Figure - Non-ideal Conditions';
    
    % Plot filename
    figName = 'RotNeurite_45DegOffsetAOP_Plane';
    
elseif simulation == 13.2
    
    % 45 degree offset electrode simulation
    
    % Define neurite rotations in xz plane
    d_rot = pi/6;
    rot = [0 0:d_rot:pi+pi/10000];
    
    % Define rotation of NFL
    rotNFL = 45/180*pi;
    rotNFL_mat = [cos(rotNFL) -sin(rotNFL); ...
        sin(rotNFL)  cos(rotNFL)];
    
    % Define source locations
    Xip = [-200e-6 -200e-6 0e-6 0e-6 200e-6 200e-6];
    Yi = [100e-6 100e-6 100e-6 100e-6 100e-6 100e-6];
    Zip = [-300e-6 -100e-6 -100e-6 100e-6 100e-6 300e-6];
    
    % Apply NFL rotation to electrodes
    XiZi = rotNFL_mat*[Xip; Zip];
    
    Xi_V = {[0e-6 0e-6 0e-6 0e-6],XiZi(1,:)};
    Yi_V = {[100e-6 100e-6 100e-6 100e-6],Yi};
    Zi_V = {[-300e-6 -100e-6 100e-6 300e-6],XiZi(2,:)};
    
    % Define electrode radius
    Ri_V = {[50e-6 50e-6 50e-6 50e-6], ...
        [50e-6 50e-6 50e-6 50e-6 50e-6 50e-6]};
    
    % Define axon locations
    h_F = 100e-6;
    Ya = [-10e-6 -110e-6*ones(1,length(rot)-1)];
    
    % Define source amplitudes and pulse durations
    I_M_V = {[-1e-6 -1e-6 -1e-6 -1e-6], ...
        [-1e-6 -1e-6 -1e-6 -1e-6 -1e-6 -1e-6]};
    I_D_V = {[200e-6 200e-6 200e-6 200e-6], ...
        [200e-6 200e-6 200e-6 200e-6 200e-6 200e-6]};
    
    % Define location on neurite to calculate the membrane potential
    phi = 'magnitude';
    
    % Define post-simulation resampling
    d_theta_upsamp = pi/720;
    
    %Plot scale
    pScale = 'linear';
    crossing = 0.5;
    
    % Folder name
    folderName = 'Figure - Non-ideal Conditions';
    
    % Plot filename
    figName = 'CurrentVsActAndRadius_45DegOffsetAOP';
    
elseif simulation == 14
    
    % Offset electrode simulation
    
    % Define neurite rotations in xz plane
    d_rot = pi/18;
    rot = [0 0:d_rot:pi/2+pi/10000];
    
    % Define rotation of NFL
    rotNFL = 0/180*pi;
    rotNFL_mat = [cos(rotNFL) -sin(rotNFL); ...
        sin(rotNFL)  cos(rotNFL)];
    
    % Define source locations
    Xip = [50e-6 50e-6 50e-6 50e-6 -150e-6 -150e-6 -150e-6 -150e-6];
    Yi = [100e-6 100e-6 100e-6 100e-6 100e-6 100e-6 100e-6 100e-6];
    Zip = [-300e-6 -100e-6 100e-6 300e-6 -300e-6 -100e-6 100e-6 300e-6];
    
    % Apply NFL rotation to electrodes
    XiZi = rotNFL_mat*[Xip; Zip];
    Xi = XiZi(1,:);
    Zi = XiZi(2,:);
    
    % Define electrode radius
    Ri = [50e-6 50e-6 50e-6 50e-6 50e-6 50e-6 50e-6 50e-6];
    
    % Define axon locations
    h_F = 100e-6;
    Ya = [-10e-6 -110e-6*ones(1,length(rot)-1)];
    
    % Define source amplitudes and pulse durations
    I_M = [-67e-6 -67e-6 -67e-6 -67e-6 -33e-6 -33e-6 -33e-6 -33e-6];
    I_D = [200e-6 200e-6 200e-6 200e-6 200e-6 200e-6 200e-6 200e-6];
    
    % Define location on neurite to calculate the membrane potential
    phi = 'magnitude';
    
    % Should we scale the output so the AOP just hits threshold?
    scaleOut = true;
    
    % Folder name
    folderName = 'Figure - Non-ideal Conditions';
    
    % Plot filename
    figName = 'RotNeurite_OffsetAOP';
    
elseif simulation == 14.1
    
    % Offset electrode simulation
    
    % Define neurite rotations in xz plane
    d_rot = pi/6;
    rot = [0 0:d_rot:pi+pi/10000];
    
    % Define rotation of NFL
    rotNFL = 0/180*pi;
    rotNFL_mat = [cos(rotNFL) -sin(rotNFL); ...
        sin(rotNFL)  cos(rotNFL)];
    
    % Define source locations
    Xip = [50e-6 50e-6 50e-6 50e-6 -150e-6 -150e-6 -150e-6 -150e-6];
    Yi = [100e-6 100e-6 100e-6 100e-6 100e-6 100e-6 100e-6 100e-6];
    Zip = [-300e-6 -100e-6 100e-6 300e-6 -300e-6 -100e-6 100e-6 300e-6];
    
    % Apply NFL rotation to electrodes
    XiZi = rotNFL_mat*[Xip; Zip];
    Xi = XiZi(1,:);
    Zi = XiZi(2,:);
    
    % Define electrode radius
    Ri = [50e-6 50e-6 50e-6 50e-6 50e-6 50e-6 50e-6 50e-6];
    
    % Define axon locations
    h_F = 100e-6;
    Ya = [-10e-6 -110e-6*ones(1,length(rot)-1)];
    
    % Define source amplitudes and pulse durations
    I_M = [-67e-6 -67e-6 -67e-6 -67e-6 -33e-6 -33e-6 -33e-6 -33e-6];
    I_D = [200e-6 200e-6 200e-6 200e-6 200e-6 200e-6 200e-6 200e-6];
    
    % Define location on neurite to calculate the membrane potential
    phi = 'magnitude';
    
    % Define post-simulation resampling
    d_theta_upsamp = pi/90;
    
    % Should we scale the output so the AOP just hits threshold?
    scaleOut = 0.5;
    
    % Folder name
    folderName = 'Figure - Non-ideal Conditions';
    
    % Plot filename
    figName = 'RotNeurite_OffsetAOP_Plane';
    
elseif simulation == 14.2
    
    % Offset electrode simulation
    
    % Define neurite rotations in xz plane
    d_rot = pi/6;
    rot = [0 0:d_rot:pi+pi/10000];
    
    % Define rotation of NFL
    rotNFL = 0/180*pi;
    rotNFL_mat = [cos(rotNFL) -sin(rotNFL); ...
        sin(rotNFL)  cos(rotNFL)];
    
    % Define source locations
    Xip = [50e-6 50e-6 50e-6 50e-6 -150e-6 -150e-6 -150e-6 -150e-6];
    Yi = [100e-6 100e-6 100e-6 100e-6 100e-6 100e-6 100e-6 100e-6];
    Zip = [-300e-6 -100e-6 100e-6 300e-6 -300e-6 -100e-6 100e-6 300e-6];
    
    % Apply NFL rotation to electrodes
    XiZi = rotNFL_mat*[Xip; Zip];
    
    Xi_V = {[0e-6 0e-6 0e-6 0e-6],XiZi(1,:)};
    Yi_V = {[100e-6 100e-6 100e-6 100e-6],Yi};
    Zi_V = {[-300e-6 -100e-6 100e-6 300e-6],XiZi(2,:)};
    
    % Define electrode radius
    Ri_V = {[50e-6 50e-6 50e-6 50e-6], ...
        [50e-6 50e-6 50e-6 50e-6 50e-6 50e-6 50e-6 50e-6]};
    
    % Define axon locations
    h_F = 100e-6;
    Ya = [-10e-6 -110e-6*ones(1,length(rot)-1)];
    
    % Define source amplitudes and pulse durations
    I_M_V = {[-1e-6 -1e-6 -1e-6 -1e-6], ...
        [-1.33e-6 -1.33e-6 -1.33e-6 -1.33e-6 -0.67e-6 -0.67e-6 -0.67e-6 -0.67e-6]};
    I_D_V = {[200e-6 200e-6 200e-6 200e-6], ...
        [200e-6 200e-6 200e-6 200e-6 200e-6 200e-6 200e-6 200e-6]};
    
    % Define location on neurite to calculate the membrane potential
    phi = 'magnitude';
    
    % Define post-simulation resampling
    d_theta_upsamp = pi/720;
    
    %Plot scale
    pScale = 'linear';
    crossing = 0.5;
    
    % Folder name
    folderName = 'Figure - Non-ideal Conditions';
    
    % Plot filename
    figName = 'CurrentVsActAndRadius_OffsetAOP';
    
elseif simulation == 15
    
    %% Four Electrode Current Optimisation
    
    % Define neurite rotations in xz plane
    d_rot = pi/6;
    rot = [0 0:d_rot:pi+pi/10000];
    
    % Define current ratio
    currentRatio = 0:0.1:1;
    
    % Define source locations
    Xi_V = mat2cell(repmat([0e-6 0e-6 0e-6 0e-6],length(currentRatio),1), ...
        ones(1,length(currentRatio)),4);
    Yi_V = mat2cell(repmat([100e-6 100e-6 100e-6 100e-6],length(currentRatio),1), ...
        ones(1,length(currentRatio)),4);
    Zi_V = mat2cell(repmat([-300e-6 -100e-6 100e-6 300e-6],length(currentRatio),1), ...
        ones(1,length(currentRatio)),4);
    
    % Define electrode radius
    Ri_V = mat2cell(repmat([50e-6 50e-6 50e-6 50e-6],length(currentRatio),1), ...
        ones(1,length(currentRatio)),4);
    
    % Define axon locations
    h_F = 100e-6;
    Ya = [-10e-6 -110e-6*ones(1,length(rot)-1)];
    
    % Define source amplitudes and pulse durations
    I_M_V = repmat([1e-6 1e-6 1e-6 1e-6],length(currentRatio),1).* ...
        [currentRatio' ones(length(currentRatio),2) currentRatio'];
    I_M_V = -I_M_V./repmat(sum(I_M_V,2),1,4)*1e-6;
    I_M_V = mat2cell(I_M_V,ones(1,length(currentRatio)),4);
    I_D_V = mat2cell(repmat([200e-6 200e-6 200e-6 200e-6],length(currentRatio),1), ...
        ones(1,length(currentRatio)),4);
    
    % Define location on neurite to calculate the membrane potential
    phi = 'magnitude';
    
    % Define post-simulation resampling
    d_theta_upsamp = pi/720;
    
    %Plot scale
    pScale = 'linear';
    crossing = 0.5;
    
    % Folder name
    folderName = 'Figure - Current Optim';
    
    % Plot filename
    figName = 'CurrentVsActAndRadius_Yi100um';
    
elseif simulation == 16
    
    %% Four Electrode Current Optimisation
    
    % Define neurite rotations in xz plane
    d_rot = pi/6;
    rot = [0 0:d_rot:pi+pi/10000];
    
    % Define current ratio
    currentRatio = 0:0.1:1;
    
    % Define source locations
    Xi_V = mat2cell(repmat([0e-6 0e-6 0e-6 0e-6],length(currentRatio),1), ...
        ones(1,length(currentRatio)),4);
    Yi_V = mat2cell(repmat([200e-6 200e-6 200e-6 200e-6],length(currentRatio),1), ...
        ones(1,length(currentRatio)),4);
    Zi_V = mat2cell(repmat([-300e-6 -100e-6 100e-6 300e-6],length(currentRatio),1), ...
        ones(1,length(currentRatio)),4);
    
    % Define electrode radius
    Ri_V = mat2cell(repmat([50e-6 50e-6 50e-6 50e-6],length(currentRatio),1), ...
        ones(1,length(currentRatio)),4);
    
    % Define axon locations
    h_F = 100e-6;
    Ya = [-10e-6 -110e-6*ones(1,length(rot)-1)];
    
    % Define source amplitudes and pulse durations
    I_M_V = repmat([1e-6 1e-6 1e-6 1e-6],length(currentRatio),1).* ...
        [currentRatio' ones(length(currentRatio),2) currentRatio'];
    I_M_V = -I_M_V./repmat(sum(I_M_V,2),1,4)*1e-6;
    I_M_V = mat2cell(I_M_V,ones(1,length(currentRatio)),4);
    I_D_V = mat2cell(repmat([200e-6 200e-6 200e-6 200e-6],length(currentRatio),1), ...
        ones(1,length(currentRatio)),4);
    
    % Define location on neurite to calculate the membrane potential
    phi = 'magnitude';
    
    % Define post-simulation resampling
    d_theta_upsamp = pi/720;
    
    %Plot scale
    pScale = 'linear';
    crossing = 0.5;
    
    % Folder name
    folderName = 'Figure - Current Optim';
    
    % Plot filename
    figName = 'CurrentVsActAndRadius_Yi200um';
    
elseif simulation == 17
    
    %% One Electrode Simulation
    
    % Define source locations
    Xi = [0e-6];
    Yi = [100e-6];
    Zi = [0e-6];
    
    % Define electrode radius
    Ri = 50e-6;
    
    % Define axon locations
    h_F = 100e-6;
    Ya = (100:-10:-200)*1e-6;
    Ya(1) = Ya(1) - 1e-6;
    
    % Define neurite rotations in xz plane
    rot = zeros(1,length(Ya));
    
    % Define source amplitudes and pulse durations
    I_M = [-50e-6];
    I_D = [200e-6];
    
    % Define location on neurite to calculate the membrane potential
    phi = 'magnitude';
    
    % Should we scale the output so the AOP just hits threshold?
    scaleOut = true;
    
    % Folder name
    folderName = 'Figure - Depth';
    
    % Plot filename
    figName = 'RotNeurite_OneElectrode';
    
elseif simulation== 17.8
    % Eight Electrode Simulation
    
    
    
    % Define source locations
    
    Xi = [-500e-6 -200e-6 -500e-6 -500e-6 500e-6 500e-6 500e-6 500e-6];
    
    Yi = [100e-6 100e-6 100e-6 100e-6 100e-6 100e-6 100e-6 100e-6]; %y =depth plane
    
    Zi = [-500e-6 -200e-6 200e-6 500e-6 -500e-6 -200e-6 200e-6 500e-6];
    
    
    
    % Define source amplitudes and pulse durations
    
    I_M = [-50e-6 -50e-6 -50e-6 -50e-6 -50e-6 -50e-6 -50e-6 -50e-6];
    
    I_D = [200e-6 200e-6 200e-6 200e-6 200e-6 200e-6 200e-6 200e-6];
    
    
    
    % Define electrode radius
    
    Ri = [50e-6 50e-6 50e-6 50e-6 50e-6 50e-6 50e-6 50e-6];
    
    
    
    % Define axon locations
    
    h_F = 100e-6;
    
    Ya = (100:-10:-200)*1e-6;
    
    Ya(1) = Ya(1) - 1e-6;
    
    
    
    % Define neurite rotations in xz plane
    
    rot = zeros(1,length(Ya));
    
    
    
    % Define location on neurite to calculate the membrane potential
    
    phi = 'magnitude';
    
    
    
    % Should we scale the output so the AOP just hits threshold?
    
    scaleOut = true;
    
    
    
    % Folder name
    
    folderName = 'Figure - Depth';
    
    
    
    % Plot filename
    
    figName = 'RotNeurite_8Electrode';
elseif simulation ==17.6
    %1 row simulation of electrodes according to Argus II devices
    
    % Define source locations
    
    Xi = [ -1725e-6 -575e-6 0e-6 575e-6 1150-6 1725-6 ];
    
    Yi = [100e-6 100e-6 100e-6 100e-6 100e-6 100e-6 ]; %y =depth plane
    %Yi = Yi.*4;
    % Zi = [-1150e-6 -575e-6 0e-6 575e-6 1150-6 1725-6 ];%
    Zi = [-0e-6 -0e-6 0e-6 0e-6 0e-6 0e-6 ];
    
    
    
    % Define source amplitudes and pulse durations
    
    %pairs
    %     I_M = [-50e-6 50-6 -50e-6 -50e-6 -50e-6 -50e-6];
    %     I_D = [200e-6 200e-6 200e-6 200e-6 200e-6 200e-6];
    
    %one electrode stimulation
    I_M = [0e-6 0e-6 50e-6 0e-6 0e-6 0e-6];
    I_D = [0e-6  0e-6 500e-6 0e-6  0e-6  0e-6 ];
    
    
    
    % Define electrode radius
    
    Ri = [200e-6 200e-6 200e-6 200e-6 200e-6 200e-6];
    
    
    
    % Define axon locations
    
    h_F = 100e-6;
    
    Ya = (100:-10:-200)*1e-6;
    
    Ya(1) = Ya(1) - 1e-6;
    
    
    
    % Define neurite rotations in xz plane
    
    rot = zeros(1,length(Ya));
    
    
    
    % Define location on neurite to calculate the membrane potential
    
    phi = 'magnitude';
    
    
    
    % Should we scale the output so the AOP just hits threshold?
    
    scaleOut = true;
    
    
    
    % Folder name
    
    folderName = 'Figure - Depth';
    
    
    
    % Plot filename
    
    figName = 'RotNeurite_ArgusOneElectrode';
elseif simulation 1017.6
    % 6 electrode simulation with varying rotations
    %1 row simulation of electrodes according to Argus II devices
    
    % Define source locations
    
    Xi = [ -1725e-6 -575e-6 0e-6 575e-6 1150-6 1725-6 ];
    
    Yi = [100e-6 100e-6 100e-6 100e-6 100e-6 100e-6 ]; %y =depth plane
    %Yi = Yi.*4;
    % Zi = [-1150e-6 -575e-6 0e-6 575e-6 1150-6 1725-6 ];%
    Zi = [-0e-6 -0e-6 0e-6 0e-6 0e-6 0e-6 ];
    
    
    
    % Define source amplitudes and pulse durations
    
    %pairs
    %     I_M = [-50e-6 50-6 -50e-6 -50e-6 -50e-6 -50e-6];
    %     I_D = [200e-6 200e-6 200e-6 200e-6 200e-6 200e-6];
    
    %one electrode stimulation
    I_M = [0e-6 0e-6 50e-6 0e-6 0e-6 0e-6];
    I_D = [0e-6  0e-6 500e-6 0e-6  0e-6  0e-6 ];
    
    
    
    % Define electrode radius
    
    Ri = [200e-6 200e-6 200e-6 200e-6 200e-6 200e-6];
    
    
    
    % Define axon locations
    
    h_F = 100e-6;
    
    Ya = (100:-10:-200)*1e-6;
    
    Ya(1) = Ya(1) - 1e-6;
    
    
    
    % Define neurite rotations in xz plane
    
    
    % Define neurite rotations in xz plane
    d_rot = pi/12;
    rot = [0 0:d_rot:pi+pi/10000];
%     rot= linspace(-pi,pi, 10);
    
    
    % Define location on neurite to calculate the membrane potential
    
    phi = 'magnitude';
    
    
    
    % Should we scale the output so the AOP just hits threshold?
    
    scaleOut = true;
    
    
    
    % Folder name
    
    folderName = 'Figure - Depth';
    
    
    
    % Plot filename
    
    figName = 'RotNeurite_ArgusOneElectrodeOfOneRow';
elseif simulation == 18
    
    %% One Electrode Simulation
    
    % Define source locations
    Xi = [0e-6];
    Yi = [50e-6];
    Zi = [0e-6];
    
    % Define electrode radius
    Ri = 50e-6;
    
    % Define axon locations
    h_F = 150e-6;
    Ya = (50:-10:-250)*1e-6;
    Ya(1) = Ya(1) - 1e-6;
    
    % Define neurite rotations in xz plane
    rot = zeros(1,length(Ya));
    
    % Define source amplitudes and pulse durations
    I_M = [-50e-6];
    I_D = [200e-6];
    
    % Define location on neurite to calculate the membrane potential
    phi = 'magnitude';
    
    % Should we scale the output so the AOP just hits threshold?
    scaleOut = true;
    
    % Folder name
    folderName = 'Figure - Depth';
    
    % Plot filename
    figName = 'RotNeurite_OneElectrode';
    
elseif simulation == 19
    
    %% One Electrode Simulation
    
    % Define source locations
    Xi = [0e-6];
    Yi = [100e-6];
    Zi = [0e-6];
    
    % Define electrode radius
    Ri = 50e-6;
    
    % Define axon locations
    h_F = 10e-6;
    Ya = (100:-10:-200)*1e-6;
    Ya(1) = Ya(1) - 1e-6;
    
    % Define neurite rotations in xz plane
    rot = zeros(1,length(Ya));
    
    % Define source amplitudes and pulse durations
    I_M = [-50e-6];
    I_D = [200e-6];
    
    % Define location on neurite to calculate the membrane potential
    phi = 'magnitude';
    
    % Should we scale the output so the AOP just hits threshold?
    scaleOut = true;
    
    % Folder name
    folderName = 'Figure - Depth';
    
    % Plot filename
    figName = 'RotNeurite_OneElectrode';
    
elseif simulation== 20
    %% Four Electrode Depth Simulation
    % Define neurite rotations in xz plane
    d_rot = pi/6;
    rot = [0 0:d_rot:pi+pi/10000];
    
    
    % Define source locations
    Xi = [-100e-6 -100e-6 -100e-6 -100e-6 100e-6 100e-6 100e-6 100e-6];
    Yi = [100e-6 100e-6 100e-6 100e-6 100e-6 100e-6 100e-6 100e-6];
    Zi = [-300e-6 -100e-6 100e-6 300e-6 -300e-6 -100e-6 100e-6 300e-6];
    
    % Define electrode radius
    Ri = [50e-6 50e-6 50e-6 50e-6 50e-6 50e-6 50e-6 50e-6];
    
    % Define axon locations
    h_F = 100e-6;
    Ya = [-10e-6 -110e-6*ones(1,length(rot)-1)];
    
    % Define source amplitudes and pulse durations
    I_M = [-50e-6 -50e-6 -50e-6 -50e-6 -50e-6 -50e-6 -50e-6 -50e-6]/2;
    I_D = [200e-6 200e-6 200e-6 200e-6 200e-6 200e-6 200e-6 200e-6];
    
    % Define location on neurite to calculate the membrane potential
    phi = 'magnitude';
    
    % Define post-simulation resampling
    d_theta_upsamp = pi/90;
    
    %Plot scale
    pScale = 'linear';
    crossing = 0.5;
    
    
    
    % Define axon locations
    h_F = 100e-6;
    Ya = (100:-10:-200)*1e-6;
    Ya(1) = Ya(1) - 1e-6;
    
    % Define neurite rotations in xz plane
    rot = zeros(1,length(Ya));
    
    % Define source amplitudes and pulse durations
    I_M = [-50e-6];
    I_D = [200e-6];
    
    % Define location on neurite to calculate the membrane potential
    phi = 'magnitude';
    
    % Should we scale the output so the AOP just hits threshold?
    scaleOut = true;
    
    % Folder name
    folderName = 'Figure - Depth';
    
    % Plot filename
    figName = 'RotNeurite_FourElectrode'
elseif simulation == 1001
    
    %% One Electrode Simulation
    
    % Define neurite rotations in xz plane
    d_rot = pi/18;
    rot = [0 0:d_rot:pi/2+pi/10000];
    
    % Define source locations
    Xi = [0e-6];
    Yi = [200e-6];
    Zi = [0e-6];
    
    % Define electrode radius
    Ri = 50e-6;
    
    % Define axon locations
    h_F = 100e-6;
    Ya = [-10e-6 -110e-6*ones(1,length(rot)-1)];
    
    % Define source amplitudes and pulse durations
    I_M = [-300e-6];
    I_D = [200e-6];
    
    % Define location on neurite to calculate the membrane potential
    phi = 'magnitude';
    
    % Should we scale the output so the AOP just hits threshold?
    scaleOut = false;
    
    % Folder name
    folderName = 'Collaboration';
    
    % Sketchup filename
    filename = 'SimGeom_OneElectrode.jpg';
    
    % Plot filename
    figName = 'RotNeurite_OneElectrode';
    
elseif simulation == 1002
    
    %% One Electrode Plane Simulations
    
    % Define neurite rotations in xz plane
    d_rot = pi/6;
    rot = [0 0:d_rot:pi+pi/10000];
    
    % Define source locations
    Xi = [0e-6];
    Yi = [200]*1e-6;
    Zi = [0e-6];
    
    % Define electrode radius
    Ri = [50e-6];
    
    % Define axon locations
    h_F = 20e-6;
    Ya = [-10e-6 -110e-6*ones(1,length(rot)-1)];
    
    % Define source amplitudes and pulse durations
    I_M = [-300e-6];
    I_D = [200e-6];
    
    % Define location on neurite to calculate the membrane potential
    phi = 'magnitude';
    
    % Define post-simulation resampling
    d_theta_upsamp = pi/90;
    
    % Should we scale the output so the AOP just hits threshold?
    scaleOut = false;
    
    % Folder name
    folderName = 'Collaboration';
    
    figName = 'Plane200um_Normal';
end
