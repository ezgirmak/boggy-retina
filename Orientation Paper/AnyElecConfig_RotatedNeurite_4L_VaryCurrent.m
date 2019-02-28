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

% simulation = 10;

SimulationParameters

% Define membrane thresholds
Th = [12.09e-3 6.30e-3*ones(1,length(rot)-1)];

%% Error checking

if max(rot) ~= pi
    error('"rot" should vary between 0 and pi')
end

%% Calculate the longitudinal and transverse components separately

h = waitbar(0,'Calculating Rotated Membrane Potential');

Vm_max = zeros(length(rot),length(Xi_V),length(X),length(Z));
iterations = length(rot)*length(Xi_V);
disp(['Iterations: ' num2str(iterations)])
count = 0;

for i = 1:length(Xi_V)
    Xi = Xi_V{i};
    Yi = Yi_V{i};
    Zi = Zi_V{i};
    Ri = Ri_V{i};
    I_M = I_M_V{i};
    I_D = I_D_V{i};
    
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
        
        Vm_max(k,i,:,:) = squeeze(max(abs(Vm_rot),[],3));
        
        count = count + 1;
        waitbar(count/iterations,h)
        
    end
end
close(h)

%% Pull out AOP and AIS components

Vm_max_AOP = Vm_max(1,:,:,:);

Vm_max = Vm_max(2:end,:,:,:);

rot_AOP = rot(1);
Ya_AOP = Ya(1);
rot = rot(2:end);
Ya = Ya(2:end);

%% Upsample w.r.t. orientation values

% Interestingly, Vm_max vs rot is predicted (almost) perfectly by a cubic
% polynomial, so we use this to upsample
% TO-DO: it is implemented perfectly by a sum of sines and cosines, so use
% that
rot_upsamp = 0:d_theta_upsamp:max(rot);

Vm_max_upsamp = interp1(rot,Vm_max,rot_upsamp,'pchip');

clear Vm_max

Vm_max_max_upsamp = squeeze(max(max(Vm_max_upsamp,[],4),[],3));
Vm_max_max_AOP = squeeze(max(max(Vm_max_AOP,[],4),[],3));

%% Get mean preferential activation and radius of activation

currentVec = 0:10e-6:2000e-6;

meanPrefAct = zeros(length(Xi_V),length(currentVec));
actRadius = zeros(length(Xi_V),length(currentVec));

h = waitbar(0,'Get pref. act.');
iterations = length(Xi_V)*length(currentVec);
count = 0;

for j = 1:length(Xi_V)
    for i = 1:length(currentVec)
        meanPrefAct(j,i) = squeeze(mean(Vm_max_max_upsamp(:,j)*currentVec(i)/sum(-I_M_V{j}) > Th(2),1));
        
        meanPrefActPlane = squeeze(mean(Vm_max_upsamp(:,j,:,:)*currentVec(i)/sum(-I_M_V{j}) > Th(2),1));
        [nzRows, nzCols] = find(meanPrefActPlane);
        if ~isempty(nzRows)
            actRadius(j,i) = max(sqrt(X(nzRows).^2 + Z(nzCols).^2));
        end
        
        count = count + 1;
        waitbar(count/iterations,h)
    end
end
close(h)

%% Create stacked plot

% smooth points
for i = 1:length(Xi_V)
    
    indsgt0 = meanPrefAct(i,:) > 0;
    indslt1 = meanPrefAct(i,:) < 1;
    meanPrefAct(i,indsgt0 & indslt1) = smooth(meanPrefAct(i,indsgt0 & indslt1),5);
    actRadius(i,indsgt0) = smooth(actRadius(i,indsgt0),5);
    
end

fig = figure('Units','centimeters','Position',[2 5 25 20],'Color','w');
set(fig, 'PaperPosition',[2 5 25 20])

if length(Xi_V) == 6
    ColorSet = linspecer(length(Xi_V)/2,'qualitative');
else
    ColorSet = linspecer(length(Xi_V),'qualitative');
end

ax1 = subplot(2,1,1);
set(ax1, 'ColorOrder', ColorSet);
set(ax1,'XScale',pScale)
hold on

AOPActLevel = zeros(1,length(Xi_V));
legendStr = cell(1,length(Xi_V));
crossingPointX = zeros(1, length(Xi_V));
crossingPointY = zeros(1, length(Xi_V));

for i = 1:length(Xi_V)
    if length(Xi_V) == 6
        colInd = ceil(i/2-0.1);
    else
        colInd = i;
    end
    ax1.ColorOrderIndex = colInd;
    
    AOPActLevel(i) = Th(1)./Vm_max_max_AOP(i)*sum(-I_M_V{i});
    
    meanPrefActTemp = meanPrefAct(i,(find(meanPrefAct(i,:) > 0,1)-1):find(meanPrefAct(i,:) == 1,1));
    currentVecTemp = currentVec((find(meanPrefAct(i,:) > 0,1)-1):find(meanPrefAct(i,:) == 1,1));

    h(i) = plot(currentVecTemp*1e6, ...
        meanPrefActTemp*100,'LineWidth', 3);

    plot(currentVecTemp(currentVecTemp > AOPActLevel(i))*1e6, ...
        meanPrefActTemp(currentVecTemp > AOPActLevel(i))*100, 'w','LineWidth', 3)
    ax1.ColorOrderIndex = colInd;
    plot(currentVecTemp(currentVecTemp > AOPActLevel(i))*1e6, ...
        meanPrefActTemp(currentVecTemp > AOPActLevel(i))*100, '--','LineWidth', 1)
    
    crossUpInd = find(meanPrefAct(i,:) >= crossing,1);
    crossDownInd = length(currentVec) + 1 - find(meanPrefAct(i,end:-1:1) < crossing,1);
    crossUpDist = meanPrefAct(i,crossUpInd) - crossing;
    crossDownDist = crossing - meanPrefAct(i,crossDownInd);
    crossingPointX(i) = (currentVec(crossUpInd)*crossDownDist + ...
        currentVec(crossDownInd)*crossUpDist)/(crossUpDist +crossDownDist);
    crossingPointY(i) = (actRadius(i,crossUpInd)*crossDownDist + ...
        actRadius(i,crossDownInd)*crossUpDist)/(crossUpDist +crossDownDist);
   
    if length(Xi_V) == 6
        legendStr{i} = [num2str(length(Xi_V{i})) ' electrode'];
    else
        legendStr{i} = [num2str(I_M_V{i}(1)/I_M_V{i}(2)) ' current ratio'];
    end
end
for i = 1:length(Xi_V)
    if uint8(mean(Yi_V{i})*1e6) == 100
        %         ln1 = line([0 crossingPointX(i)*1e6], crossing*100*ones(1,2));
        %         set(ln1,'LineStyle','--','Color','k')
        %         ln2 = line(crossingPointX(i)*1e6*ones(1,2), [crossing*100 0]);
        %         set(ln2,'LineStyle','--','Color','k')
        plot(crossingPointX(i)*1e6, crossing*100, 'ks', 'MarkerFaceColor','k', ...
            'MarkerSize',3)
    end
end
% plot([0 1e4],[100 100],'k')
ylabel('GCL activation level (%)')
xlim([currentVec(2) currentVec(end)]*1e6)
set(gca,'fontsize',16)
if length(Xi_V) == 6
    l = legend(h(2:2:end),legendStr(2:2:end),'Location','SouthEast');
else
    l = legend(h,legendStr,'Location','SouthEast');
    xlim([currentVec(2) currentVec(end)/2]*1e6)
end
set(l,'fontsize',12)
set(l,'Box','off')

ax2 = subplot(2,1,2);
set(ax2, 'ColorOrder', ColorSet);
set(ax2,'XScale',pScale)

hold on
for i = 1:length(Xi_V)
    if length(Xi_V) == 6
        colInd = ceil(i/2-0.1);
    else
        colInd = i;
    end
    ax2.ColorOrderIndex = colInd;
    
    actRadiusTemp = actRadius(i,(find(meanPrefAct(i,:) > 0,1)-1):end);
    currentVecTemp = currentVec((find(meanPrefAct(i,:) > 0,1)-1):end);
    
    h(i) = plot(currentVecTemp*1e6, actRadiusTemp*1e6, 'LineWidth', 3);
    plot(currentVecTemp(currentVecTemp > AOPActLevel(i))*1e6, ...
        actRadiusTemp(currentVecTemp > AOPActLevel(i))*1e6, 'w', 'LineWidth', 3)
    ax2.ColorOrderIndex = colInd;
    plot(currentVecTemp(currentVecTemp > AOPActLevel(i))*1e6, ...
        actRadiusTemp(currentVecTemp > AOPActLevel(i))*1e6, '--', 'LineWidth', 1.5)
end
for i = 1:length(Xi_V)
    if uint8(mean(Yi_V{i})*1e6) == 100
        %         ln1 = line(crossingPointX(i)*1e6*ones(1,2), [10000 crossingPointY(i)*1e6]);
        %         set(ln1,'LineStyle','--','Color','k')
        %         ln2 = line([0 crossingPointX(i)*1e6], crossingPointY(i)*1e6*ones(1,2));
        %         set(ln2,'LineStyle','--','Color','k')
        plot(crossingPointX(i)*1e6, crossingPointY(i)*1e6, 'ks', 'MarkerFaceColor','k', ...
            'MarkerSize',3)
    end
end
% plot(currentVec(end)*ones(1,2)*1e6,[0 2000],'k')
ylabel('Activation radius (\mum)')
xlabel('Total stimulus current (\muA)')
xlim([currentVec(2) currentVec(end)]*1e6)
if length(Xi_V) ~= 6
    xlim([currentVec(2) currentVec(end)/2]*1e6)
end
ylim([0 1499])
set(gca,'fontsize',16)

ax1Pos = get(ax1,'Position');
ax2Pos = get(ax2,'Position');

% set(ax1,'Box','on')
% set(ax2,'Box','on')
set(ax1,'XTickLabel','')
set(ax1,'Position', [ax1Pos(1) (ax1Pos(2)+ax2Pos(2)+ax2Pos(4))/2 ...
    ax1Pos(3) ax1Pos(4)+(ax1Pos(2)-ax2Pos(2)-ax2Pos(4))/2]);
set(ax2,'Position', [ax2Pos(1) ax2Pos(2) ...
    ax2Pos(3) ax2Pos(4)+(ax1Pos(2)-ax2Pos(2)-ax2Pos(4))/2]);

saveas(fig, ['Figures/' folderName '/' figName],'fig')
saveas(fig, ['Figures/' folderName '/' figName],'jpeg')
saveas(fig, ['Figures/' folderName '/' figName],'epsc')

fig = figure('Units','centimeters','Position',[2 5 25 10.7],'Color','w');
set(fig, 'PaperPosition',[2 5 25 10.7])

ax = axes;
set(ax, 'ColorOrder', ColorSet);
set(ax,'XScale',pScale)
hold on

for i = 1:length(Xi_V)
    if length(Xi_V) == 6
        colInd = ceil(i/2-0.1);
    else
        colInd = i;
    end
    ax.ColorOrderIndex = colInd;
    
    meanPrefActTemp = meanPrefAct(i,1:find(meanPrefAct(i,:) == 1,1));
    actRadiusTemp = actRadius(i,1:find(meanPrefAct(i,:) == 1,1));
    currentVecTemp = currentVec(1:find(meanPrefAct(i,:) == 1,1));
    
    h(i) = plot(actRadiusTemp*1e6, meanPrefActTemp*100, 'LineWidth', 3);
    uistack(h(i),'down');
    plot(actRadiusTemp(currentVecTemp > AOPActLevel(i))*1e6, ...
        meanPrefActTemp(currentVecTemp > AOPActLevel(i))*100, 'w', 'LineWidth', 3)
    ax.ColorOrderIndex = colInd;
    plot(actRadiusTemp(currentVecTemp > AOPActLevel(i))*1e6, ...
        meanPrefActTemp(currentVecTemp > AOPActLevel(i))*100, '--', 'LineWidth', 1.5)
end
for i = 1:length(Xi_V)
    if mean(Yi_V{i}) == 100e-6 || mean(Yi_V{i}) == 200e-6
        ln1 = line([0 crossingPointY(i)*1e6], crossing*100*ones(1,2));
        set(ln1,'LineStyle','--','Color','k')
        ln2 = line(crossingPointY(i)*1e6*ones(1,2), [crossing*100 0]);
        set(ln2,'LineStyle','--','Color','k')
        plot(crossingPointY(i)*1e6, crossing*100, 'ks', 'MarkerFaceColor','k', ...
            'MarkerSize',3)
    end
end
xlabel('Activation radius (\mum)')
ylabel('GCL activation level (%)')
% title('Activation level vs. radius')
% ax.Box = 'on';
set(gca,'fontsize',16)

saveas(fig, ['Figures/' folderName '/Inset' figName],'fig')
saveas(fig, ['Figures/' folderName '/Inset' figName],'jpeg')
saveas(fig, ['Figures/' folderName '/Inset' figName],'epsc')

