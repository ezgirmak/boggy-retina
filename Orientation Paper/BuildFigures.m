addpath([fileparts(pwd()) '\Utilities'])
addpath([fileparts(pwd()) '\Utilities\altmany-export_fig-04ca93c'])
addpath([fileparts(pwd()) '\Utilities\altmany-export_fig-04ca93c\xpdf-tools-win-4.00'])
addpath(fileparts(pwd()))

set(0,'defaultfigurecolor',[1 1 1])

%% Figure - Orientation Distribution

% Params
fontSize = 12;

figOrientDist = openfig('Figures\Figure - Orientation Distribution\LineOrientationDistribution.fig');
figOrientDist.Children.Children(1).HorizontalAlignment = 'center';
figOrientDist.Children.Children(1).Position = [68 0.26 0];
figOrientDist.Children.Children(1).String = {'Uniform distribution', 'of orientations'};
annotation('arrow',[0.78 0.78],[0.47 0.41], 'Color', [0 0 0], ...
    'LineWidth',1)

set(figOrientDist,'Name','Figure - Orientation Distribution')
xlabel('Orientation difference (degrees)')
colors = [0 0 0; 0 0 1; 1 0 0; 0 0.5 0];
count1 = 0;
count2 = 0;
for i = 1:numel(figOrientDist.CurrentAxes.Children)
    if strcmp(figOrientDist.CurrentAxes.Children(i).Type,'line')
        count1 = count1 + 1;
        figOrientDist.CurrentAxes.Children(i).Color = colors(count1,:);
        figOrientDist.CurrentAxes.Children(i).MarkerFaceColor = colors(count1,:);
        figOrientDist.CurrentAxes.Children(i).MarkerFaceColor = colors(count1,:);
        figOrientDist.CurrentAxes.Children(i).LineStyle = '-';
    end
    if strcmp(figOrientDist.CurrentAxes.Children(i).Type,'text')
        count2 = count2 + 1;
        figOrientDist.CurrentAxes.Children(i).Color = colors(count2,:);
    end
end
figOrientDist.Position = figOrientDist.Position.*[1 1 1 1.15];
figOrientDist.PaperPosition = figOrientDist.Position;
figOrientDist.CurrentAxes.Position = figOrientDist.CurrentAxes.Position.*[1 0.9 1 0.93];
axInvisible = axes('color','none','Ytick',0:0.1:0.5,'Position', ...
    figOrientDist.CurrentAxes.Position,'XTick',0:15:90,'XTickLabel',[], ...
    'TickDir','in','YTickLabel',[]);
set(axInvisible,'Xlim',get(figOrientDist.Children(2),'Xlim'), ...
    'Ylim',get(figOrientDist.Children(2),'Ylim'))
set(figOrientDist.Children(2),'XTick',[7.5:15:82.5], ...
    'XTickLabel',{'0-15','15-30','30-45','45-60','60-75','75-90'}, ...
    'TickLength',[0 0.0250])
text(-14,0.54,'(a)','FontSize',fontSize)

% print('Figures\Figure - Orientation Distribution\Final', '-depsc')
print('Figures\Figure - Orientation Distribution\Final','-djpeg','-r100')
export_fig 'Figures\Figure - Orientation Distribution\Final' -pdf -m8 -q101

figOrientFit = openfig('Figures\Figure - Orientation Distribution\DistOrientationDistribution.fig');
figOrientFit.Children.Children(1).String = '';
set(figOrientFit,'Name','Figure - Orientation Distribution')
xlabel('Orientation difference (degrees)')
colors = [0 0 0; 0 0 1; 1 0 0; 0 0.5 0; 0 0 1; 0 0 1; 1 0 0; 1 0 0; 0 0.5 0; 0 0.5 0];
count1 = 0;
count2 = 0;
for i = 1:numel(figOrientFit.CurrentAxes.Children)
    if strcmp(figOrientFit.CurrentAxes.Children(i).Type,'line')
        count1 = count1 + 1;
        figOrientFit.CurrentAxes.Children(i).Color = colors(count1,:);
        figOrientFit.CurrentAxes.Children(i).MarkerFaceColor = colors(count1,:);
    end
    if strcmp(figOrientFit.CurrentAxes.Children(i).Type,'text')
        count2 = count2 + 1;
        figOrientFit.CurrentAxes.Children(i).Color = colors(count2,:);
    end
    if strcmp(figOrientFit.CurrentAxes.Children(i).Type,'patch')
        figOrientFit.CurrentAxes.Children(i).FaceColor = colors(count2,:);
    end
end
for i = numel(figOrientFit.CurrentAxes.Children):-1:1
    if strcmp(figOrientFit.CurrentAxes.Children(i).Type,'patch')
        cl = figOrientFit.CurrentAxes.Children(i).FaceColor;
    else
        cl = figOrientFit.CurrentAxes.Children(i).Color;
    end
    if isequal(colors(1,:), cl)
        uistack(figOrientFit.CurrentAxes.Children(i), 'up', 100)
    elseif isequal(colors(2,:), cl)
        uistack(figOrientFit.CurrentAxes.Children(i), 'up', 50)
    elseif isequal(colors(4,:), cl)
        uistack(figOrientFit.CurrentAxes.Children(i), 'down', 50)
    end
end
figOrientFit.Position = figOrientFit.Position.*[1 1 1 1.15];
figOrientFit.PaperPosition = figOrientFit.Position;
figOrientFit.CurrentAxes.Position = figOrientFit.CurrentAxes.Position.*[1 0.9 1 0.93];
axInvisible = axes('color','none','Ytick',0:0.1:0.5,'Position', ...
    figOrientFit.CurrentAxes.Position,'XTick',0:15:90,'XTickLabel',[], ...
    'TickDir','in','YTickLabel',[]);
set(axInvisible,'Xlim',get(figOrientFit.Children(2),'Xlim'), ...
    'Ylim',get(figOrientFit.Children(2),'Ylim'))
set(figOrientFit.Children(2),'XTick',[7.5:15:82.5], ...
    'XTickLabel',{'0-15','15-30','30-45','45-60','60-75','75-90'}, ...
    'TickLength',[0 0.0250])
text(-14,0.54,'(b)','FontSize',fontSize)
uistack(figOrientFit.Children(2).Children(5), 'top')
uistack(figOrientFit.Children(2).Children(4), 'down', 1)
figOrientFit.Children(2).Children(1).LineStyle = '-';

% print('Figures\Figure - Orientation Distribution\FinalFit', '-depsc')
print('Figures\Figure - Orientation Distribution\FinalFit','-djpeg','-r100')
figOrientFit.PaperPosition(1:2) = [0 0];
print('Figures\Figure - Orientation Distribution\FinalFit','-dpdf', '-r1200')

figOrientDistZ = openfig('Figures\Figure - Orientation Distribution\LineZOrientationDistribution.fig');
set(figOrientDistZ,'Name','Figure - Orientation Distribution')
xlabel('Orientation difference (degrees)')
colors = [0 0 1; 1 0 0; 0 0.5 0];
count1 = 0;
count2 = 0;
for i = 1:numel(figOrientDistZ.CurrentAxes.Children)
    if strcmp(figOrientDistZ.CurrentAxes.Children(i).Type,'line')
        count1 = count1 + 1;
        figOrientDistZ.CurrentAxes.Children(i).Color = colors(count1,:);
        figOrientDistZ.CurrentAxes.Children(i).MarkerFaceColor = colors(count1,:);
    end
    if strcmp(figOrientDistZ.CurrentAxes.Children(i).Type,'text')
        count2 = count2 + 1;
        figOrientDistZ.CurrentAxes.Children(i).Color = colors(count2,:);
    end
end
figOrientDistZ.Position = figOrientDist.Position;
figOrientDistZ.PaperPosition = figOrientDist.Position;
figOrientDistZ.CurrentAxes.Position = figOrientDist.CurrentAxes.Position;
axInvisible = axes('color','none','Ytick',0:0.1:0.5,'Position', ...
    figOrientDistZ.CurrentAxes.Position,'XTick',0:15:90,'XTickLabel',[], ...
    'TickDir','in','YTickLabel',[]);
set(axInvisible,'Xlim',get(figOrientDistZ.Children(2),'Xlim'), ...
    'Ylim',get(figOrientDistZ.Children(2),'Ylim'))
set(figOrientDistZ.Children(2),'XTick',[7.5:15:82.5], ...
    'XTickLabel',{'0-15','15-30','30-45','45-60','60-75','75-90'}, ...
    'TickLength',[0 0.0250])
text(-14,1.08,'(c)','FontSize',fontSize)

% print('Figures\Figure - Orientation Distribution\FinalZ', '-depsc')
print('Figures\Figure - Orientation Distribution\FinalZ','-djpeg','-r100')
export_fig 'Figures\Figure - Orientation Distribution\FinalZ' -pdf -m8 -q101

%% Figure - Single Voltage Traces

% Params
figSize = [1 1 18 13]*2;
fontSize = 8*2;
marTop = 0.04;
marBot = 0.07;
marLef = 0.05;
marRig = 0.02;
spacHor = 0.07;
spacVer = 0.12;

Th = [12.09e-3 6.30e-3];

% Load 3D model image
% Model3D = imread('Figures\Figure - Single Voltage Traces\SimGeom_AllElectrodes.jpg', ...
%     'jpeg');

% Model3D = Model3D(1780:7060,320:6890,:);

% Load elec traces
figTrace1 = openfig('Figures\Figure - Single Voltage Traces\RotNeurite_OneElectrode.fig');
figTrace2 = openfig('Figures\Figure - Single Voltage Traces\RotNeurite_TwoElectrode.fig');
figTrace3 = openfig('Figures\Figure - Single Voltage Traces\RotNeurite_FourElectrode.fig');

% Combine

figSVT = figure('Units','centimeters','Position',figSize,'Name',...
    'Figure - Single Voltage Traces','Color','w');
set(figSVT, 'PaperPosition',figSize)

ax(1) = subaxis(2,2,1,'SpacingHoriz',0,'SpacingVert',0,'MarginRight',0,'MarginTop',0, ...
    'MarginLeft',0,'MarginBottom',0);
% image(Model3D)
axis equal
axis off

ax(2) = subaxis(2,2,2,'SpacingHoriz',spacHor,'SpacingVert',spacVer,'MarginRight',marRig,'MarginTop',marTop, ...
    'MarginLeft',marLef,'MarginBottom',marBot);
set(figTrace1.CurrentAxes.Children,'LineWidth',1.5)
copyobj(figTrace1.CurrentAxes.Children,ax(2))

ax(3) = subaxis(2,2,3,'SpacingHoriz',spacHor,'SpacingVert',spacVer,'MarginRight',marRig,'MarginTop',marTop, ...
    'MarginLeft',marLef,'MarginBottom',marBot);
set(figTrace2.CurrentAxes.Children,'LineWidth',1.5)
copyobj(figTrace2.CurrentAxes.Children,ax(3))

ax(4) = subaxis(2,2,4,'SpacingHoriz',spacHor,'SpacingVert',spacVer,'MarginRight',marRig,'MarginTop',marTop, ...
    'MarginLeft',marLef,'MarginBottom',marBot);
set(figTrace3.CurrentAxes.Children,'LineWidth',1.5)
copyobj(figTrace3.CurrentAxes.Children,ax(4))

letters = {'(a)','(b)','(c)','(d)'};
set(figSVT, 'CurrentAxes', ax(3))
text(-1140, 41.6, letters{1},'FontSize',fontSize);
set(ax(1),'fontsize',fontSize)
for i = 2:numel(ax)
    xlabel(ax(i),'Neurite axes (\mum)')
    ylabel(ax(i),'V_{m} (mV)')
    xlim(ax(i),[-1000 1000])
    ylim(ax(i),[-6 14])
    set(ax(i),'XTick',[-1000 0 1000])
    set(ax(i),'YTick',[-4 0 4 8 12])
    set(ax(i),'fontsize',fontSize)
    ax(i).Box = 'off';
    title(ax(i),letters{i},'Position',[-1110 14.5],'FontWeight','normal');
end

set(figSVT, 'CurrentAxes', ax(2))
text(-950,Th(1)*1e3+1,'V_{th,AOP}','FontSize',fontSize)
text(-950,Th(2)*1e3+1,'V_{th,AIS}','FontSize',fontSize)

close(figTrace1)
close(figTrace2)
close(figTrace3)

set(figSVT, 'CurrentAxes', ax(2))
hold on
rectangle('Position',[550 7.5 400 4*figSize(3)/figSize(4)],'FaceColor','w')
active = ['w' 'w' 'w' 'w' 'k' 'w' 'w' 'w' 'w'];
xPeriod = 400/3;
yPeriod = 4*figSize(3)/figSize(4)/3;
x0 = [1; 1; 1]*(550+(xPeriod/2:xPeriod:xPeriod*3));
y0 = (7.5+(yPeriod/2:yPeriod:yPeriod*3))'*[1 1 1];
ra = xPeriod/4;
rb = yPeriod/4;
h = ellipse(ra,rb,0,x0(:),y0(:),'k');
for i = 1:numel(h)
    x=get(h(i),'Xdata');
    y=get(h(i),'Ydata');
    patch(x,y,active(i));
end
hold off

set(figSVT, 'CurrentAxes', ax(3))
hold on
rectangle('Position',[550-xPeriod 7.5 xPeriod*4 4*figSize(3)/figSize(4)],'FaceColor','w')
active = ['w' 'w' 'w' 'w' 'k' 'w' 'w' 'k' 'w' 'w' 'w' 'w'];
xPeriod = 400/3;
yPeriod = 4*figSize(3)/figSize(4)/3;
x0 = [1; 1; 1]*(550-xPeriod+(xPeriod/2:xPeriod:xPeriod*4));
y0 = (7.5+(yPeriod/2:yPeriod:yPeriod*3))'*[1 1 1 1];
ra = xPeriod/4;
rb = yPeriod/4;
h = ellipse(ra,rb,0,x0(:),y0(:),'k');
for i = 1:numel(h)
    x=get(h(i),'Xdata');
    y=get(h(i),'Ydata');
    patch(x,y,active(i));
end
hold off

set(figSVT, 'CurrentAxes', ax(4))
hold on
rectangle('Position',[550-xPeriod 7.5 xPeriod*4 4*figSize(3)/figSize(4)],'FaceColor','w')
active = ['w' 'k' 'w' 'w' 'k' 'w' 'w' 'k' 'w' 'w' 'k' 'w'];
xPeriod = 400/3;
yPeriod = 4*figSize(3)/figSize(4)/3;
x0 = [1; 1; 1]*(550-xPeriod+(xPeriod/2:xPeriod:xPeriod*4));
y0 = (7.5+(yPeriod/2:yPeriod:yPeriod*3))'*[1 1 1 1];
ra = xPeriod/4;
rb = yPeriod/4;
h = ellipse(ra,rb,0,x0(:),y0(:),'k');
for i = 1:numel(h)
    x=get(h(i),'Xdata');
    y=get(h(i),'Ydata');
    patch(x,y,active(i));
end
hold off

print('Figures\Figure - Single Voltage Traces\Final.jpg','-djpeg','-r600')
export_fig 'Figures\Figure - Single Voltage Traces\Final' -eps -m8 -q101

%% Figure - Param Search (2 x 3)

% Params
figSize = [1 1 18 10]*2;
fontSize = 8*2;
marTop = 0.05;
marBot = 0.09;
marLef = 0.06;
marRig = 0.09;
spacHor = 0.07;
spacVer = 0.13;

Th = [12.09e-3 6.30e-3];

% Load subplots
fig1 = openfig('Figures\Figure - Param Search\ParamSearch_OneElectrode.fig');
fig2 = openfig('Figures\Figure - Param Search\ParamSearch_TwoElectrode.fig');
fig3 = openfig('Figures\Figure - Param Search\ParamSearch_FourElectrode.fig');
fig4 = openfig('Figures\Figure - Param Search\RotNeurite_FourElectrode_Y100D100.fig');
fig5 = openfig('Figures\Figure - Param Search\RotNeurite_FourElectrode_Y300D100.fig');
fig6 = openfig('Figures\Figure - Param Search\RotNeurite_FourElectrode_Y100D300.fig');

% Load plots for sensitivity analysis
fig7 = openfig('Figures\Figure - Param Search\ParamSearch_FourElectrode_ratio1.fig');
fig8 = openfig('Figures\Figure - Param Search\ParamSearch_FourElectrode_ratio1p5.fig');

% Conduct sensitivity analysis
CData_2p0 = fig3.CurrentAxes.Children(2).CData;
CData_1p5 = fig7.CurrentAxes.Children(2).CData;
CData_1p0 = fig8.CurrentAxes.Children(2).CData;
disp(['Proportion over 10% pref. act. for 2.0 is ', ...
    num2str(sum(CData_2p0(:) > 0.1) / numel(CData_2p0(:)))])
disp(['Proportion over 10% pref. act. for 1.5 is ', ...
    num2str(sum(CData_1p5(:) > 0.1) / numel(CData_1p5(:)))])
disp(['Proportion over 10% pref. act. for 1.0 is ', ...
    num2str(sum(CData_1p0(:) > 0.1) / numel(CData_1p0(:)))])

% Combine
figPS = figure('Units','centimeters','Position',figSize,'Name',...
    'Figure - Param Search','Color','w');
set(figPS, 'PaperPosition',figSize)

ax(1) = subaxis(2,3,1,'SpacingHoriz',spacHor,'SpacingVert',spacVer,'MarginRight',marRig,'MarginTop',marTop, ...
    'MarginLeft',marLef,'MarginBottom',marBot);
copyobj(fig1.CurrentAxes.Children,ax(1))
hold on
plot(200,200,'ko','MarkerSize',6,'MarkerFaceColor','w')
% text(200+8,200+20,'Fig. 3(b)','FontSize',fontSize*0.8,'Color','w')
text(450,465,'10%','FontSize',fontSize*0.8,'Color','k')
hold off

ax(2) = subaxis(2,3,2,'SpacingHoriz',spacHor,'SpacingVert',spacVer,'MarginRight',marRig,'MarginTop',marTop, ...
    'MarginLeft',marLef,'MarginBottom',marBot);
copyobj(fig2.CurrentAxes.Children,ax(2))
hold on
plot(200,200,'ko','MarkerSize',6,'MarkerFaceColor','w')
% text(200+8,200+20,'Fig. 3(c)','FontSize',fontSize*0.8,'Color','w')
text(450,310,'10%','FontSize',fontSize*0.8,'Color','k')
text(450,460,'40%','FontSize',fontSize*0.8,'Color','k')
hold off

ax(3) = subaxis(2,3,3,'SpacingHoriz',spacHor,'SpacingVert',spacVer,'MarginRight',marRig,'MarginTop',marTop, ...
    'MarginLeft',marLef,'MarginBottom',marBot);
copyobj(fig3.CurrentAxes.Children,ax(3))
origSize = get(ax(3),'Position');
axis tight
axis xy
cb = colorbar;
ylabel(cb,'Preferential AIS activation (%)')
cb.Ticks = 0:0.25:1;
cb.TickLabels = {'0', '25', '50', '75', '100'};
set(ax(3),'Position',origSize)
hold on
plot([100 100 300],[100 300 100],'ks','MarkerSize',6,'MarkerFaceColor','k')
text([100 100 300]+8,[100 300 100]+20,{'(d)','(e)','(f)'},'FontSize',fontSize*0.8)
plot(200,200,'ko','MarkerSize',6,'MarkerFaceColor','w')
% text(200+8,200+20,'Fig. 3(d)','FontSize',fontSize*0.8,'Color','k')
text(450,30,'10%','FontSize',fontSize*0.8,'Color','k')
text(450,65,'40%','FontSize',fontSize*0.8,'Color','k')
hold off

ax(4) = subaxis(2,3,4,'SpacingHoriz',spacHor,'SpacingVert',spacVer,'MarginRight',marRig,'MarginTop',marTop, ...
    'MarginLeft',marLef,'MarginBottom',marBot);
set(fig4.CurrentAxes.Children,'LineWidth',1.5)
copyobj(fig4.CurrentAxes.Children,ax(4))

ax(5) = subaxis(2,3,5,'SpacingHoriz',spacHor,'SpacingVert',spacVer,'MarginRight',marRig,'MarginTop',marTop, ...
    'MarginLeft',marLef,'MarginBottom',marBot);
set(fig5.CurrentAxes.Children,'LineWidth',1.5)
copyobj(fig5.CurrentAxes.Children,ax(5))

ax(6) = subaxis(2,3,6,'SpacingHoriz',spacHor,'SpacingVert',spacVer,'MarginRight',marRig,'MarginTop',marTop, ...
    'MarginLeft',marLef,'MarginBottom',marBot);
set(fig6.CurrentAxes.Children,'LineWidth',1.5)
copyobj(fig6.CurrentAxes.Children,ax(6))

axis(ax(1:3),'tight')
axis(ax(1:3),'xy')
xlabel(ax(2),'Pulse duration (\mus)')
ylabel(ax(1),'Electrode height (\mum)')
xlabel(ax(5),'Neurite axes (\mum)')
ylabel(ax(4),'V_{m} (mV)')

letters = {'(a)','(b)','(c)','(d)','(e)','(f)'};
for i = 1:3
    rectangle(ax(i),'Position',[5 5 500 500])
    set(ax(i),'fontsize',fontSize)
    title(ax(i),letters{i},'Position',[-62 525],'FontWeight','normal');
    colormap(ax(i),flipud(hot))
    caxis(ax(i),[0 1])
    for j = 1:numel(ax(i).Children)
        if strcmp(ax(i).Children(j).Type, 'contour')
            ax(i).Children(j).LineColor = 'k';
        end
    end
end
for i = 4:6
    xlim(ax(i),[-1000 1000])
    ylim(ax(i),[-6 14])
    set(ax(i),'XTick',[-1000 0 1000])
    set(ax(i),'YTick',[-4 0 4 8 12])
    set(ax(i),'fontsize',fontSize)
    ax(i).Box = 'off';
    title(ax(i),letters{i},'Position',[-1250 14.2],'FontWeight','normal');
end

set(figPS, 'CurrentAxes', ax(4))
text(-950,Th(1)*1e3+1.25,'V_{th,AOP}','FontSize',fontSize)
text(-950,Th(2)*1e3+1.25,'V_{th,AIS}','FontSize',fontSize)

xPeriod = 25;
yPeriod = xPeriod*1.05;
set(figPS, 'CurrentAxes', ax(1))
hold on
rectangle('Position',[445-xPeriod 20 xPeriod*3 yPeriod*3], ...
    'FaceColor','w','EdgeColor','k')
active = ['w' 'w' 'w' 'w' 'k' 'w' 'w' 'w' 'w'];
x0 = [1; 1; 1]*(445-xPeriod+(xPeriod/2:xPeriod:xPeriod*3));
y0 = (20+(yPeriod/2:yPeriod:yPeriod*3))'*[1 1 1];
ra = xPeriod/4;
rb = yPeriod/4;
h = ellipse(ra,rb,0,x0(:),y0(:),'k');
for i = 1:numel(h)
    x=get(h(i),'Xdata');
    y=get(h(i),'Ydata');
    patch(x,y,active(i));
end
hold off

set(figPS, 'CurrentAxes', ax(2))
hold on
rectangle('Position',[420-xPeriod 20 xPeriod*4 yPeriod*3], ...
    'FaceColor','w','EdgeColor','k')
active = ['w' 'w' 'w' 'w' 'k' 'w' 'w' 'k' 'w' 'w' 'w' 'w'];
x0 = [1; 1; 1]*(420-xPeriod+(xPeriod/2:xPeriod:xPeriod*4));
y0 = (20+(yPeriod/2:yPeriod:yPeriod*3))'*[1 1 1 1];
ra = xPeriod/4;
rb = yPeriod/4;
h = ellipse(ra,rb,0,x0(:),y0(:),'k');
for i = 1:numel(h)
    x=get(h(i),'Xdata');
    y=get(h(i),'Ydata');
    patch(x,y,active(i));
end
hold off

set(figPS, 'CurrentAxes', ax(3))
hold on
rectangle('Position',[420-xPeriod 415 xPeriod*4 yPeriod*3], ...
    'FaceColor','w','EdgeColor','k')
active = ['w' 'k' 'w' 'w' 'k' 'w' 'w' 'k' 'w' 'w' 'k' 'w'];
x0 = [1; 1; 1]*(420-xPeriod+(xPeriod/2:xPeriod:xPeriod*4));
y0 = (415+(yPeriod/2:yPeriod:yPeriod*3))'*[1 1 1 1];
ra = xPeriod/4;
rb = yPeriod/4;
h = ellipse(ra,rb,0,x0(:),y0(:),'k');
for i = 1:numel(h)
    x=get(h(i),'Xdata');
    y=get(h(i),'Ydata');
    patch(x,y,active(i));
end
hold off


xPeriod = 400/4;
yPeriod = 1.125;
for j = 4:6
    set(figPS, 'CurrentAxes', ax(j))
    hold on
    rectangle('Position',[650-xPeriod 9.5 xPeriod*4 yPeriod*3], ...
        'FaceColor','w','EdgeColor','k')
    active = ['w' 'k' 'w' 'w' 'k' 'w' 'w' 'k' 'w' 'w' 'k' 'w'];
    x0 = [1; 1; 1]*(650-xPeriod+(xPeriod/2:xPeriod:xPeriod*4));
    y0 = (9.5+(yPeriod/2:yPeriod:yPeriod*3))'*[1 1 1 1];
    ra = xPeriod/4;
    rb = yPeriod/4;
    h = ellipse(ra,rb,0,x0(:),y0(:),'k');
    for i = 1:numel(h)
        x=get(h(i),'Xdata');
        y=get(h(i),'Ydata');
        patch(x,y,active(i));
    end
    hold off
end

for i = 1:3
    ax_tmp = axes('Position',get(ax(i),'position'), ...
        'Xlim',get(ax(i),'Xlim'),'Ylim',get(ax(i),'Ylim'),'color','none');
    ax_tmp.XTick = [0 ax(i).XTick]+5;
    ax_tmp.XTickLabel = {};
    ax(i).XTickLabel = ['0'; ax(i).XTickLabel];
    ax(i).XTick = ax_tmp.XTick;
    ax_tmp.YTick = [0 ax(i).YTick]+5;
    ax_tmp.YTickLabel = {};
    ax(i).YTickLabel = ['0'; ax(i).YTickLabel];
    ax(i).YTick = ax_tmp.YTick;
end

close(fig1)
close(fig2)
close(fig3)
close(fig4)
close(fig5)
close(fig6)
close(fig7)
close(fig8)

print('Figures\Figure - Param Search\Final.jpg','-djpeg','-r100')
export_fig 'Figures\Figure - Param Search\Final' -eps -m8 -q101
% saveas(figSVT, 'Figures\Figure - Single Voltage Traces\Final_2x3','fig')

%% Figure - Current vs Activation

% Params
figSize = [1 1 18 18]*2;
fontSize = 8*2;
marTop = 0.03;
marBot = 0.06;
marLef = 0.06;
marRig = 0.027;
spacHor = 0.05;
spacVer = 0.09;
elecArea = 0.00007853981; % cm^2

axLim = 700;

% Load subplots
fig1 = openfig('Figures\Figure - Current vs Activation\CurrentVsActAndRadius.fig');
fig2 = openfig('Figures\Figure - Current vs Activation\InsetCurrentVsActAndRadius.fig');
fig3 = openfig('Figures\Figure - Current vs Activation\RotNeurite_OneElectrode_Plane100um.fig');
fig4 = openfig('Figures\Figure - Current vs Activation\RotNeurite_TwoElectrode_Plane100um.fig');
fig5 = openfig('Figures\Figure - Current vs Activation\RotNeurite_FourElectrode_Plane100um.fig');

% Combine
figCVA = figure('Units','centimeters','Position',figSize,'Name',...
    'Figure - Current vs Activation','Color','w');
set(figCVA, 'PaperPosition',figSize)

ax(1) = subaxis(3,2,1,'SpacingHoriz',spacHor,'SpacingVert',spacVer,'MarginRight',marRig,'MarginTop',marTop, ...
    'MarginLeft',marLef,'MarginBottom',marBot);
copyobj(fig1.Children(3).Children,ax(1))
copyobj(fig1.Children(2).Children,ax(1))
for i = [1 4:9]
    ax(1).Children(i).XData = ax(1).Children(i).XData/(4*elecArea)*1e-6*0.0002*1e3;
end
for i = [2 10:15]
    ax(1).Children(i).XData = ax(1).Children(i).XData/(2*elecArea)*1e-6*0.0002*1e3;
end
for i = [3 16:21]
    ax(1).Children(i).XData = ax(1).Children(i).XData/(elecArea)*1e-6*0.0002*1e3;
end
yl1 = ylabel('GCL activation level (%)');
xlabel('Stimulus charge density per phase (mC/cm^{2})')
xlim([0 3])
set(ax(1),'YTick',0:25:100)
legendStr = ['1 electrode'; '2 electrode'; '4 electrode'];
l = legend(ax(1).Children(18:-6:1),legendStr,'Location','SouthEast');
set(l,'fontsize',fontSize)
set(l,'Box','off')

ax(2) = subaxis(3,2,2,'SpacingHoriz',spacHor,'SpacingVert',spacVer,'MarginRight',marRig,'MarginTop',marTop, ...
    'MarginLeft',marLef,'MarginBottom',marBot);
copyobj(fig3.CurrentAxes.Children,ax(2))

ax(3) = subaxis(3,2,3,'SpacingHoriz',spacHor,'SpacingVert',spacVer,'MarginRight',marRig,'MarginTop',marTop, ...
    'MarginLeft',marLef,'MarginBottom',marBot);
copyobj(fig1.Children(1).Children,ax(3))
for i = [1 4:9]
    ax(3).Children(i).XData = ax(3).Children(i).XData/(4*elecArea)*1e-6*0.0002*1e3;
end
for i = [2 10:15]
    ax(3).Children(i).XData = ax(3).Children(i).XData/(2*elecArea)*1e-6*0.0002*1e3;
end
for i = [3 16:21]
    ax(3).Children(i).XData = ax(3).Children(i).XData/(elecArea)*1e-6*0.0002*1e3;
end
yl2 = ylabel('Activation radius (\mum)');
xlabel('Stimulus charge density per phase (mC/cm^{2})')
xlim([0 3])
ylim([0 800])
set(ax(3),'YTick',0:200:800)

ax(4) = subaxis(3,2,4,'SpacingHoriz',spacHor,'SpacingVert',spacVer,'MarginRight',marRig,'MarginTop',marTop, ...
    'MarginLeft',marLef,'MarginBottom',marBot);
copyobj(fig4.CurrentAxes.Children,ax(4))

ax(5) = subaxis(3,2,5,'SpacingHoriz',spacHor,'SpacingVert',spacVer,'MarginRight',marRig,'MarginTop',marTop, ...
    'MarginLeft',marLef,'MarginBottom',marBot);
copyobj(fig2.CurrentAxes.Children,ax(5))
xlabel('Activation radius (\mum)')
yl3 = ylabel('GCL activation level (%)');
set(ax(5),'YTick',0:25:100)

ax(6) = subaxis(3,2,6,'SpacingHoriz',spacHor,'SpacingVert',spacVer,'MarginRight',marRig,'MarginTop',marTop, ...
    'MarginLeft',marLef,'MarginBottom',marBot);
copyobj(fig5.CurrentAxes.Children,ax(6))
line([400 600],[-410 -410],'LineWidth',4,'Color','k')
text(380,-300,'200\mum','FontSize',fontSize,'Color','k')
text(-630,-300,'AOP orientation','FontSize',fontSize,'Color',[1 0.5 0])
annotation('arrow',[0.74 0.81],[0.079 0.079], 'Color', [1 0.5 0], ...
    'LineWidth',4)

for i = 1:2:5
    set(ax(i),'fontsize',fontSize)
end

for i = 2:2:6
    set(figCVA, 'CurrentAxes', ax(i))
    axis xy
    colormap(flipud(hot))
    caxis([0 1])
    xlim([-axLim axLim])
    ylim([-axLim axLim]*2/3)
    rectangle('Position',[-axLim -axLim*2/3 2*axLim 2*axLim*2/3])
    set(ax(i),'XTickLabel',[])
    set(ax(i),'YTickLabel',[])
    set(ax(i),'fontsize',fontSize)
    line([0 0],[-axLim axLim],'Color',[0 0.4470 0.7410],'LineStyle','--')
end
set(ax,'Box','off')

axPos = zeros(numel(ax),4);
for i = 1:numel(ax)
    axPos(i,:) = ax(i).Position;
end

yAdj = (axPos(1,2)-axPos(3,2)-axPos(3,4))/8;
xAdj = axPos(1,3)*0.4;
yAdj2 = (axPos(2,2)-axPos(4,2)-axPos(4,4))/11;
ax(1).Position = axPos(1,:) + [0 -yAdj xAdj 0];
ax(3).Position = axPos(3,:) + [0 -yAdj xAdj 0];
ax(5).Position = axPos(5,:) + [0 -yAdj xAdj 0];
ax(2).Position = axPos(2,:) + [xAdj*0.95 -3*yAdj2 -xAdj*0.9 2*yAdj2];
ax(4).Position = axPos(4,:) + [xAdj*0.95 -1*yAdj2 -xAdj*0.9 2*yAdj2];

cb = colorbar('peer',ax(6),'SouthOutside');
ylabel(cb,'GCL activation level (%)')
set(cb,'XTick',0:0.25:1)
set(cb,'XTickLabel',0:25:100)
set(cb,'fontsize',fontSize)

ax(6).Position = axPos(6,:) + [xAdj*0.95 1*yAdj2 -xAdj*0.9 2*yAdj2];
cb.Position(4) = cb.Position(4)*0.85;
cb.Position(2) = cb.Position(2) + cb.Position(4)*3/7;

pause(1)
yl2.Position(1) = yl2.Position(1)*0.95;
yl1.Position(1) = yl2.Position(1);
yl3.Position(1) = yl3.Position(1)*1.16;


for j = 1:2:5
    for i = numel(ax(j).Children):-1:1
        if strcmp(ax(j).Children(i).Marker,'square')
            ax(j).Children(i).MarkerSize = 6;
        end
    end
end
for i = 2:2:6
    ax(i).Position = ax(i).Position.*[1 1 1 2/3];
end

for i = 1:numel(ax(5).Children)
    if strcmp(ax(5).Children(i).Type,'line') && ax(5).Children(i).LineWidth == 0.5
        ax(5).Children(i).LineStyle = ':';
        ax(5).Children(i).LineWidth = 1.5;
    end
end

set(figCVA, 'CurrentAxes', ax(1))
text([ax(1).Children(3).XData ax(1).Children(2).XData ax(1).Children(1).XData]-0.4*3/14, ...
    [ax(1).Children(3).YData ax(1).Children(2).YData ax(1).Children(1).YData]+5, ...
    {'(d)','(e)','(f)'},'FontSize',fontSize*0.8)
set(figCVA, 'CurrentAxes', ax(3))
text([ax(3).Children(3).XData ax(3).Children(2).XData ax(3).Children(1).XData]-0.4*3/14, ...
    [ax(3).Children(3).YData ax(3).Children(2).YData ax(3).Children(1).YData]+30, ...
    {'(d)','(e)','(f)'},'FontSize',fontSize*0.8)
set(figCVA, 'CurrentAxes', ax(5))
text([ax(5).Children(7).XData ax(5).Children(4).XData ax(5).Children(1).XData]-24, ...
    [ax(5).Children(7).YData ax(5).Children(4).YData ax(5).Children(1).YData]+4, ...
    {'(d)','(e)','(f)'},'FontSize',fontSize*0.8)

title(ax(1),'(a)','Position',[-0.185 107],'FontWeight','normal');
title(ax(3),'(b)','Position',[-0.185 860],'FontWeight','normal');
title(ax(5),'(c)','Position',[-50 107],'FontWeight','normal');
title(ax(2),'(d)','Position',[-axLim*1.0 axLim*1.42],'FontWeight','normal');
title(ax(4),'(e)','Position',[-axLim*1.0 axLim*1.42],'FontWeight','normal');
title(ax(6),'(f)','Position',[-axLim*1.0 axLim*1.42],'FontWeight','normal');

set(figCVA, 'CurrentAxes', ax(1))
line([0.28 0.28 0.7; 0.28 0.7 0.7],[95 100 95; 100 100 100],'Color','k')
line([0.93 0.93 2.71; 0.93 2.71 2.71],[95 100 95; 100 100 100],'Color','k')
text(0.2,106,['d_{ER} = 100\mum (' char(9642) ')'],'FontSize',fontSize)
text(1.45,106,['d_{ER} = 300\mum (' char(9653) ')'],'FontSize',fontSize)

set(figCVA, 'CurrentAxes', ax(3))
text(0.3692,820,char(9642),'FontSize',fontSize)
text(0.8403,820,char(9642),'FontSize',fontSize)
text(1.757,820,char(9642),'FontSize',fontSize)
text(1.216,820,char(9653),'FontSize',fontSize)
text(2.546,820,char(9653),'FontSize',fontSize)
text(3.03,515,char(9653),'FontSize',fontSize)

set(figCVA, 'CurrentAxes', ax(5))
count = 1;
for i = numel(ax(5).Children):-1:1
    if strcmp(ax(5).Children(i).Type,'line') && strcmp(ax(5).Children(i).LineStyle,'--') && ...
            ~isequal(ax(5).Children(i).Color,[0 0 0])
        if mod(count,2) == 1
            text(ax(5).Children(i).XData(end)-5,104,char(9642),'FontSize',fontSize)
        else
            text(ax(5).Children(i).XData(end)-5,103.7,char(9653),'FontSize',fontSize)
        end
        count = count + 1;
    end
end

set(figCVA, 'CurrentAxes', ax(2))
ellipse(50,50,0,0,0,'k');
set(figCVA, 'CurrentAxes', ax(4))
ellipse(50,50,0,[0 0],[-100 100],'k');
set(figCVA, 'CurrentAxes', ax(6))
ellipse(50,50,0,[0 0 0 0],[-300 -100 100 300],'k');

for j = 2:2:6
    for i = 1:numel(ax(j).Children)
        if strcmp(ax(j).Children(i).Type,'image')
            ax(j).Children(i).CData = ax(j).Children(i).CData';
        elseif strcmp(ax(j).Children(i).Type,'contour')
            ax(j).Children(i).ZData = ax(j).Children(i).ZData';
            ax(j).Children(i).LineColor = 'k';
            ax(j).Children(i).LineStyle = 'none';
        elseif strcmp(ax(j).Children(i).Type,'line') && ax(j).Children(i).LineWidth < 3
            XData = ax(j).Children(i).YData;
            YData = ax(j).Children(i).XData;
            ax(j).Children(i).XData = XData;
            ax(j).Children(i).YData = YData;
        end
    end
end

close(fig1)
close(fig2)
close(fig3)
close(fig4)
close(fig5)

for i = 2:2:6
    axes('Position',[ax(i).Position(1) ax(i).Position(2)+ax(i).Position(4) ...
        ax(i).Position(3) ax(i).Position(4)/2],'Box','on');
    for j = 1:numel(ax(i).Children)
        if strcmp(ax(i).Children(j).Type,'image')
            yD = ax(i).Children(j).YData;
            cD = ax(i).Children(j).CData(201,:);
            cD(cD > 0) = smooth(cD(cD > 0),3);
            cD = cD/max(cD)/2;
            plot(yD(abs(yD) < 700), cD(abs(yD) < 700), 'LineWidth', 1.5)
            yl = ylabel('%','FontSize',fontSize);
            yl.Rotation = 0;
            yl.Position = yl.Position + [150 -0.05 0];
            xlim([-700 700])
            set(gca,'XTick',[])
            set(gca,'YTick',[0 0.5],'YTickLabel',{'0', '50'})
            set(gca,'FontSize',fontSize)
        end
    end
end
% set(gcf, 'PaperPosition', get(gcf, 'Position'))
set(gcf, 'PaperSize', get(gcf, 'PaperSize')*2)


print('Figures\Figure - Current vs Activation\Final.jpg','-djpeg','-r100')
print('Figures\Figure - Current vs Activation\Final', '-dpdf')

%% Figure - Non-ideal conditions 2

% Params
figSize = [1 1 18 18]*2;
fontSize = 8*2;
marTop = 0.02;
marBot = 0.01;
marLef = 0.07;
marRig = 0.02;
spacHor = 0.08;
spacVer = 0.07;

axLim = 1000;

Th = [12.09e-3 6.30e-3];

% Load subplots
fig1 = openfig('Figures\Figure - Non-ideal Conditions\RotNeurite_EightElectrode.fig');
fig2 = openfig('Figures\Figure - Non-ideal Conditions\RotNeurite_EightElectrode_Plane.fig');
fig3 = openfig('Figures\Figure - Non-ideal Conditions\RotNeurite_225DegAOP.fig');
fig4 = openfig('Figures\Figure - Non-ideal Conditions\RotNeurite_225DegAOP_Plane.fig');
fig5 = openfig('Figures\Figure - Non-ideal Conditions\InsetCurrentVsActAndRadius_EightElectrode.fig');
fig6 = openfig('Figures\Figure - Non-ideal Conditions\InsetCurrentVsActAndRadius_225DegAOP.fig');

fig0 = openfig('Figures\Figure - Current vs Activation\RotNeurite_FourElectrode_Plane100um.fig');

% Combine
figNIC = figure('Units','centimeters','Position',figSize,'Name',...
    'Figure - Non-ideal conditions','Color','w');
set(figNIC, 'PaperPosition',figSize)

ax(1) = subaxis(3,2,1,'SpacingHoriz',spacHor,'SpacingVert',spacVer,'MarginRight',marRig,'MarginTop',marTop, ...
    'MarginLeft',marLef,'MarginBottom',marBot);
set(fig1.CurrentAxes.Children,'LineWidth',1.5)
copyobj(fig1.CurrentAxes.Children,ax(1))
ylabel(ax(1),'V_{m} (mV)')

ax(2) = subaxis(3,2,2,'SpacingHoriz',spacHor,'SpacingVert',spacVer,'MarginRight',marRig,'MarginTop',marTop, ...
    'MarginLeft',marLef,'MarginBottom',marBot);
set(fig3.CurrentAxes.Children,'LineWidth',1.5)
copyobj(fig3.CurrentAxes.Children,ax(2))

ax(3) = subaxis(3,2,3,'SpacingHoriz',spacHor,'SpacingVert',spacVer,'MarginRight',marRig,'MarginTop',marTop, ...
    'MarginLeft',marLef,'MarginBottom',marBot);
copyobj(fig5.CurrentAxes.Children,ax(3))
ylabel('GCL activation level (%)');
legendStr = {'Non-ideal','Ideal'};
l = legend(ax(3).Children([10 12]),legendStr,'Location','NorthWest');
set(l,'fontsize',fontSize)
set(l,'Box','off')

ax(4) = subaxis(3,2,4,'SpacingHoriz',spacHor,'SpacingVert',spacVer,'MarginRight',marRig,'MarginTop',marTop, ...
    'MarginLeft',marLef,'MarginBottom',marBot);
copyobj(fig6.CurrentAxes.Children,ax(4))

ax(5) = subaxis(3,2,5,'SpacingHoriz',spacHor,'SpacingVert',spacVer,'MarginRight',marRig,'MarginTop',marTop, ...
    'MarginLeft',marLef,'MarginBottom',marBot);
copyobj(fig2.CurrentAxes.Children,ax(5))

ax(6) = subaxis(3,2,6,'SpacingHoriz',spacHor,'SpacingVert',spacVer,'MarginRight',marRig,'MarginTop',marTop, ...
    'MarginLeft',marLef,'MarginBottom',marBot);
copyobj(fig4.CurrentAxes.Children,ax(6))
line([650 850],[-360 -360],'LineWidth',4,'Color','k')
text(630,-295,'200\mum','FontSize',fontSize,'Color','k')

active = ['w' 'k' 'k' 'w' 'w' 'k' 'k' 'w' 'w' 'k' 'k' 'w' 'w' 'k' 'k' 'w'; ...
    'w' 'k' 'w' 'w' 'w' 'w' 'k' 'w' 'w' 'w' 'k' 'w' 'w' 'w' 'w' 'k'];
rot = [0 pi/8];
target = [0 0; 0 100];
x0 = repmat(-300:200:300,4,1);
y0 = x0';

for i = 1:2
    xlim(ax(i),[-axLim axLim])
    ylim(ax(i),[-6 16])
    set(ax(i),'XTick',[-axLim 0 axLim])
    set(ax(i),'YTick',[-4 0 4 8 12 16])
    set(ax(i),'fontsize',fontSize)
    xlabel(ax(i),'Neurite axes (\mum)')
    ax(i).Box = 'off';
    
    axes('Position',[ax(i).Position(1)+ax(i).Position(3)*0.7 ax(i).Position(2)+ax(i).Position(4)*0.65 ...
        ax(i).Position(3)*0.35 ax(i).Position(4)*0.4],'Box','on');
    
    hold on
    line([-400 400], ...
        [-400*tan(rot(i)) 400*tan(rot(i))]+target(i,2), ...
        'Color',[1 0.5 0],'LineWidth',1.5)
    plot(target(i,1),target(i,2),'ro','MarkerFaceColor','r')
    
    x0 = x0(:);
    y0 = y0(:);
    scatter(x0(active(i,:) == 'k'),y0(active(i,:) == 'k'),130,'k','filled','MarkerEdgeColor','k')
    scatter(x0(active(i,:) == 'w'),y0(active(i,:) == 'w'),130,'w','filled','MarkerEdgeColor','k')
    axis equal
    xlim([-400 400])
    ylim([-400 400])
    xticks([])
    yticks([])
end

set(figNIC, 'CurrentAxes', ax(1))
text(-950,Th(1)*1e3+1.3,'V_{th,AOP}','FontSize',fontSize)
text(-950,Th(2)*1e3+1.3,'V_{th,AIS}','FontSize',fontSize)

for i = 3:4
    set(ax(i),'YTick',0:25:100)
    set(ax(i),'fontsize',fontSize)
    xlabel(ax(i),'Activation radius (\mum)')
    ax(i).Box = 'off';
    
    for j = numel(ax(i).Children):-1:1
        if strcmp(ax(i).Children(j).Marker,'square')
            ax(i).Children(j).MarkerSize = 6;
        end
        if isequal(round(ax(i).Children(j).Color,4),[0.9153 0.2816 0.2878])
            ax(i).Children(j).Color = [0 0.4470 0.7410];
        end
        if isequal(round(ax(i).Children(j).Color,4),[0.3467 0.5360 0.6907])
            ax(i).Children(j).Color = [1 0 0];
        end
    end
end

yD0 = fig0.CurrentAxes.Children(2).YData;
cD0 = fig0.CurrentAxes.Children(2).CData(:,201);
cD0(cD0 > 0) = smooth(cD0(cD0 > 0),3);
cD0 = cD0/max(cD0)/2;
for i = 5:6
    set(figNIC, 'CurrentAxes', ax(i))
    axis xy
    colormap(flipud(hot))
    caxis([0 1])
    xlim([-axLim axLim])
    ylim([-axLim axLim]*3.7/8)
    rectangle('Position',[-axLim -axLim*3.7/8 2*axLim 2*axLim*3.7/8])
    set(ax(i),'XTickLabel',[])
    set(ax(i),'YTickLabel',[])
    set(ax(i),'fontsize',fontSize)
    line([-axLim*tan(rot(i-4)) axLim*tan(rot(i-4))], ...
        [-axLim axLim], ...
        'Color',[0 0.4470 0.7410],'LineStyle','--')
    
    x0 = x0(:);
    y0 = y0(:);
    usedElec = active(i-4,:) == 'k';
    hold on
    h = ellipse(50,50,0,y0(usedElec)-target(i-4,2),x0(usedElec)-target(i-4,1),'k');
    hold off
    
    ax(i).Position = ax(i).Position.*[1 1.01 1 2/3];
    axes('Position',[ax(i).Position(1) ax(i).Position(2)+ax(i).Position(4) ...
        ax(i).Position(3) ax(i).Position(4)/2],'Box','on');
    for j = 1:numel(ax(i).Children)
        if strcmp(ax(i).Children(j).Type,'image')
            yD = ax(i).Children(j).YData;
            cD = ax(i).Children(j).CData(:,201);
            cD(cD > 0) = smooth(cD(cD > 0),3);
            cD = cD/max(cD)/2;
            plot(yD0, cD0, 'Color', [1 0 0], 'LineWidth', 1.5)
            hold on
            plot(yD, cD, 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5)
            hold off
            yl(i-4) = ylabel('%','FontSize',fontSize);
            xlim([-axLim axLim])
            ylim([0 0.5])
            set(gca,'XTick',[])
            set(gca,'YTick',[0 0.5],'YTickLabel',{'0', '50'})
            set(gca,'FontSize',fontSize)
        end
    end
    
    for j = 1:numel(ax(i).Children)
        if strcmp(ax(i).Children(j).Type,'image')
            ax(i).Children(j).CData = imrotate(ax(i).Children(j).CData', ...
                -rot(i-4)*180/pi,'bilinear','crop');
        elseif strcmp(ax(i).Children(j).Type,'contour')
            ax(i).Children(j).ZData = ax(i).Children(j).ZData';
            ax(i).Children(j).LineColor = 'k';
            ax(i).Children(j).LineStyle = 'none';
        elseif strcmp(ax(i).Children(j).Type,'line') && ax(i).Children(j).LineWidth < 3
            XData = ax(i).Children(j).YData;
            YData = ax(i).Children(j).XData;
            ax(i).Children(j).XData = XData;
            ax(i).Children(j).YData = YData;
        end
    end
end

for i = 1:2
    yl(i).Position = [-1.08e3 0.25 -1];
end

labels = ['(a)';'(c)';'(e)';'(b)';'(d)';'(f)'];

annotation('textbox',[0.01 0.905 0.1 0.1],'String',labels(1,:), ...
    'FontSize',fontSize,'EdgeColor','w')
annotation('textbox',[0.01 0.56 0.1 0.1],'String',labels(2,:), ...
    'FontSize',fontSize,'EdgeColor','w')
annotation('textbox',[0.01 0.217 0.1 0.1],'String',labels(3,:), ...
    'FontSize',fontSize,'EdgeColor','w')
annotation('textbox',[0.51 0.905 0.1 0.1],'String',labels(4,:), ...
    'FontSize',fontSize,'EdgeColor','w')
annotation('textbox',[0.51 0.56 0.1 0.1],'String',labels(5,:), ...
    'FontSize',fontSize,'EdgeColor','w')
annotation('textbox',[0.51 0.217 0.1 0.1],'String',labels(6,:), ...
    'FontSize',fontSize,'EdgeColor','w')

close(fig1)
close(fig2)
close(fig3)
close(fig4)
close(fig5)
close(fig6)
close(fig0)

set(gcf, 'PaperSize', get(gcf, 'PaperSize')*2)
print('Figures\Figure - Non-ideal Conditions\Final.jpg','-djpeg','-r100')
print('Figures\Figure - Non-ideal Conditions\Final', '-dpdf')

%% Figure - Non-ideal conditions supplementary

% Params
figSize = [1 1 18 18]*2;
fontSize = 8*2;
marTop = 0.02;
marBot = 0.01;
marLef = 0.07;
marRig = 0.02;
spacHor = 0.08;
spacVer = 0.07;

axLim = 1000;

Th = [12.09e-3 6.30e-3];

% Load subplots
fig1 = openfig('Figures\Figure - Non-ideal Conditions\RotNeurite_45DegOffsetAOP.fig');
fig2 = openfig('Figures\Figure - Non-ideal Conditions\RotNeurite_45DegOffsetAOP_Plane.fig');
fig3 = openfig('Figures\Figure - Non-ideal Conditions\RotNeurite_OffsetAOP.fig');
fig4 = openfig('Figures\Figure - Non-ideal Conditions\RotNeurite_OffsetAOP_Plane.fig');
fig5 = openfig('Figures\Figure - Non-ideal Conditions\InsetCurrentVsActAndRadius_45DegOffsetAOP.fig');
fig6 = openfig('Figures\Figure - Non-ideal Conditions\InsetCurrentVsActAndRadius_OffsetAOP.fig');

fig0 = openfig('Figures\Figure - Current vs Activation\RotNeurite_FourElectrode_Plane100um.fig');

% Combine
figNIC = figure('Units','centimeters','Position',figSize,'Name',...
    'Figure - Non-ideal conditions','Color','w');
set(figNIC, 'PaperPosition',figSize)

ax(1) = subaxis(3,2,1,'SpacingHoriz',spacHor,'SpacingVert',spacVer,'MarginRight',marRig,'MarginTop',marTop, ...
    'MarginLeft',marLef,'MarginBottom',marBot);
set(fig1.CurrentAxes.Children,'LineWidth',1.5)
copyobj(fig1.CurrentAxes.Children,ax(1))
ylabel(ax(1),'V_{m} (mV)')

ax(2) = subaxis(3,2,2,'SpacingHoriz',spacHor,'SpacingVert',spacVer,'MarginRight',marRig,'MarginTop',marTop, ...
    'MarginLeft',marLef,'MarginBottom',marBot);
set(fig3.CurrentAxes.Children,'LineWidth',1.5)
copyobj(fig3.CurrentAxes.Children,ax(2))

ax(3) = subaxis(3,2,3,'SpacingHoriz',spacHor,'SpacingVert',spacVer,'MarginRight',marRig,'MarginTop',marTop, ...
    'MarginLeft',marLef,'MarginBottom',marBot);
copyobj(fig5.CurrentAxes.Children,ax(3))
ylabel('GCL activation level (%)');
legendStr = {'Non-ideal','Ideal'};
l = legend(ax(3).Children([10 12]),legendStr,'Location','NorthWest');
set(l,'fontsize',fontSize)
set(l,'Box','off')

ax(4) = subaxis(3,2,4,'SpacingHoriz',spacHor,'SpacingVert',spacVer,'MarginRight',marRig,'MarginTop',marTop, ...
    'MarginLeft',marLef,'MarginBottom',marBot);
copyobj(fig6.CurrentAxes.Children,ax(4))

ax(5) = subaxis(3,2,5,'SpacingHoriz',spacHor,'SpacingVert',spacVer,'MarginRight',marRig,'MarginTop',marTop, ...
    'MarginLeft',marLef,'MarginBottom',marBot);
copyobj(fig2.CurrentAxes.Children,ax(5))

ax(6) = subaxis(3,2,6,'SpacingHoriz',spacHor,'SpacingVert',spacVer,'MarginRight',marRig,'MarginTop',marTop, ...
    'MarginLeft',marLef,'MarginBottom',marBot);
copyobj(fig4.CurrentAxes.Children,ax(6))
line([650 850],[-360 -360],'LineWidth',4,'Color','k')
text(630,-295,'200\mum','FontSize',fontSize,'Color','k')

active = ['w' 'k' 'w' 'w' 'w' 'k' 'k' 'w' 'w' 'w' 'k' 'k' 'w' 'w' 'w' 'k'; ...
    'w' 'g' 'k' 'w' 'w' 'g' 'k' 'w' 'w' 'g' 'k' 'w' 'w' 'g' 'k' 'w'];
rot = [pi/4 0];
target = [0 100; 0 50];
x0 = repmat(-300:200:300,4,1);
y0 = x0';

for i = 1:2
    xlim(ax(i),[-axLim axLim])
    ylim(ax(i),[-6 16])
    set(ax(i),'XTick',[-axLim 0 axLim])
    set(ax(i),'YTick',[-4 0 4 8 12 16])
    set(ax(i),'fontsize',fontSize)
    xlabel(ax(i),'Neurite axes (\mum)')
    ax(i).Box = 'off';
    
    axes('Position',[ax(i).Position(1)+ax(i).Position(3)*0.7 ax(i).Position(2)+ax(i).Position(4)*0.65 ...
        ax(i).Position(3)*0.35 ax(i).Position(4)*0.4],'Box','on');
    
    hold on
    line([-400 400], ...
        [-400*tan(rot(i)) 400*tan(rot(i))]+target(i,2), ...
        'Color',[1 0.5 0],'LineWidth',1.5)
    plot(target(i,1),target(i,2),'ro','MarkerFaceColor','r')
    
    x0 = x0(:);
    y0 = y0(:);
    scatter(x0(active(i,:) == 'k'),y0(active(i,:) == 'k'),130,'k','filled','MarkerEdgeColor','k')
    scatter(x0(active(i,:) == 'w'),y0(active(i,:) == 'w'),130,'w','filled','MarkerEdgeColor','k')
    scatter(x0(active(i,:) == 'g'),y0(active(i,:) == 'g'),130,[0.8 0.8 0.8],'filled','MarkerEdgeColor','k')
    axis equal
    xlim([-400 400])
    ylim([-400 400])
    xticks([])
    yticks([])
end

set(figNIC, 'CurrentAxes', ax(1))
text(-950,Th(1)*1e3+1.3,'V_{th,AOP}','FontSize',fontSize)
text(-950,Th(2)*1e3+1.3,'V_{th,AIS}','FontSize',fontSize)

for i = 3:4
    set(ax(i),'YTick',0:25:100)
    set(ax(i),'fontsize',fontSize)
    xlabel(ax(i),'Activation radius (\mum)')
    ax(i).Box = 'off';
    
    for j = numel(ax(i).Children):-1:1
        if strcmp(ax(i).Children(j).Marker,'square')
            ax(i).Children(j).MarkerSize = 6;
        end
        if isequal(round(ax(i).Children(j).Color,4),[0.9153 0.2816 0.2878])
            ax(i).Children(j).Color = [0 0.4470 0.7410];
        end
        if isequal(round(ax(i).Children(j).Color,4),[0.3467 0.5360 0.6907])
            ax(i).Children(j).Color = [1 0 0];
        end
    end
end

yD0 = fig0.CurrentAxes.Children(2).YData;
cD0 = fig0.CurrentAxes.Children(2).CData(:,201);
cD0(cD0 > 0) = smooth(cD0(cD0 > 0),3);
cD0 = cD0/max(cD0)/2;
for i = 5:6
    set(figNIC, 'CurrentAxes', ax(i))
    axis xy
    colormap(flipud(hot))
    caxis([0 1])
    xlim([-axLim axLim])
    ylim([-axLim axLim]*3.7/8)
    rectangle('Position',[-axLim -axLim*3.7/8 2*axLim 2*axLim*3.7/8])
    set(ax(i),'XTickLabel',[])
    set(ax(i),'YTickLabel',[])
    set(ax(i),'fontsize',fontSize)
    line([-axLim*tan(rot(i-4)) axLim*tan(rot(i-4))], ...
        [-axLim axLim], ...
        'Color',[0 0.4470 0.7410],'LineStyle','--')
    
    x0 = x0(:);
    y0 = y0(:);
    usedElec = active(i-4,:) == 'k' | active(i-4,:) == 'g';
    hold on
    h = ellipse(50,50,0,y0(usedElec)-target(i-4,2),x0(usedElec)-target(i-4,1),'k');
    hold off
    
    ax(i).Position = ax(i).Position.*[1 1.01 1 2/3];
    axes('Position',[ax(i).Position(1) ax(i).Position(2)+ax(i).Position(4) ...
        ax(i).Position(3) ax(i).Position(4)/2],'Box','on');
    for j = 1:numel(ax(i).Children)
        if strcmp(ax(i).Children(j).Type,'image')
            yD = ax(i).Children(j).YData;
            cD = ax(i).Children(j).CData(:,201);
            cD(cD > 0) = smooth(cD(cD > 0),3);
            cD = cD/max(cD)/2;
            plot(yD0, cD0, 'Color', [1 0 0], 'LineWidth', 1.5)
            hold on
            plot(yD, cD, 'Color', [0 0.4470 0.7410], 'LineWidth', 1.5)
            hold off
            yl(i-4) = ylabel('%','FontSize',fontSize);
            xlim([-axLim axLim])
            ylim([0 0.5])
            set(gca,'XTick',[])
            set(gca,'YTick',[0 0.5],'YTickLabel',{'0', '50'})
            set(gca,'FontSize',fontSize)
        end
    end
    
    for j = 1:numel(ax(i).Children)
        if strcmp(ax(i).Children(j).Type,'image')
            ax(i).Children(j).CData = imrotate(ax(i).Children(j).CData', ...
                -rot(i-4)*180/pi,'bilinear','crop');
        elseif strcmp(ax(i).Children(j).Type,'contour')
            ax(i).Children(j).ZData = ax(i).Children(j).ZData';
            ax(i).Children(j).LineColor = 'k';
            ax(i).Children(j).LineStyle = 'none';
        elseif strcmp(ax(i).Children(j).Type,'line') && ax(i).Children(j).LineWidth < 3
            XData = ax(i).Children(j).YData;
            YData = ax(i).Children(j).XData;
            ax(i).Children(j).XData = XData;
            ax(i).Children(j).YData = YData;
        end
    end
end

for i = 1:2
    yl(i).Position = [-1.08e3 0.25 -1];
end

labels = ['(a)';'(c)';'(e)';'(b)';'(d)';'(f)'];

annotation('textbox',[0.01 0.905 0.1 0.1],'String',labels(1,:), ...
    'FontSize',fontSize,'EdgeColor','w')
annotation('textbox',[0.01 0.56 0.1 0.1],'String',labels(2,:), ...
    'FontSize',fontSize,'EdgeColor','w')
annotation('textbox',[0.01 0.217 0.1 0.1],'String',labels(3,:), ...
    'FontSize',fontSize,'EdgeColor','w')
annotation('textbox',[0.51 0.905 0.1 0.1],'String',labels(4,:), ...
    'FontSize',fontSize,'EdgeColor','w')
annotation('textbox',[0.51 0.56 0.1 0.1],'String',labels(5,:), ...
    'FontSize',fontSize,'EdgeColor','w')
annotation('textbox',[0.51 0.217 0.1 0.1],'String',labels(6,:), ...
    'FontSize',fontSize,'EdgeColor','w')

close(fig1)
close(fig2)
close(fig3)
close(fig4)
close(fig5)
close(fig6)
close(fig0)

set(gcf, 'PaperSize', get(gcf, 'PaperSize')*2)
print('Figures\Figure - Non-ideal Conditions\FinalSupp.jpg','-djpeg','-r100')
print('Figures\Figure - Non-ideal Conditions\FinalSupp', '-dpdf')

%% Figure - Depth

% Params
figSize = [1 1 9 10]*2;
fontSize = 8*2;
marTop = 0.04;
marBot = 0.22;
marLef = 0.12;
marRig = 0.05;
spacHor = 0.06;
spacVer = 0.11;

set(0,'defaultAxesFontName', 'Arial')

% Load subplots
fig1 = openfig('Figures\Figure - Depth\RotNeurite_OneElectrodez.fig');
fig2 = openfig('Figures\Figure - Depth\RotNeurite_OneElectrodex.fig');

% Combine
figDepth = figure('Units','centimeters','Position',figSize,'Name',...
    'Figure - Depth','Color','w');
set(figDepth, 'PaperPosition',figSize)

ax(1) = subaxis(2,1,1,'SpacingHoriz',spacHor,'SpacingVert',spacVer,'MarginRight',marRig,'MarginTop',marTop, ...
    'MarginLeft',marLef,'MarginBottom',marBot);
copyobj(fig1.CurrentAxes.Children,ax(1))
xlabel(ax(1),'Y (\mum)')
% title('Electrode','FontWeight','normal','FontSize',fontSize)

ax(2) = subaxis(2,1,2,'SpacingHoriz',spacHor,'SpacingVert',spacVer,'MarginRight',marRig,'MarginTop',marTop, ...
    'MarginLeft',marLef,'MarginBottom',marBot);
copyobj(fig2.CurrentAxes.Children,ax(2))
xlabel(ax(2),'X (\mum)')

labels = ['(a)'; '(b)'];

for i = 1:numel(ax)
    set(ax(i),'fontsize',fontSize)
    ylabel(ax(i),'Z (\mum)')
    set(figDepth, 'CurrentAxes', ax(i))
    colormap(flipud(hot))
    caxis([0 1.2])
    xlim([-1000 1000])
    ylim([-250 50])
    ax(i).Children(1).LineWidth = 1.5;
    line([-1000 1000], [0 0], 'LineStyle','--','LineWidth',0.5,'Color','k')
    line([-1000 1000], [-150 -150], 'LineStyle','--','LineWidth',0.5,'Color','k')
    text(-980,35,'Vitreous','FontSize',fontSize)
    text(-980,-15,'NFL','FontSize',fontSize)
    text(-980,-165,'GCL','FontSize',fontSize)
    text(-1280,120,labels(i,:),'FontSize',fontSize)
    
    ax_tmp = axes('Position',get(ax(i),'position'), ...
        'Xlim',get(ax(i),'Xlim'),'Ylim',get(ax(i),'Ylim'),'color','none', ...
        'Box','on');
    ax_tmp.XTick = ax(i).XTick;
    ax_tmp.XTickLabel = {};
    ax(i).XTickLabel = ax(i).XTickLabel;
    ax(i).XTick = ax_tmp.XTick;
    ax_tmp.YTick = ax(i).YTick;
    ax_tmp.YTickLabel = {};
    ax(i).YTickLabel = ax(i).YTickLabel;
    ax(i).YTick = ax_tmp.YTick;
end

origSize = get(ax(2),'Position');
cb = colorbar(ax(2), 'SouthOutside');
ylabel(cb,'Normalized extracellular potential')
cb.Ticks = 0:0.25:1;
cb.YLim = [0 1];
set(ax(2),'Position',origSize)
set(cb,'fontsize',fontSize)

close(fig1)
close(fig2)

print('Figures\Figure - Depth\Final.jpg','-djpeg','-r100')
export_fig 'Figures\Figure - Depth\Final' -eps -m8 -q101