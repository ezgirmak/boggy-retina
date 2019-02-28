addpath([fileparts(pwd()) '\Utilities'])
addpath([fileparts(pwd()) '\Utilities\altmany-export_fig-04ca93c'])
addpath([fileparts(pwd()) '\Utilities\altmany-export_fig-04ca93c\xpdf-tools-win-4.00'])
addpath(fileparts(pwd()))

set(0,'defaultfigurecolor',[1 1 1])

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

axis_limits = [-200 100];
NFL_GCL_boundary = -100;

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
    ylim(axis_limits)
    ax(i).Children(1).LineWidth = 1.5;
    line([-1000 1000], [0 0], 'LineStyle','--','LineWidth',0.5,'Color','k')
    line([-1000 1000], [NFL_GCL_boundary NFL_GCL_boundary], 'LineStyle','--','LineWidth',0.5,'Color','k')
    text(-980,35,'Vitreous','FontSize',fontSize)
    text(-980,-15,'NFL','FontSize',fontSize)
    text(-980,NFL_GCL_boundary-15,'GCL','FontSize',fontSize)
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