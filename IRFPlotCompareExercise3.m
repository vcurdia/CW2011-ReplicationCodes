function IRFPlotCompareExercise3(varargin)

% IRFPlotCompareExercise3
%
% Generate plots of IRFs comparing different specifications, or policies (not both)
%
% .........................................................................
%
% Created: February 24, 2010 by Vasco Curdia
% Updated: April 21, 2014 by Vasco Curdia
%
% Copyright 2010-2014 by Vasco Curdia

%% ------------------------------------------------------------------------

%% preamble
nsteps = 21;
yMaxSlack = []; % 0.001
yMinScale = 1e-2; % 0.01
FigShow = 1;
FigPrint = 0;
KeepEPS = 1;
OpenPDF = 0;
SaveFigData = 0;

ShowLegend = 1;

XTickStep = 4;

Shocks2Plot = {'hchitil','hXitil','hchitiladd','hXitiladd'};
Shocks2PlotPretty = {'Mult \chi','Mult \Xi','Additive \chi','Additive \Xi'};
% Shocks2Plot = {'hchitil','hXitil','hchitiladd'};
% Shocks2PlotPretty = {'Mult \chi','Mult \Xi','Additive \chi'};

LineStyle = {'-','--','--+','--x','--o','--s'};
MarkerSize = {1,1,7,7,3,3};
LineColor = {'b','r','k',[0,0.5,0],[0,0.5,0.5],[0.87,0.49,0]};
LineWidth = 1;

FileNameSuffix = '_dSP_12_Pers_90'; 
FigShape = {2,2};
FigPrefix = '';

%% Update options
if ~isempty(varargin)
  nOptions = length(varargin);
  if nOptions==1 && isstruct(varargin{1})
    Options = fieldnames(varargin{1});
    for jO=1:length(Options)
      eval(sprintf('%1$s = varargin{1}.%1$s;',Options{jO}))
    end
  elseif mod(nOptions,2)
    error('Incorrect number of optional arguments.')
  else
    for jO=1:nOptions/2
      eval(sprintf('%s = varargin{%.0f};',varargin{(jO-1)*2+1},jO*2))
    end
  end
end

%% designate and label the variables to plot and scale
if all([FigShape{:}]==[1,1])
    var_plot = {'L_LHS'};
    var_label = {''};
    scale = [4]; % annualize inflation and interest rates
elseif all([FigShape{:}]==[1,2])
    var_plot = {'L_LHS','L_LHS'};
    var_label = {'Ignoring the ZLB','Accounting for the ZLB'};
    scale = [4,4]; % annualize inflation and interest rates
elseif all([FigShape{:}]==[2,1])
    var_plot = {'L_LHS','L_LHS'};
    var_label = {'Optimal Interest Rate Policy','Taylor Rule'};
    scale = [4,4]; % annualize inflation and interest rates
elseif all([FigShape{:}]==[2,2])
    var_plot = {'L_LHS','L_LHS','L_LHS','L_LHS'};
    var_label = {'','','',''};
    ColLabel = {'Ignoring the ZLB','Accounting for the ZLB'};
    RowLabel = {'Optimal Interest Rate Policy','Taylor Rule'};
    scale = [4,4,4,4]; % annualize inflation and interest rates
elseif all([FigShape{:}]==[4,3])
    var_plot = {'Y','Pi','RdLevel','omegatil','b','gammacb','cs','cb','Omega',...
        'zetalevel','L_LHS','Upsilon'};
    var_label = {'Y','\pi','i^d (level)','\omega','b','l^{cb}','c^s','c^b','\Omega',...
        '\zeta (level)','LLHS','\Upsilon'};
    scale = [1,4,1,4,1,1,1,1,1,1,4,1]; % annualize inflation and interest rates
end
nPlots = length(var_plot);

%% load mat file
load(['Output_Exercise3',FileNameSuffix],'zz','nzz','IRF','csi',...
    'Rd_ss','zeta_ss','LMNN_ss');
nSP = length(Shocks2Plot);

%% for panel only
if all([FigShape{:}]==[1,2])
    IgnoreZLB = load(['Output_Exercise3',FileNameSuffix,'_IgnoreZLB'],'IRF');
elseif all([FigShape{:}]==[2,1])
    OptPol = load(['Output_Exercise2',FileNameSuffix],'IRF');
elseif all([FigShape{:}]==[2,2])
    IgnoreZLB = load(['Output_Exercise3',FileNameSuffix,'_IgnoreZLB'],'IRF');
    OptPol = load(['Output_Exercise2',FileNameSuffix],'IRF');
    OptPolIgnoreZLB = load(['Output_Exercise2',FileNameSuffix,'_IgnoreZLB'],'IRF');
end

%% Plot IRFs
tid = 0:1:nsteps-1; ntid = length(tid);
nS = length(Shocks2Plot);
nIRF = nS;
if FigShow
    figure
else
    figure('Visible','off')
end
FigData = cell(nPlots,1);
for jj=1:nPlots
    hsubp(jj) = subplot(FigShape{:},jj);
    IRF2Plot = NaN(nIRF,ntid);
    for jS=1:nS
        [tfS,shock_pos] = ismember(Shocks2Plot{jS},csi);
        [tfV,var_pos] = ismember(var_plot{jj},zz);
        if tfS && tfV
            if all([FigShape{:}]==[1,2]) && jj==1
                IRF2Plot(jS,:) = scale(jj)*IgnoreZLB.IRF.(csi{shock_pos})(var_pos,1:nsteps);
            elseif all([FigShape{:}]==[2,1]) && jj==1
                IRF2Plot(jS,:) = scale(jj)*OptPol.IRF.(csi{shock_pos})(var_pos,1:nsteps);
            elseif all([FigShape{:}]==[2,2])
                if jj==1
                    IRF2Plot(jS,:) = scale(jj)*OptPolIgnoreZLB.IRF.(csi{shock_pos})(var_pos,1:nsteps);
                elseif jj==2
                    IRF2Plot(jS,:) = scale(jj)*OptPol.IRF.(csi{shock_pos})(var_pos,1:nsteps);
                elseif jj==3
                    IRF2Plot(jS,:) = scale(jj)*IgnoreZLB.IRF.(csi{shock_pos})(var_pos,1:nsteps);
                else
                    IRF2Plot(jS,:) = scale(jj)*IRF.(csi{shock_pos})(var_pos,1:nsteps);
                end
            else
                IRF2Plot(jS,:) = scale(jj)*IRF.(csi{shock_pos})(var_pos,1:nsteps);
            end
        end
    end
    for jS=1:nS
        plot(tid,IRF2Plot(jS,:),LineStyle{jS},...
            'Color',LineColor{jS},'LineWidth',LineWidth,...
            'MarkerSize',MarkerSize{jS},'MarkerFaceColor',LineColor{jS})
        hold on
    end
    if ismember('\varphi',var_label{jj})
        h=title('');
        set(h,'Interpreter','latex');
        set(h,'String',['$',var_label{jj},'$']);
    else
        title(var_label{jj})
    end
    if all([FigShape{:}]==[2,2])
        if ismember(jj,[1,2])
            title(ColLabel{jj})
        end
        if ismember(jj,[1,3])
            ylabel(RowLabel{(jj-1)/2+1})
        end
    end
    if strcmp(var_plot{jj},'RdLevel')
        RefLevel = (Rd_ss^4-1)*100;
    elseif strcmp(var_plot{jj},'L_LHS')
        RefLevel = LMNN_ss*400;
    elseif strcmp(var_plot{jj},'zetalevel')
        RefLevel = zeta_ss*100;
    else
        RefLevel = 0;
    end
    plot(tid,RefLevel*ones(size(tid)),'k:')
    FigData{jj} = [tid;IRF2Plot;RefLevel*ones(size(tid))]';
    xlim([0 nsteps-1])
    set(gca,'XTick',0:XTickStep:nsteps)
    yMax = max([RefLevel,max(IRF2Plot(:))]);
    yMin = min([RefLevel,min(IRF2Plot(:))]);
    ySlack = max([0.05*(yMax-yMin),yMaxSlack]);
    ylim([min([yMin-ySlack,RefLevel-yMinScale]) max([yMax+ySlack,RefLevel+yMinScale])])
    if all([FigShape{:}]==[1,2]) || all([FigShape{:}]==[2,1]) || all([FigShape{:}]==[2,2])
        yBounds(jj,:) = [min([yMin-ySlack,RefLevel-yMinScale]) max([yMax+ySlack,RefLevel+yMinScale])];
    end
end
if all([FigShape{:}]==[1,2]) || all([FigShape{:}]==[2,1]) || all([FigShape{:}]==[2,2])
    yBoundsCommon = [min(yBounds(:,1)),max(yBounds(:,2))];
    for jj=1:nPlots
        set(hsubp(jj),'YLim',yBoundsCommon)
    end
end
if ShowLegend
    if all([FigShape{:}]==[4,3])
        hleg = legend(Shocks2PlotPretty{:},'Orientation','horizontal');
        legPos = get(hleg,'Position');
        xL = get(hsubp((FigShape{1}-1)*FigShape{2}+1),'Position');
        xR = get(hsubp(FigShape{1}*FigShape{2}),'Position');
        legPos(1) = xL(1)+(xR(1)-xL(1))/2+(xL(3)-legPos(3))/2;
        legPos(2) = 0;
        set(hleg,'Position',legPos)
    else
        legend(Shocks2PlotPretty{:},'Location','NE')
    end
end
%     legend('boxoff')
% convert plot to an eps file
if FigPrint
    if all([FigShape{:}]==[1,1])
        FigName = [FigPrefix,'Plot_LagMultiplier',FileNameSuffix,'_Taylor'];
    elseif all([FigShape{:}]==[1,2])
        FigName = [FigPrefix,'Plot_LagMultiplier',FileNameSuffix,'_Taylor_IgnoreZLB'];
    elseif all([FigShape{:}]==[2,1])
        FigName = [FigPrefix,'Plot_LagMultiplier',FileNameSuffix,'_Taylor_OptIntPol'];
    elseif all([FigShape{:}]==[2,2])
        FigName = [FigPrefix,'Plot_LagMultiplier',FileNameSuffix,'_Taylor_OptIntPol_IgnoreZLB'];
    else
        FigName = [FigPrefix,'Plot_LagMultiplier',FileNameSuffix,'_Taylor_AllVars'];
    end
%     print('-depsc2',[FigName,'.eps'])
    vcPrintPDF(FigName,KeepEPS,OpenPDF)
end
if SaveFigData
    if all([FigShape{:}]==[1,1])
        FigName = [FigPrefix,'Data_Plot_LagMultiplier',FileNameSuffix,'_Taylor'];
    elseif all([FigShape{:}]==[1,2])
        FigName = [FigPrefix,'Data_Plot_LagMultiplier',FileNameSuffix,'_Taylor_IgnoreZLB'];
    elseif all([FigShape{:}]==[2,1])
        FigName = [FigPrefix,'Data_Plot_LagMultiplier',FileNameSuffix,'_Taylor_OptIntPol'];
    elseif all([FigShape{:}]==[2,2])
        FigName = [FigPrefix,'Data_Plot_LagMultiplier',FileNameSuffix,'_Taylor_OptIntPol_IgnoreZLB'];
    else
        FigName = [FigPrefix,'Data_Plot_LagMultiplier',FileNameSuffix,'_Taylor_AllVars'];
    end
    save(FigName,'FigData')
end

%% ------------------------------------------------------------------------
