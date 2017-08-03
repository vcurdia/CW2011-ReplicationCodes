function IRFPlotCompareExercise4(varargin)

% IRFPlotCompareExercise4
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
nsteps = 41;
yMaxSlack = []; % 0.001
yMinScale = 1e-2; % 0.01
FigShow = 1;
NewFig = 1;
FigPrint = 0;
KeepEPS = 1;
OpenPDF = 0;
SaveFigData = 0;

ShowLegend = 1;
XTickStep = 8; %(nsteps-1)/10; %16

% Shocks2Plot = {'hchitil','hXitil','hchitiladd','hXitiladd'};
Shocks2Plot = {'hchitil'};

Pol2Plot = {'OptCP','NoCP'};
Pol2PlotPretty = {'Optimal Credit Policy','No Credit Policy'};
nPol = length(Pol2Plot);


LineStyle = {'-','--','--+','--x','--o','--s'};
MarkerSize = {1,1,5,5,3,3};
LineColor = {'b','r','k',[0,0.5,0],[0,0.5,0.5],[0.87,0.49,0]};
LineWidth = 1;

FileNameSuffix = '_dSP_4_Pers_90'; 
FigShape = {3,3};
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
if all([FigShape{:}]==[3,3])
    var_plot = {'Y','Pi','RdLevel','omegatil','b','gammacb','cs','cb','Omega'};
    var_label = {'Y','\pi','i^d (level)','\omega','b','L^{cb}','c^s','c^b','\Omega'};
    scale = [1,4,1,4,1,1,1,1,1]; % annualize inflation and interest rates
elseif all([FigShape{:}]==[2,3])
    var_plot = {'Y','Pi','Rd','omegatil','b','gammacb'};
    var_label = {'Y','\pi','i^d','\omega','b','L^{cb}'};
    scale = [1,4,4,4,1,1]; % annualize inflation and interest rates
elseif all([FigShape{:}]==[4,3])
    var_plot = {'Y','Pi','RdLevel','omegatil','b','gammacb','cs','cb','Omega',...
        'zetalevel','varphi_Xi_level','varphi_omega_level'};
    var_label = {'Y','\pi','i^d (level)','\omega','b','L^{cb}','c^s','c^b','\Omega',...
        '\zeta (level)','\varphi_\Xi (level)','\varphi_\omega (level)'};
    scale = [1,4,1,4,1,1,1,1,1,1,1,1]; % annualize inflation and interest rates
elseif all([FigShape{:}]==[2,2])
    var_plot = {'RdLevel','Upsilon','gammacb','zetalevel'};
    var_label = {'i^d (level)','\Upsilon','L^{cb}','\zeta (level)'};
    scale = [1,1,1,1]; % annualize inflation and interest rates
end
nPlots = length(var_plot);

%% load mat file
load(['Output_Exercise4',FileNameSuffix],'zz','nzz','IRF','csi',...
    'Rd_ss','zeta_ss','FLM2_ss','FLM5_ss');

%% Plot IRFs
tid = 0:1:nsteps-1; ntid = length(tid);
nS = length(Shocks2Plot);
nIRF = nPol;
for j=1:nS, Sj = Shocks2Plot{j};
    if ~ismember(Sj,csi), continue, end
    if NewFig
        if FigShow
            figure('Name',sprintf('Responses to a shock in %s',Sj))
        else
            figure('Visible','off')
        end
    end
    FigData = cell(nPlots,1);
    for jj=1:nPlots
        hsubp(jj) = subplot(FigShape{:},jj);
        IRF2Plot = NaN(nIRF,ntid);
        for jP=1:nPol
            [tf,var_pos] = ismember(var_plot{jj},zz);
            if tf
                IRF2Plot(jP,:) = scale(jj)*IRF.(Sj).(Pol2Plot{jP})(var_pos,1:nsteps);
            end
        end
        for jP=1:nPol
            plot(tid,IRF2Plot(jP,:),LineStyle{jP},...
                'Color',LineColor{jP},'LineWidth',LineWidth,...
                'MarkerSize',MarkerSize{jP},'MarkerFaceColor',LineColor{jP})
            hold on
        end
        if ismember('\varphi',var_label{jj})
            h=title('');
            set(h,'Interpreter','latex');
            set(h,'String',strrep(['$',var_label{jj},'$'],' (level)$','$ (level)'));
        else
            title(var_label{jj})
        end
        xlim([0 nsteps-1])
        set(gca,'XTick',0:XTickStep:nsteps-1)
%         if strcmp(var_plot{jj},'RdLevel')
%             RefLevel = (Rd_ss^4-1)*100;
%         elseif strcmp(var_plot{jj},'zetalevel')
%             RefLevel = 100*zeta_ss;
        if strcmp(var_plot{jj},'varphi_Xi_level')
            RefLevel = -100*FLM2_ss;
        elseif strcmp(var_plot{jj},'varphi_omega_level')
            RefLevel = -100*FLM5_ss;
        else
            RefLevel = 0;
        end
        plot(tid,RefLevel*ones(size(tid)),'k:')
        hold off
        FigData{jj} = [tid;IRF2Plot;RefLevel*ones(size(tid))]';
        yMax = max([RefLevel,max(IRF2Plot(:))]);
        yMin = min([RefLevel,min(IRF2Plot(:))]);
        ySlack = max([0.05*(yMax-yMin),yMaxSlack]);
        ylim([min([yMin-ySlack,RefLevel-yMinScale]) max([yMax+ySlack,RefLevel+yMinScale])])
    end
    if ShowLegend
        hleg = legend(Pol2PlotPretty{:},'Orientation','horizontal');
        legPos = get(hleg,'Position');
        xL = get(hsubp((FigShape{1}-1)*FigShape{2}+1),'Position');
        xR = get(hsubp(FigShape{1}*FigShape{2}),'Position');
        legPos(1) = xL(1)+(xR(1)-xL(1))/2+(xL(3)-legPos(3))/2;
        legPos(2) = 0;
        set(hleg,'Position',legPos)
    end
    % convert plot to an eps file
    if FigPrint
        if all([FigShape{:}]==[4,3])
            FigName = strrep(sprintf('%sIRF_OptCP_OptIntPol%s_%s_WithMultipliers',...
                FigPrefix,FileNameSuffix,Sj),'_h','_');
        else
            FigName = strrep(sprintf('%sIRF_OptCP_OptIntPol%s_%s',...
                FigPrefix,FileNameSuffix,Sj),'_h','_');
        end
%         print('-depsc2',[FigName,'.eps'])
        vcPrintPDF(FigName,KeepEPS,OpenPDF)
    end
    if SaveFigData
        if all([FigShape{:}]==[4,3])
            FigName = strrep(sprintf('%sFigData_IRF_OptCP_OptIntPol%s_%s_WithMultipliers',...
                FigPrefix,FileNameSuffix,Sj),'_h','_');
        else
            FigName = strrep(sprintf('%sFigData_IRF_OptCP_OptIntPol%s_%s',...
                FigPrefix,FileNameSuffix,Sj),'_h','_');
        end
        save(FigName,'FigData')
    end
end

%% ------------------------------------------------------------------------
