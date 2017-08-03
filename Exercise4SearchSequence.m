function Exercise4SearchSequence(varargin)

% Exercise4SearchSequence
%
% Searches for optimal sequence in Exercise 4.
%
% Convention for refering to variables in equations:
%   x_t   refers to x{t}
%   x_tF  refers to x{t+1}
%   x_tL refers to x{t-1}
%   x_ss  refers to the steady state level of 'x'
%   (where x refers to some variable with name 'x')
%
% Required m-files:
%   - symbolic toolbox
%   - LQ package
%   - csolve.m, available in Chris Sims's website
%
% See also:
% LQ, LQSolveREE, LQCheckSOC, LQCheckSOCold, LQGenSymVar  
%
% .........................................................................
%
% Created: March 10, 2010 by Vasco Curdia
% Updated: April 22, 2014 by Vasco Curdia
%
% Copyright 2010-2014 by Vasco Curdia

%% ------------------------------------------------------------------------

%% preamble
% clear all
% tic
format short g

FileNameSuffix = '_dSP_4_Pers_90_Bind_10bp';
Shocks2Plot = {'hchitil'};

nsteps = 101;

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

%% Initialize:
psi = [];

%% Load data
load(['Output_Exercise4',FileNameSuffix])

TZLB.hchitil.NoCP = 0;
TZLB.hchitil.OptCP = 0;
if strcmp(FileNameSuffix,'_dSP_4_Pers_90_Bind_10bp')
    nCP1.hchitil.OptCP = 2;
else
    nCP1.hchitil.OptCP = 0;
end
nNoCP1.hchitil.OptCP = 0;
nCP2.hchitil.OptCP = 0;

TZLB.hXitil.NoCP = 0;%3
TZLB.hXitil.OptCP = 0;
if strcmp(FileNameSuffix,'_dSP_4_Pers_90')
    nCP1.hXitil.OptCP = 50;
elseif strcmp(FileNameSuffix,'_dSP_4_Pers_90_Bind_10bp')
    nCP1.hXitil.OptCP = 10;
else
    nCP1.hXitil.OptCP = 0;
end
nNoCP1.hXitil.OptCP = 0;
nCP2.hXitil.OptCP = 0;

TZLB.hchitiladd.NoCP = 0;
TZLB.hchitiladd.OptCP = 0;
nCP1.hchitiladd.OptCP = 0;
nNoCP1.hchitiladd.OptCP = 0;
nCP2.hchitiladd.OptCP = 0;

TZLB.hXitiladd.NoCP = 0;
TZLB.hXitiladd.OptCP = 0;
nCP1.hXitiladd.OptCP = 0;
nNoCP1.hXitiladd.OptCP = 0;
nCP2.hXitiladd.OptCP = 0;

%% Generate IRFs
fprintf('\nGenerating IRFs...\n')
for jS=1:ncsi,Sj = csi{jS};
    fprintf('%s...\n',Sj)
    % OptCP
    TZLBj=TZLB.(Sj).OptCP;
    nCP1j=nCP1.(Sj).OptCP;
    nNoCP1j=nNoCP1.(Sj).OptCP;
    nCP2j=nCP2.(Sj).OptCP;
    TCP2j = nCP2j;
    TNoCP1j = TCP2j+nNoCP1j;
    TCP1j = TNoCP1j+nCP1j;
    % Generate REE
    REEj(nsteps) = REE.OptCP.NoZLB.NoCP;
    for t=nsteps-1:-1:1
        if t<=TZLBj,ZLBj = 'ZLB';else ZLBj = 'NoZLB';end
        if ismember(t,[1:TCP2j,TNoCP1j+1:TCP1j])
            CPj = 'CP';
        else
            CPj = 'NoCP';
        end
        G0j = Mat.OptCP.(ZLBj).(CPj).G0;
        G1j = Mat.OptCP.(ZLBj).(CPj).G1;
        Cj  = Mat.OptCP.(ZLBj).(CPj).C;
        G2j = Mat.OptCP.(ZLBj).(CPj).G2;
        G3j = Mat.OptCP.(ZLBj).(CPj).G3;
        cv = find(~all(G3j==0,2));
        Cj(cv) = Cj(cv)-G0j(cv,:)*REEj(t+1).C;
        G0j(cv,:) = G0j(cv,:)*REEj(t+1).Phi1-G1j(cv,:);
        G1j(cv,:) = 0;
        G0ji = rbinv(G0j);
        REEj(t).Phi1 = G0ji*G1j;
        REEj(t).Phi2 = G0ji*G2j;
        REEj(t).C = G0ji*Cj;
    end
    % Generate IRF
    irfj = REEj(1).C+REEj(1).Phi2*ShockSize(:,jS);
    for t=2:nsteps
        irfj(:,t) = REEj(t).C+REEj(t).Phi1*irfj(:,t-1);
    end
    clear REEj
    % check solution
    CheckZLBj = all((irfj(idxRd,:)>(log(1/Rd_ss)-NumPrecision)));
    CheckUpsilonj = all((irfj(idxUpsilon,:)>-NumPrecision));
    Checklcbj = all(irfj(idxlcb,:)>-NumPrecision);
    Checkzetaj = all(irfj(idxzeta,:)>-zeta_ss-NumPrecision);
    isFoundSequence = (CheckZLBj && CheckUpsilonj && Checklcbj && Checkzetaj);
    if ~isFoundSequence
        fprintf('WARNING: Did not find sequence for OptCP!\n')
    end
    IRF.(Sj).OptCP = irfj(1:nz+nF,:);
    CheckZLB.(Sj).OptCP = CheckZLBj;
    CheckUpsilon.(Sj).OptCP = CheckUpsilonj;
    Checklcb.(Sj).OptCP = Checklcbj;
    Checkzeta.(Sj).OptCP = Checkzetaj;
    
    % NoCP
    TZLBj=TZLB.(Sj).NoCP;
    % Generate REE
    REEj(nsteps) = REE.NoCP.NoZLB;
    for t=nsteps-1:-1:1
        if t<=TZLBj,ZLBj = 'ZLB';else ZLBj = 'NoZLB';end
        G0j = Mat.NoCP.(ZLBj).G0;
        G1j = Mat.NoCP.(ZLBj).G1;
        Cj  = Mat.NoCP.(ZLBj).C;
        G2j = Mat.NoCP.(ZLBj).G2;
        G3j = Mat.NoCP.(ZLBj).G3;
        cv = find(~all(G3j==0,2));
        Cj(cv) = Cj(cv)-G0j(cv,:)*REEj(t+1).C;
        G0j(cv,:) = G0j(cv,:)*REEj(t+1).Phi1-G1j(cv,:);
        G1j(cv,:) = 0;
        G0ji = rbinv(G0j);
        REEj(t).Phi1 = G0ji*G1j;
        REEj(t).Phi2 = G0ji*G2j;
        REEj(t).C = G0ji*Cj;
    end
    % Generate IRF
    irfj = REEj(1).C+REEj(1).Phi2*ShockSize(:,jS);
    for t=2:nsteps
        irfj(:,t) = REEj(t).C+REEj(t).Phi1*irfj(:,t-1);
    end
    clear REEj
    % check solution
    CheckZLBj = all((irfj(idxRd,:)>(log(1/Rd_ss)-NumPrecision)));
    CheckUpsilonj = all((irfj(idxUpsilon,:)>-NumPrecision));
    isFoundSequence = CheckZLBj && CheckUpsilonj;
    if ~isFoundSequence
        fprintf('WARNING: Did not find sequence for OptCP!\n')
    end
    IRF.(Sj).NoCP = irfj(1:nz+nF,:);
    CheckZLB.(Sj).NoCP = CheckZLBj;
    CheckUpsilon.(Sj).NoCP = CheckUpsilonj;
end
clear irfj

%% Add more variables
for jS=1:ncsi,Sj = csi{jS};
    for jPol=1:nPol,Polj = PolList{jPol};
        irf = IRF.(Sj).(Polj);
        % add borrowing rate
        irf(idxRb,:) = irf(idxRd,:)+irf(idxomegatil,:);
        % add Long borrowing rate (RBaa)
        irf(idxRBaa,:) = NaN(1,nsteps);
        for t=1:nsteps-Duration_Baa+1
            irf(idxRBaa,t) = sum(irf(idxRb,t-1+(1:Duration_Baa)))/Duration_Baa;
        end
        % add Long savings rate (RT10y)
        irf(idxRT10y,:) = NaN(1,nsteps);
        for t=1:nsteps-Duration_T10y+1
            irf(idxRT10y,t) = sum(irf(idxRd,t-1+(1:Duration_T10y)))/Duration_T10y;
        end
        % add Long Spread
        irf(idxLongSP,:) = irf(idxRBaa,:)-irf(idxRT10y,:);
        % add deposit rate level (annualized pct)
        irf(idxRdLevel,:) = (Rd_ss*exp(irf(idxRd,:))).^4-1;
        % add real deposit rate
        irf(idxRrd,:) = irf(idxRd,:)-cat(2,irf(idxPi,2:end),NaN(1,1));
        % add cb
        irf(idxcb,:) = -sigma_b*irf(idxlambda_b,:);
        % add cs
        irf(idxcs,:) = -sigma_s*irf(idxlambda_s,:);
        % add w
        irf(idxw,:) = -(pi_b*(psi/psi_b*lambda_b_ss/lambdatil_ss)^(1/nu)*irf(idxlambda_b,:)+...
            (1-pi_b)*(psi/psi_s*lambda_s_ss/lambdatil_ss)^(1/nu)*irf(idxlambda_s,:))+...
            (1+omega_y)*irf(idxY,:)+irf(idxDelta,:);
        % add wb
        irf(idxwb,:) = irf(idxw,:)+...
            (1-pi_b)/nu*(psi/psi_s*lambda_s_ss/lambdatil_ss)^(1/nu)*...
            (irf(idxlambda_s,:)-irf(idxlambda_b,:));
        % add ws
        irf(idxws,:) = irf(idxw,:)-...
            pi_b/nu*(psi/psi_b*lambda_b_ss/lambdatil_ss)^(1/nu)*...
            (irf(idxlambda_s,:)-irf(idxlambda_b,:));
        % add Omega
        irf(idxOmega,:) = irf(idxlambda_b,:)-irf(idxlambda_s,:);
        % add CB credit as fraction og steady state debt
        irf(idxgammacb,:) = 1/b_ss*irf(idxlcb,:);
        % add fraction of CB credit (level)
        irf(idxgammacbLevel,:) = irf(idxgammacb,:);
        % add Total credit
        irf(idxlTot,:) = b_ss.*exp(irf(idxb,:));
        % add CB credit (level)
        irf(idxlCB,:) = lcb_ss+irf(idxlcb,:);
        % add Private credit
        irf(idxlPriv,:) = irf(idxlTot,:)-irf(idxlCB,:);
        % add level of zeta
        irf(idxzetalevel,:) = irf(idxzeta,:)+zeta_ss;
        % add level of varphi_omega_level
        irf(idxvarphi_omega_level,:) = -irf(idxFLM5,:)-FLM5_ss;
        % add level of varphi_Xi_level
        irf(idxvarphi_Xi_level,:) = -irf(idxFLM2,:)-FLM2_ss;
        % save the irf
        IRF.(Sj).(Polj) = 100*irf;
        clear irf
    end
end

%% display checks and number of periods
for jS=1:ncsi,Sj = csi{jS};
    if ~ismember(Sj,Shocks2Plot),continue,end
    fprintf('\nShock: %s\n',Sj)
    fprintf('\n  CheckZLB:\n')
    disp(CheckZLB.(Sj))
    fprintf('\n  CheckUpsilon:\n')
    disp(CheckUpsilon.(Sj))
    fprintf('  Checklcb:\n')
    disp(Checklcb.(Sj))
    fprintf('  Checkzeta:\n')
    disp(Checkzeta.(Sj))
    if nMax.CP2>0
        fprintf('  nCP2:\n')
        disp(nCP2.(Sj))
    end
    if nMax.NoCP1>0
        fprintf('  nNoCP1:\n')
        disp(nNoCP1.(Sj))
    end
    fprintf('  nCP1:\n')
    disp(nCP1.(Sj))
    fprintf('  TZLB:\n')
    disp(TZLB.(Sj))
end

%% ------------------------------------------------------------------------

%% save information
fprintf('\nMAT file: %s\n',ExerciseName)
save(ExerciseName)

%% Plot
IRFPlotCompareExercise4(...
    'FileNameSuffix',FileNameSuffix,...
    'Shocks2Plot',Shocks2Plot,...
    'FigShape',{2,2},...
    'NewFig',0)

%% Elapsed time
% disp(' '), vctoc(), disp(' ')
