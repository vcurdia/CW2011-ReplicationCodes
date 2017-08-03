function [RMSmat,RMSk_t]=Ramsey(RMSy_t,RMScsi_t,RMSU,RMSF,RMSG,RMSBETTA,RMSS,...
    RMSy_ss,FLM_ss,GLM_ss,yLogIdx)

% Ramsey
%
% This function solves for Ramsey optimal policy.
% 
%   [RMSmat,RMSk_t]=Ramsey(RMSy_t,RMScsi_t,RMSU,RMSF,RMSG,RMSBETTA,RMSS,...
%       RMSy_ss,RMSFLM_ss,RMSGLM_ss,yLogIdx)
%
% Model:
%   max sum_{t=0,...,inf} E_0[U(y{t},csi{t})]
%   s.t.  0 = F(y{t},csi{t},y{t-1})
%         0 = E_t[G(y{t},csi{t},y{t+1})]
%   with
%         csi{t} = S*csi{t-1} + eps{t}
%
% log-linear FOC is given by
%   0 = A0*E_t[hy{t+1}] + A1*hy{t} + RMSBETTA^(-1)*A0'*hy{t-1} 
%       + B0*csi{t} + B1*csi{t-1}
%       + C0'*FLM{t} + RMSBETTA*C1'*E_t[FLM{t+1}]
%       + D1'*GLM{t} + RMSBETTA^(-1)*D0'*GLM{t-1}
%   0 = C0*hy{t} + C1*hy{t-1} + C2*csi{t}
%   0 = D0*E_t[hy{t+1}] + D1*hy{t} + D2*csi{t}
%   0 = S*csi{t-1} + eps{t} - csi{t}
%
% Inputs:
%   - RMSy_t: symbolic array with endogenous state variables
%   - RMScsi_t: symbolic array with exogenous variables
%   - RMSU: symbolic period utility function to be maximized by authorities
%   - RMSF: symbolic system of possibly backward looking constraints
%   - RMSG: symbolic system of forward looking constraints
%   - RMSBETTA: symbolic expression for the time discount factor
%   - RMSS: symbolic autocorrelation matrix for exogenous variables
%   - RMSy_ss: symbolic array of steady state endogenous state variables
%   - FLM_ss: symbolic array of steady state lagrange multipliers of F
%     constraints
%   - GLM_ss: symbolic array of steady state lagrange multipliers of G
%     constraints
%   - yLogIdx: [optional] index of endogenous state variables to be 
%     converted to logs (binary vector of same size of y_t). If ommited all
%     endogenous variables are log-linearized.
%
% Outputs:
%   - RMSmat: structure with symbolic matrices used in the Ramsey solution
%     described above.
%       * RMSmat.A0, RMSmat.A1
%       * RMSmat.B0, RMSmat.B1
%       * RMSmat.C0, RMSmat.C1, RMSmat.C2
%       * RMSmat.D0, RMSmat.D1, RMSmat.D2
%       * RMSmat.G0, RMSmat.G1, RMSmat.G2, RMSmat.G3 
%         are the matrices such that
%           G0*k_tF = G1*k_t + G2*eps_tF + G3*eta_tF
%         with
%           k_t = (hy_t,csi_t,FLM_t,GLM_t,hyL_t,csiL_t,GLML_t)
%	- RMSmat.S: autocorrelation matrix for shocks
%	- RMSmat.beta: time discount factor
% 
% Required Matlab routines:
%   - Symbolic Toolbox
%
% See also:
% RMS, RMSSolveREE, LQ, LQSolveREE, LQCheckSOC, LQCheckSOCold, LQGenSymVar, 
% LQAltRule, LQWEval, MonFrictions, MonFrictionsIRFPlot, 
% MonFrictionsIRFPlotComp, MonFrictionsIRFPlotAltRule, SmetsWouters
%   
% .........................................................................
%
% Created: July 25, 2009
% Updated: September 18, 2009
% by Vasco Curdia

% -------------------------------------------------------------------------

% The previous information above can be accessed issuing the following
% command:
%    help Ramsey
% or
%    doc Ramsey

%% ------------------------------------------------------------------------

%% Setup some background information

% vpa inputs
inputList = {'RMSy_t','RMScsi_t','RMSU','RMSF','RMSG','RMSBETTA','RMSS','RMSy_ss','FLM_ss','GLM_ss'};
for j=1:length(inputList),eval([inputList{j},'=vpa(',inputList{j},');']),end

% sizes of variables/functions
ny=size(RMSy_t,2);
ncsi=size(RMScsi_t,2);
nF=size(RMSF,1);
nG=size(RMSG,1);
RMScsi_ss = sym(zeros(1,ncsi));

% generate tt and t_1 variables
RMSy_tF = RMSy_t;
RMSy_tL = RMSy_t;
for j=1:ny
    RMSy_tF = subs(RMSy_tF,RMSy_tF(j),[char(RMSy_tF(j)),'F'],0);
    RMSy_tL = subs(RMSy_tL,RMSy_tL(j),[char(RMSy_tL(j)),'L'],0);
end

%% ------------------------------------------------------------------------

%% Substitute variables with logs

% check if yLogIdx was not specified
if nargin<11
    yLogIdx = true(size(RMSy_t));
end

% create log-variables
hRMSy_t = RMSy_t;
hRMSy_tF = RMSy_tF;
hRMSy_tL = RMSy_tL;
hRMSy_ss = RMSy_ss;
for j=1:ny
    if yLogIdx(j)
        hRMSy_t = subs(hRMSy_t,RMSy_t(j),['h' char(RMSy_t(j))],0);
        hRMSy_tF = subs(hRMSy_tF,RMSy_tF(j),['h' char(RMSy_tF(j))],0);
        hRMSy_tL = subs(hRMSy_tL,RMSy_tL(j),['h' char(RMSy_tL(j))],0);
        hRMSy_ss(j) = log(RMSy_ss(j));
        % Plug in transformed variables into the model
        RMSF = subs(RMSF, [RMSy_t(j),RMSy_tL(j)], exp([hRMSy_t(j),hRMSy_tL(j)]),0);
        RMSG = subs(RMSG, [RMSy_t(j),RMSy_tF(j)], exp([hRMSy_t(j),hRMSy_tF(j)]),0);
        RMSU = subs(RMSU, [RMSy_t(j)], exp([hRMSy_t(j)]),0);
    end
end

%% ------------------------------------------------------------------------

%% Compute derivatives of U

RMSU_y = vpa(jacobian(RMSU,hRMSy_t));
RMSU_yy = vpa(jacobian(RMSU_y,hRMSy_t));
RMSU_ycsi = vpa(jacobian(RMSU_y,RMScsi_t));

idxSubs = find(RMSU_y~=0);
RMSU_y(idxSubs) = subs(RMSU_y(idxSubs),[hRMSy_t,RMScsi_t],[hRMSy_ss,RMScsi_ss],0);

idxSubs = find(RMSU_yy~=0);
RMSU_yy(idxSubs) = subs(RMSU_yy(idxSubs),[hRMSy_t,RMScsi_t],[hRMSy_ss,RMScsi_ss],0);

idxSubs = find(RMSU_ycsi~=0);
RMSU_ycsi(idxSubs) = subs(RMSU_ycsi(idxSubs),[hRMSy_t,RMScsi_t],[hRMSy_ss,RMScsi_ss],0);

%% Compute derivatives of F

RMSF_y = vpa(jacobian(RMSF,hRMSy_t));
RMSF_csi = vpa(jacobian(RMSF,RMScsi_t));
RMSF_yb = vpa(jacobian(RMSF,hRMSy_tL));

RMSF_yy = vpa(jacobian(RMSF_y(:),hRMSy_t));
RMSF_ycsi = vpa(jacobian(RMSF_y(:),RMScsi_t));
RMSF_ybcsi = vpa(jacobian(RMSF_yb(:),RMScsi_t));
RMSF_ybyb = vpa(jacobian(RMSF_yb(:),hRMSy_tL));
RMSF_yby = vpa(jacobian(RMSF_yb(:),hRMSy_t));

idxSubs = find(RMSF_y~=0);
RMSF_y(idxSubs) = subs(RMSF_y(idxSubs),[hRMSy_t,RMScsi_t,hRMSy_tL],[hRMSy_ss,RMScsi_ss,hRMSy_ss],0);

idxSubs = find(RMSF_csi~=0);
RMSF_csi(idxSubs) = subs(RMSF_csi(idxSubs),[hRMSy_t,RMScsi_t,hRMSy_tL],[hRMSy_ss,RMScsi_ss,hRMSy_ss],0);

idxSubs = find(RMSF_yb~=0);
RMSF_yb(idxSubs) = subs(RMSF_yb(idxSubs),[hRMSy_t,RMScsi_t,hRMSy_tL],[hRMSy_ss,RMScsi_ss,hRMSy_ss],0);

idxSubs = find(RMSF_yy~=0);
RMSF_yy(idxSubs) = subs(RMSF_yy(idxSubs),[hRMSy_t,RMScsi_t,hRMSy_tL],[hRMSy_ss,RMScsi_ss,hRMSy_ss],0);
RMSF_yy = reshape(RMSF_yy,nF,ny,ny);

idxSubs = find(RMSF_ycsi~=0);
RMSF_ycsi(idxSubs) = subs(RMSF_ycsi(idxSubs),[hRMSy_t,RMScsi_t,hRMSy_tL],[hRMSy_ss,RMScsi_ss,hRMSy_ss],0);
RMSF_ycsi = reshape(RMSF_ycsi,nF,ny,ncsi);

idxSubs = find(RMSF_ybcsi~=0);
RMSF_ybcsi(idxSubs) = subs(RMSF_ybcsi(idxSubs),[hRMSy_t,RMScsi_t,hRMSy_tL],[hRMSy_ss,RMScsi_ss,hRMSy_ss]);
RMSF_ybcsi = reshape(RMSF_ybcsi,nF,ny,ncsi);

idxSubs = find(RMSF_ybyb~=0);
RMSF_ybyb(idxSubs) = subs(RMSF_ybyb(idxSubs),[hRMSy_t,RMScsi_t,hRMSy_tL],[hRMSy_ss,RMScsi_ss,hRMSy_ss],0);
RMSF_ybyb = reshape(RMSF_ybyb,nF,ny,ny);

idxSubs = find(RMSF_yby~=0);
RMSF_yby(idxSubs) = subs(RMSF_yby(idxSubs),[hRMSy_t,RMScsi_t,hRMSy_tL],[hRMSy_ss,RMScsi_ss,hRMSy_ss],0);
RMSF_yby = reshape(RMSF_yby,nF,ny,ny);

%% Compute derivatives of G

if nG>0
    RMSG_y = vpa(jacobian(RMSG,hRMSy_t));
    RMSG_csi = vpa(jacobian(RMSG,RMScsi_t));
    RMSG_yf = vpa(jacobian(RMSG,hRMSy_tF));

    RMSG_yy = vpa(jacobian(RMSG_y(:),hRMSy_t));
    RMSG_ycsi = vpa(jacobian(RMSG_y(:),RMScsi_t));
    RMSG_yyf = vpa(jacobian(RMSG_y(:),hRMSy_tF));
    RMSG_yfcsi = vpa(jacobian(RMSG_yf(:),RMScsi_t));
    RMSG_yfyf = vpa(jacobian(RMSG_yf(:),hRMSy_tF));

    idxSubs = find(RMSG_y~=0);
    RMSG_y(idxSubs) = subs(RMSG_y(idxSubs),[hRMSy_t,RMScsi_t,hRMSy_tF],[hRMSy_ss,RMScsi_ss,hRMSy_ss],0);

    idxSubs = find(RMSG_csi~=0);
    RMSG_csi(idxSubs) = subs(RMSG_csi(idxSubs),[hRMSy_t,RMScsi_t,hRMSy_tF],[hRMSy_ss,RMScsi_ss,hRMSy_ss],0);

    idxSubs = find(RMSG_yf~=0);
    RMSG_yf(idxSubs) = subs(RMSG_yf(idxSubs),[hRMSy_t,RMScsi_t,hRMSy_tF],[hRMSy_ss,RMScsi_ss,hRMSy_ss],0);

    idxSubs = find(RMSG_yy~=0);
    RMSG_yy(idxSubs) = subs(RMSG_yy(idxSubs),[hRMSy_t,RMScsi_t,hRMSy_tF],[hRMSy_ss,RMScsi_ss,hRMSy_ss],0);
    RMSG_yy = reshape(RMSG_yy,nG,ny,ny);

    idxSubs = find(RMSG_ycsi~=0);
    RMSG_ycsi(idxSubs) = subs(RMSG_ycsi(idxSubs),[hRMSy_t,RMScsi_t,hRMSy_tF],[hRMSy_ss,RMScsi_ss,hRMSy_ss],0);
    RMSG_ycsi = reshape(RMSG_ycsi,nG,ny,ncsi);

    idxSubs = find(RMSG_yfcsi~=0);
    RMSG_yfcsi(idxSubs) = subs(RMSG_yfcsi(idxSubs),[hRMSy_t,RMScsi_t,hRMSy_tF],[hRMSy_ss,RMScsi_ss,hRMSy_ss],0);
    RMSG_yfcsi = reshape(RMSG_yfcsi,nG,ny,ncsi);

    idxSubs = find(RMSG_yfyf~=0);
    RMSG_yfyf(idxSubs) = subs(RMSG_yfyf(idxSubs),[hRMSy_t,RMScsi_t,hRMSy_tF],[hRMSy_ss,RMScsi_ss,hRMSy_ss],0);
    RMSG_yfyf = reshape(RMSG_yfyf,nG,ny,ny);

    idxSubs = find(RMSG_yyf~=0);
    RMSG_yyf(idxSubs) = subs(RMSG_yyf(idxSubs),[hRMSy_t,RMScsi_t,hRMSy_tF],[hRMSy_ss,RMScsi_ss,hRMSy_ss],0);
    RMSG_yyf = reshape(RMSG_yyf,nG,ny,ny);
else
    RMSG_y = [];
    RMSG_csi = [];
    RMSG_yf = [];

    RMSG_yy = [];
    RMSG_ycsi = [];
    RMSG_yyf = [];
    RMSG_yfcsi = [];
    RMSG_yfyf = [];
end

%% note about the derivatives
% all second order derivatives are stored in 3-dim arrays, corresponding to
%   dim_function x dim_vars_D1 x dim_vars_D2
% This means that for each call of squeeze(...(j,:,:)) yields
%   dim_vars_D1 x dim_vars_D2

%% ------------------------------------------------------------------------

%% Compute RMS approximation matrices

RMSA0 = zeros(ny);
RMSA1 = RMSU_yy; 
RMSB0 = RMSU_ycsi;
RMSB1 = zeros(ny,ncsi);
for j=1:nF
    RMSA0 = RMSA0+RMSBETTA*FLM_ss(j)*squeeze(RMSF_yby(j,:,:));
    RMSA1 = RMSA1+FLM_ss(j)*(squeeze(RMSF_yy(j,:,:))+RMSBETTA*squeeze(RMSF_ybyb(j,:,:)));
    RMSB0 = RMSB0+FLM_ss(j)*(squeeze(RMSF_ycsi(j,:,:))+RMSBETTA*squeeze(RMSF_ybcsi(j,:,:)));
end
for j=1:nG
    RMSA0 = RMSA0+GLM_ss(j)*squeeze(RMSG_yyf(j,:,:));
    RMSA1 = RMSA1+GLM_ss(j)*(squeeze(RMSG_yy(j,:,:))+RMSBETTA^(-1)*squeeze(RMSG_yfyf(j,:,:)));
    RMSB0 = RMSB0+GLM_ss(j)*squeeze(RMSG_ycsi(j,:,:));
    RMSB1 = RMSB1+GLM_ss(j)*RMSBETTA^(-1)*squeeze(RMSG_yfcsi(j,:,:));
end

RMSC0 = RMSF_y;
RMSC1 = RMSF_yb;
RMSC2 = RMSF_csi;

RMSD0 = RMSG_yf;
RMSD1 = RMSG_y;
RMSD2 = RMSG_csi;

%% ------------------------------------------------------------------------

%% Generate variables for the state space
% We now need to set up the state space system in the following form:
%   G0*k_tF.' = G1*k_t.' + G3*eps_tF.' + G4*eta_tF.'
% with
%   k_t = (hy_t,csi_t,FLM_t,GLM_t,hyL_t,csiL_t,GLML_t)

hRMSyL_t=hRMSy_t; hRMSyL_tF=hRMSyL_t;
for j=1:ny
    hRMSyL_t = subs(hRMSyL_t,hRMSyL_t(j),strrep(char(hRMSyL_t(j)),'_t','L_t'),0);
    hRMSyL_tF = subs(hRMSyL_tF,hRMSyL_tF(j),strrep(char(hRMSyL_tF(j)),'_t','L_tF'),0);
end

RMScsi_tF=RMScsi_t; RMScsiL_t=RMScsi_t; RMScsiL_tF=RMScsi_t; RMSeps_tF=RMScsi_t;
for j=1:ncsi
    RMScsi_tF = subs(RMScsi_tF,RMScsi_tF(j),[char(RMScsi_tF(j)),'t'],0);
    RMScsiL_t = subs(RMScsiL_t,RMScsiL_t(j),strrep(char(RMScsiL_t(j)),'_t','L_t'),0);
    RMScsiL_tF = subs(RMScsiL_tF,RMScsiL_tF(j),strrep(char(RMScsiL_tF(j)),'_t','L_tF'),0);
    RMSeps_tF = subs(RMSeps_tF,RMSeps_tF(j),['eps_',char(RMScsi_tF(j)),'t'],0);
end

FLM_t=[]; FLM_tF=[];
for j=1:nF
    eval(['syms FLM',int2str(j),'_t FLM',int2str(j),'_tF'])
    FLM_t = [FLM_t, eval(['FLM',int2str(j),'_t'])];
    FLM_tF = [FLM_tF, eval(['FLM',int2str(j),'_tF'])];
end

GLM_t=[]; GLM_tF=[];GLML_t=[]; GLML_tF=[];
for j=1:nG
    eval(['syms GLM',int2str(j),'_t GLM',int2str(j),'_tF'])
    GLM_t = [GLM_t, eval(['GLM',int2str(j),'_t'])];
    GLM_tF = [GLM_tF, eval(['GLM',int2str(j),'_tF'])];
    eval(['syms GLM',int2str(j),'L_t GLM',int2str(j),'L_tF'])
    GLML_t = [GLML_t, eval(['GLM',int2str(j),'L_t'])];
    GLML_tF = [GLML_tF, eval(['GLM',int2str(j),'L_tF'])];
end

RMSk_t = [hRMSy_t, RMScsi_t, FLM_t, GLM_t, hRMSyL_t, RMScsiL_t, GLML_t];
RMSk_tF = [hRMSy_tF, RMScsi_tF, FLM_tF, GLM_tF, hRMSyL_tF, RMScsiL_tF, GLML_tF];

%% Generate the state space system with the FOC

% FOC
SYS1 = RMSA0*hRMSy_tF.'+RMSA1*hRMSy_t.'+RMSBETTA^(-1)*RMSA0.'*hRMSyL_t.'+...
        RMSB0*RMScsi_t.'+RMSB1*RMScsiL_t.'+...
        RMSC0.'*FLM_t.'+RMSBETTA*RMSC1.'*FLM_tF.';
if nG>0
    SYS1 = SYS1 + RMSD1.'*GLM_t.'+RMSBETTA^(-1)*RMSD0.'*GLML_t.';
end

% F constraints
SYS2 = RMSC0*hRMSy_tF.'+RMSC1*hRMSy_t.'+RMSC2*RMScsi_tF.';

% G constraints
if nG>0
    SYS3 = RMSD0*hRMSy_tF.'+RMSD1*hRMSy_t.'+RMSD2*RMScsi_t.';
else
    SYS3 = [];
end

% Law of motion of exogenous variables
SYS4 = RMScsi_tF.'-RMSS*RMScsi_t.'-RMSeps_tF.';

% Identification of artificial variables
SYS5 = hRMSyL_tF.'-hRMSy_t.';
SYS6 = RMScsiL_tF.'-RMScsi_t.';
SYS7 = GLML_tF.'-GLM_t.';

SYS = vpa([SYS1;SYS2;SYS3;SYS4;SYS5;SYS6;SYS7]);

G0 = -jacobian(SYS,RMSk_tF);
G1 = jacobian(SYS,RMSk_t);
G2 = jacobian(SYS,RMSeps_tF);

G3 = vpa(eye(ny+ncsi+nF+nG+ny+ncsi+nG));
G3(:,ny+1:ny+nF) = [];
G3(:,ny+nG+1:end) = [];

% % Simplify the system by checking for expectation terms:
% %   In the equations that might be forward looking but not necessarily 
% %   (FOC only, the G constraints have to be forward looking)
% %   search for tt variables and if not push everything one period and erase
% %   expectation error column.
% % FOC terms
% cv = find(all(G0(1:ny,:)==0,2)~=0);
% G0(cv,:) = -G1(cv,:);
% G1(cv,:) = 0;
% G3(:,cv) = [];

%% ------------------------------------------------------------------------

%% generate output structure

RMSmat.A0 = vpa(RMSA0);
RMSmat.A1 = vpa(RMSA1);
RMSmat.B0 = vpa(RMSB0);
RMSmat.B1 = vpa(RMSB1);
RMSmat.C0 = vpa(RMSC0);
RMSmat.C1 = vpa(RMSC1);
RMSmat.C2 = vpa(RMSC2);
RMSmat.D0 = vpa(RMSD0);
RMSmat.D1 = vpa(RMSD1);
RMSmat.D2 = vpa(RMSD2);

RMSmat.G0 = vpa(G0);
RMSmat.G1 = vpa(G1);
RMSmat.G2 = vpa(G2);
RMSmat.G3 = vpa(G3);

RMSmat.S = RMSS;
RMSmat.beta = RMSBETTA;

% RMSmat.U_y = vpa(RMSU_y);
% RMSmat.U_yy = vpa(RMSU_yy);
% RMSmat.U_ycsi = vpa(RMSU_ycsi);
% 
% RMSmat.F_y = vpa(RMSF_y);
% RMSmat.F_csi = vpa(RMSF_csi);
% RMSmat.F_yb = vpa(RMSF_yb);
% RMSmat.F_yy = vpa(RMSF_yy);
% RMSmat.F_ycsi = vpa(RMSF_ycsi);
% RMSmat.F_yby = vpa(RMSF_yby);
% RMSmat.F_ybcsi = vpa(RMSF_ybcsi);
% RMSmat.F_ybyb = vpa(RMSF_ybyb);
% 
% if nG>0
%     RMSmat.G_y = vpa(RMSG_y);
%     RMSmat.G_csi = vpa(RMSG_csi);
%     RMSmat.G_yf = vpa(RMSG_yf);
%     RMSmat.G_yy = vpa(RMSG_yy);
%     RMSmat.G_ycsi = vpa(RMSG_ycsi);
%     RMSmat.G_yyf = vpa(RMSG_yyf);
%     RMSmat.G_yfcsi = vpa(RMSG_yfcsi);
%     RMSmat.G_yfyf = vpa(RMSG_yfyf);
% end

%% ------------------------------------------------------------------------

