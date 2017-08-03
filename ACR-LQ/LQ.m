function [LQmat,LQk_t]=LQ(LQy_t,LQcsi_t,LQU,LQF,LQG,LQBETTA,LQS,LQy_ss,FLM_ss,GLM_ss,yLogIdx)

% LQ
%
% This function uses the LQ approximation method proposed by Benigno and
% Woodford (2008) to solve for optimal policy.
% 
%   [LQmat,LQk_t]=LQ(LQy_t,LQcsi_t,LQU,LQF,LQG,LQBETTA,LQS,LQy_ss,...
%       LQFLM_ss,LQGLM_ss,yLogIdx)
%
% Inputs:
%   - LQy_t: symbolic array with endogenous state variables
%   - LQcsi_t: symbolic array with exogenous variables
%   - LQU: symbolic period utility function to be maximized by authorities
%   - LQF: symbolic system of possibly backward looking constraints
%   - LQG: symbolic system of forward looking constraints
%   - LQBETTA: symbolic expression for the time discount factor
%   - LQS: symbolic autocorrelation matrix for exogenous variables
%   - LQy_ss: symbolic array of steady state endogenous state variables
%   - FLM_ss: symbolic array of steady state lagrange multipliers of F
%     constraints
%   - GLM_ss: symbolic array of steady state lagrange multipliers of G
%     constraints
%   - yLogIdx: [optional] index of endogenous state variables to be 
%     converted to logs (binary vector of same size of y_t). If ommited all
%     endogenous variables are log-linearized.
%
% Outputs:
%   - LQmat: structure with symbolic matrices used in the LQ welfare 
%     approximation:
%           V_t = y_t'*A(L)*y_t + y_t'*B(L)*csi_t
%       * LQmat.A0, LQmat.A1 are the matrices such that:
%           A(L) = LQmat.A0 + LQmat.AL1*L
%       * LQmat.B0, LQmat.B1, LQmat.B2 are the matrices such that:
%         	B(L) = LQmat.B0 + LQmat.B1*L + LQmat.B2*L^2
%       * LQmat.C0, LQmat.C1 are the matrices such that:
%           C(L) = LQmat.C0 + LQmat.C1*L
%       * LQmat.C2 is also provided in order to allow for the
%         reconstruction of the linearized F constraints:
%           0 = LQmat.C0*ly_tF+LQmat.C1*ly_t-LQmat.C2*csi_tF
%       * LQmat.D0, LQmat.D1 are the matrices such that:
%           D(L) = LQmat.D0 + LQmat.D1*L
%       * LQmat.D2 is also provided in order to allow for the
%         reconstruction of the linearized G constraints:
%           0 = LQmat.D0*E_t[ly_tF]+LQmat.D1*ly_t-LQmat.D2*csi_t
%       * LQmat.G0, LQmat.G1, LQmat.G2, LQmat.G3 
%         are the matrices such that
%           G0*k_tF = G1*k_t + G2*eps_tF + G3*eta_tF
%         with
%           k_t = (hy_t,csi_t,FLM_t,GLM_t,hyL_t,csiL_t,GLML_t)
%	- LQmat.S: autocorrelation matrix for shocks
%	- LQmat.beta: time discount factor
% 
% The model must be written in the following form:
% 	max sum_0_T { beta^t * U(y_t,csi_t) }
% subject to
% 	0 = F(y_t,csi_t,y_tL)
% 	0 = E_t[ g(y_t,csi_t,y_tF) ]
% Furthermore we also require that 
% 	csi_tF = S * csi_t + eps_tF
%
% Required Matlab routines:
%   - Symbolic Toolbox
%
% See also:
% LQ, LQSolveREE, LQCheckSOC, LQCheckSOCold, LQGenSymVar, LQAltRule, 
% LQWEval, MonFrictions, MonFrictionsIRFPlot, MonFrictionsIRFPlotComp, 
% MonFrictionsIRFPlotAltRule, SmetsWouters
%   
% .........................................................................
%
% Copyright 2004-2009 by Filipo Altissimo, Vasco Curdia and Diego Rodriguez 
% Palenzuela.
% Created: July 26, 2004
% Updated: March 2, 2010

% -------------------------------------------------------------------------

% The previous information above can be accessed issuing the following
% command:
%    help LQ
% or
%    doc LQ

%% ------------------------------------------------------------------------

%% Setup some background information

% vpa inputs
inputList = {'LQy_t','LQcsi_t','LQU','LQF','LQG','LQBETTA','LQS','LQy_ss','FLM_ss','GLM_ss'};
for j=1:length(inputList),eval([inputList{j},'=vpa(',inputList{j},');']),end

% sizes of variables/functions
ny=size(LQy_t,2);
ncsi=size(LQcsi_t,2);
nF=size(LQF,1);
nG=size(LQG,1);
LQcsi_ss = sym(zeros(1,ncsi));

% generate tt and t_1 variables
LQy_tF = LQy_t;
LQy_tL = LQy_t;
for j=1:ny
    LQy_tF = subs(LQy_tF,LQy_tF(j),[char(LQy_tF(j)),'F']);
    LQy_tL = subs(LQy_tL,LQy_tL(j),[char(LQy_tL(j)),'L']);
end

%% ------------------------------------------------------------------------

%% Substitute variables with logs

% check if yLogIdx was not specified
if nargin<11
    yLogIdx = true(size(LQy_t));
end

% create log-variables
hLQy_t = LQy_t;
hLQy_tF = LQy_tF;
hLQy_tL = LQy_tL;
hLQy_ss = LQy_ss;
for j=1:ny
    if yLogIdx(j)
        hLQy_t = subs(hLQy_t,LQy_t(j),['h' char(LQy_t(j))]);
        hLQy_tF = subs(hLQy_tF,LQy_tF(j),['h' char(LQy_tF(j))]);
        hLQy_tL = subs(hLQy_tL,LQy_tL(j),['h' char(LQy_tL(j))]);
        hLQy_ss(j) = log(LQy_ss(j));
        % Plug in transformed variables into the model
        LQF = subs(LQF, [LQy_t(j),LQy_tL(j)], exp([hLQy_t(j),hLQy_tL(j)]));
        LQG = subs(LQG, [LQy_t(j),LQy_tF(j)], exp([hLQy_t(j),hLQy_tF(j)]));
        LQU = subs(LQU, [LQy_t(j)], exp([hLQy_t(j)]));
    end
end

%% ------------------------------------------------------------------------

%% Compute derivatives of U

LQU_y = vpa(jacobian(LQU,hLQy_t));
LQU_yy = vpa(jacobian(LQU_y,hLQy_t));
LQU_ycsi = vpa(jacobian(LQU_y,LQcsi_t));

idxSubs = find(LQU_y~=0);
LQU_y(idxSubs) = subs(LQU_y(idxSubs),[hLQy_t,LQcsi_t],[hLQy_ss,LQcsi_ss]);

idxSubs = find(LQU_yy~=0);
LQU_yy(idxSubs) = subs(LQU_yy(idxSubs),[hLQy_t,LQcsi_t],[hLQy_ss,LQcsi_ss]);

idxSubs = find(LQU_ycsi~=0);
LQU_ycsi(idxSubs) = subs(LQU_ycsi(idxSubs),[hLQy_t,LQcsi_t],[hLQy_ss,LQcsi_ss]);

%% Compute derivatives of F

LQF_y = vpa(jacobian(LQF,hLQy_t));
LQF_csi = vpa(jacobian(LQF,LQcsi_t));
LQF_yb = vpa(jacobian(LQF,hLQy_tL));

LQF_yy = vpa(jacobian(LQF_y(:),hLQy_t));
LQF_ycsi = vpa(jacobian(LQF_y(:),LQcsi_t));
LQF_yyb = vpa(jacobian(LQF_y(:),hLQy_tL));
LQF_ybcsi = vpa(jacobian(LQF_yb(:),LQcsi_t));
LQF_ybyb = vpa(jacobian(LQF_yb(:),hLQy_tL));

idxSubs = find(LQF_y~=0);
LQF_y(idxSubs) = subs(LQF_y(idxSubs),[hLQy_t,LQcsi_t,hLQy_tL],[hLQy_ss,LQcsi_ss,hLQy_ss]);

idxSubs = find(LQF_csi~=0);
LQF_csi(idxSubs) = subs(LQF_csi(idxSubs),[hLQy_t,LQcsi_t,hLQy_tL],[hLQy_ss,LQcsi_ss,hLQy_ss]);

idxSubs = find(LQF_yb~=0);
LQF_yb(idxSubs) = subs(LQF_yb(idxSubs),[hLQy_t,LQcsi_t,hLQy_tL],[hLQy_ss,LQcsi_ss,hLQy_ss]);

idxSubs = find(LQF_yy~=0);
LQF_yy(idxSubs) = subs(LQF_yy(idxSubs),[hLQy_t,LQcsi_t,hLQy_tL],[hLQy_ss,LQcsi_ss,hLQy_ss]);
LQF_yy = reshape(LQF_yy,nF,ny,ny);

idxSubs = find(LQF_ycsi~=0);
LQF_ycsi(idxSubs) = subs(LQF_ycsi(idxSubs),[hLQy_t,LQcsi_t,hLQy_tL],[hLQy_ss,LQcsi_ss,hLQy_ss]);
LQF_ycsi = reshape(LQF_ycsi,nF,ny,ncsi);

idxSubs = find(LQF_yyb~=0);
LQF_yyb(idxSubs) = subs(LQF_yyb(idxSubs),[hLQy_t,LQcsi_t,hLQy_tL],[hLQy_ss,LQcsi_ss,hLQy_ss]);
LQF_yyb = reshape(LQF_yyb,nF,ny,ny);

idxSubs = find(LQF_ybcsi~=0);
LQF_ybcsi(idxSubs) = subs(LQF_ybcsi(idxSubs),[hLQy_t,LQcsi_t,hLQy_tL],[hLQy_ss,LQcsi_ss,hLQy_ss]);
LQF_ybcsi = reshape(LQF_ybcsi,nF,ny,ncsi);

idxSubs = find(LQF_ybyb~=0);
LQF_ybyb(idxSubs) = subs(LQF_ybyb(idxSubs),[hLQy_t,LQcsi_t,hLQy_tL],[hLQy_ss,LQcsi_ss,hLQy_ss]);
LQF_ybyb = reshape(LQF_ybyb,nF,ny,ny);

%% Compute derivatives of G

if nG>0
    LQG_y = vpa(jacobian(LQG,hLQy_t));
    LQG_csi = vpa(jacobian(LQG,LQcsi_t));
    LQG_yf = vpa(jacobian(LQG,hLQy_tF));

    LQG_yy = vpa(jacobian(LQG_y(:),hLQy_t));
    LQG_ycsi = vpa(jacobian(LQG_y(:),LQcsi_t));
    LQG_yfy = vpa(jacobian(LQG_yf(:),hLQy_t));
    LQG_yfcsi = vpa(jacobian(LQG_yf(:),LQcsi_t));
    LQG_yfyf = vpa(jacobian(LQG_yf(:),hLQy_tF));

    idxSubs = find(LQG_y~=0);
    LQG_y(idxSubs) = subs(LQG_y(idxSubs),[hLQy_t,LQcsi_t,hLQy_tF],[hLQy_ss,LQcsi_ss,hLQy_ss]);

    idxSubs = find(LQG_csi~=0);
    LQG_csi(idxSubs) = subs(LQG_csi(idxSubs),[hLQy_t,LQcsi_t,hLQy_tF],[hLQy_ss,LQcsi_ss,hLQy_ss]);

    idxSubs = find(LQG_yf~=0);
    LQG_yf(idxSubs) = subs(LQG_yf(idxSubs),[hLQy_t,LQcsi_t,hLQy_tF],[hLQy_ss,LQcsi_ss,hLQy_ss]);

    idxSubs = find(LQG_yy~=0);
    LQG_yy(idxSubs) = subs(LQG_yy(idxSubs),[hLQy_t,LQcsi_t,hLQy_tF],[hLQy_ss,LQcsi_ss,hLQy_ss]);
    LQG_yy = reshape(LQG_yy,nG,ny,ny);

    idxSubs = find(LQG_ycsi~=0);
    LQG_ycsi(idxSubs) = subs(LQG_ycsi(idxSubs),[hLQy_t,LQcsi_t,hLQy_tF],[hLQy_ss,LQcsi_ss,hLQy_ss]);
    LQG_ycsi = reshape(LQG_ycsi,nG,ny,ncsi);

    idxSubs = find(LQG_yfy~=0);
    LQG_yfy(idxSubs) = subs(LQG_yfy(idxSubs),[hLQy_t,LQcsi_t,hLQy_tF],[hLQy_ss,LQcsi_ss,hLQy_ss]);
    LQG_yfy = reshape(LQG_yfy,nG,ny,ny);

    idxSubs = find(LQG_yfcsi~=0);
    LQG_yfcsi(idxSubs) = subs(LQG_yfcsi(idxSubs),[hLQy_t,LQcsi_t,hLQy_tF],[hLQy_ss,LQcsi_ss,hLQy_ss]);
    LQG_yfcsi = reshape(LQG_yfcsi,nG,ny,ncsi);

    idxSubs = find(LQG_yfyf~=0);
    LQG_yfyf(idxSubs) = subs(LQG_yfyf(idxSubs),[hLQy_t,LQcsi_t,hLQy_tF],[hLQy_ss,LQcsi_ss,hLQy_ss]);
    LQG_yfyf = reshape(LQG_yfyf,nG,ny,ny);
else
    LQG_y = [];
    LQG_csi = [];
    LQG_yf = [];

    LQG_yy = [];
    LQG_ycsi = [];
    LQG_yfy = [];
    LQG_yfcsi = [];
    LQG_yfyf = [];
end

%% note about the derivatives
% all second order derivatives are stored in 3-dim arrays, corresponding to
%   dim_function x dim_vars_D1 x dim_vars_D2
% This means that for each call of squeeze(...(j,:,:)) yields
%   dim_vars_D1 x dim_vars_D2

%% ------------------------------------------------------------------------

%% Compute LQ approximation matrices

LQR = vpa(zeros(ny));
LQH = vpa(zeros(ny));
LQZ0 = vpa(zeros(ny,ncsi)); LQZ1 = LQZ0; LQZ2 = LQZ0;
for j=1:nG
    LQR = LQR + vpa(LQBETTA^(-1)*GLM_ss(j)*squeeze(LQG_yfy(j,:,:)));
    LQH = LQH + vpa(GLM_ss(j)*(squeeze(LQG_yy(j,:,:))+...
                LQBETTA^(-1)*squeeze(LQG_yfyf(j,:,:))));
    LQZ1 = LQZ1 + vpa(GLM_ss(j)*reshape(LQG_ycsi(j,:,:),ny,ncsi));
    LQZ2 = LQZ2 + vpa(LQBETTA^(-1)*GLM_ss(j)*reshape(LQG_yfcsi(j,:,:),ny,ncsi));
end
for j=1:nF
    LQR = LQR + vpa(FLM_ss(j)*squeeze(LQF_yyb(j,:,:)));
    LQH = LQH + vpa(FLM_ss(j)*(squeeze(LQF_yy(j,:,:))+...
                LQBETTA*squeeze(LQF_ybyb(j,:,:))));
    LQZ0 = LQZ0 + vpa(LQBETTA*FLM_ss(j)*reshape(LQF_ybcsi(j,:,:),ny,ncsi));
    LQZ1 = LQZ1 + vpa(FLM_ss(j)*reshape(LQF_ycsi(j,:,:),ny,ncsi));
end

LQQ = LQU_yy + LQH;

LQA0 = LQQ;
LQA1 = 2*LQR;

LQB0 = LQZ0;
LQB1 = LQZ1 + LQU_ycsi;
LQB2 = LQZ2;

LQC0 = LQF_y;
LQC1 = LQF_yb;

LQC2 = -LQF_csi;

LQD0 = LQG_yf;
LQD1 = LQG_y;

LQD2 = -LQG_csi;

%% Note on polynomials defined
% At this stage we have the following polynomials defined:
%   A(L) = LQA0 + LQA1*L
%   B(L) = LQB0 + LQB1*L + LQB2*L^2
%   C(L) = LQC0 + LQC1*L
%   D(L) = LQD0 + LQD1*L
%   Z(L) = LQZ0 + LQZ1*L + LQZ2*L^2 

%% ------------------------------------------------------------------------

%% Generate variables for the state space
% We now need to set up the state space system in the following form:
%   G0*k_tF.' = G1*k_t.' + G3*eps_tF.' + G4*eta_tF.'
% with
%   k_t = (hy_t,csi_t,FLM_t,GLM_t,hyL_t,csiL_t,GLML_t)

hLQyL_t=hLQy_t; hLQyL_tF=hLQyL_t;
for j=1:ny
    hLQyL_t = subs(hLQyL_t,hLQyL_t(j),strrep(char(hLQyL_t(j)),'_t','L_t'));
    hLQyL_tF = subs(hLQyL_tF,hLQyL_tF(j),strrep(char(hLQyL_tF(j)),'_t','L_tF'));
end

LQcsi_tF=LQcsi_t; LQcsiL_t=LQcsi_t; LQcsiL_tF=LQcsi_t; LQeps_tF=LQcsi_t;
for j=1:ncsi
    LQcsi_tF = subs(LQcsi_tF,LQcsi_tF(j),[char(LQcsi_tF(j)),'t']);
    LQcsiL_t = subs(LQcsiL_t,LQcsiL_t(j),strrep(char(LQcsiL_t(j)),'_t','L_t'));
    LQcsiL_tF = subs(LQcsiL_tF,LQcsiL_tF(j),strrep(char(LQcsiL_tF(j)),'_t','L_tF'));
    LQeps_tF = subs(LQeps_tF,LQeps_tF(j),['eps_',char(LQcsi_tF(j)),'t']);
end

FLM_t=[]; FLM_tF=[];
for j=1:nF
    namej = char(FLM_ss(j)); namej = namej(1:end-3);
    eval(['syms ',namej,'_t ',namej,'_tF'])
    FLM_t = [FLM_t, eval([namej,'_t'])];
    FLM_tF = [FLM_tF, eval([namej,'_tF'])];
end

GLM_t=[]; GLM_tF=[];GLML_t=[]; GLML_tF=[];
for j=1:nG
    namej = char(GLM_ss(j)); namej = namej(1:end-3);
    eval(['syms ',namej,'_t ',namej,'_tF'])
    GLM_t = [GLM_t, eval([namej,'_t'])];
    GLM_tF = [GLM_tF, eval([namej,'_tF'])];
    eval(['syms ',namej,'L_t ',namej,'L_tF'])
    GLML_t = [GLML_t, eval([namej,'L_t'])];
    GLML_tF = [GLML_tF, eval([namej,'L_tF'])];
end

LQk_t = [hLQy_t, LQcsi_t, FLM_t, GLM_t, hLQyL_t, LQcsiL_t, GLML_t];
LQk_tF = [hLQy_tF, LQcsi_tF, FLM_tF, GLM_tF, hLQyL_tF, LQcsiL_tF, GLML_tF];

%% Generate the state space system with the FOC

% FOC
SYS1 = 1/2*(LQQ+LQQ.')*hLQy_t.'+LQR*hLQyL_t.'+LQBETTA*LQR.'*hLQy_tF.'+...
        LQB0*LQcsi_tF.'+LQB1*LQcsi_t.'+LQB2*LQcsiL_t.'+...
        LQC0.'*FLM_t.'+LQBETTA*LQC1.'*FLM_tF.';
if nG>0
    SYS1 = SYS1 + LQBETTA^(-1)*LQD0.'*GLML_t.'+LQD1.'*GLM_t.';
end

% F constraints
SYS2 = LQC0*hLQy_tF.'+LQC1*hLQy_t.'+LQF_csi*LQcsi_tF.';

% G constraints
if nG>0
    SYS3 = LQD0*hLQy_tF.'+LQD1*hLQy_t.'+LQG_csi*LQcsi_t.';
else
    SYS3 = [];
end

% Law of motion of exogenous variables
SYS4 = LQcsi_tF.'-LQS*LQcsi_t.'-LQeps_tF.';

% Identification of artificial variables
SYS5 = hLQyL_tF.'-hLQy_t.';
SYS6 = LQcsiL_tF.'-LQcsi_t.';
SYS7 = GLML_tF.'-GLM_t.';

SYS = vpa([SYS1;SYS2;SYS3;SYS4;SYS5;SYS6;SYS7]);

G0 = -jacobian(SYS,LQk_tF);
G1 = jacobian(SYS,LQk_t);
G2 = jacobian(SYS,LQeps_tF);

G3 = eye(ny+ncsi+nF+nG+ny+ncsi+nG);
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

LQmat.A0 = vpa(LQA0);
LQmat.A1 = vpa(LQA1);
LQmat.B0 = vpa(LQB0);
LQmat.B1 = vpa(LQB1);
LQmat.B2 = vpa(LQB2);
LQmat.C0 = vpa(LQC0);
LQmat.C1 = vpa(LQC1);
LQmat.C2 = vpa(LQC2);
LQmat.D0 = vpa(LQD0);
LQmat.D1 = vpa(LQD1);
LQmat.D2 = vpa(LQD2);

LQmat.G0 = vpa(G0);
LQmat.G1 = vpa(G1);
LQmat.G2 = vpa(G2);
LQmat.G3 = vpa(G3);

LQmat.S = LQS;
LQmat.beta = LQBETTA;

LQmat.U_y = vpa(LQU_y);
LQmat.U_yy = vpa(LQU_yy);
LQmat.U_ycsi = vpa(LQU_ycsi);

LQmat.F_y = vpa(LQF_y);
LQmat.F_csi = vpa(LQF_csi);
LQmat.F_yb = vpa(LQF_yb);
LQmat.F_yy = vpa(LQF_yy);
LQmat.F_ycsi = vpa(LQF_ycsi);
LQmat.F_yyb = vpa(LQF_yyb);
LQmat.F_ybcsi = vpa(LQF_ybcsi);
LQmat.F_ybyb = vpa(LQF_ybyb);

if nG>0
    LQmat.G_y = vpa(LQG_y);
    LQmat.G_csi = vpa(LQG_csi);
    LQmat.G_yf = vpa(LQG_yf);
    LQmat.G_yy = vpa(LQG_yy);
    LQmat.G_ycsi = vpa(LQG_ycsi);
    LQmat.G_yfy = vpa(LQG_yfy);
    LQmat.G_yfcsi = vpa(LQG_yfcsi);
    LQmat.G_yfyf = vpa(LQG_yfyf);
end

%% ------------------------------------------------------------------------

