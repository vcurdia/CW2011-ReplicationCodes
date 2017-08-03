function REE=LQAltRule(Rule,LQmat,LQS,LQz_t,y_ss,yLogIdx)

% LQAltRule
%
% Generates state space and REE solution for alternative rule.
%
%   REE=LQAltRule(Rule,LQmat,LQS,LQz_t,y_ss,yLogIdx)
%
% Inputs:
%   - Rule: symbolic vector with equations defining the rule(s)
%     [See below remarks about rule]
%   - LQmat: structure with numerical matrices for the LQ approx
%   - LQS: numeric autocorrelation matrix for exogenous variables
%   - LQz_t: symbolic vector with symbolic variables of state space
%   - y_ss: vector with steady state levels of endogenous variables
%   - yLogIdx: [optional] index of endogenous state variables to be 
%     converted to logs (binary vector of same size of y_t). If ommited all
%     endogenous variables are log-linearized.
%
% Outputs:
%   - REE: structure with fields:
%     * REE.Phi1 and REE.Phi2: matrices for the REE solution such that:
%       z_t = Phi1*z_tL + Phi2*eps_t
%     * REE.z_t: list of symbolic variables used in the rule
%
% The code will issue a warning message if the solution to the REE is not
% normal in any way.
%
% Remarks about rule(s):
% ----------------------
%
% This code is designed to receive as an input non-linear rule(s), based on
% the initial variables introduced in LQ.
% 
% If linear rule(s) is desired, simply write the rule in multiplicative 
% terms, with exponents corresponding to the linear slope parameters.
% 
% General form of rule:
%   Rule(y_t,y_tF,y_tL,csi_t)=0
%
% Rule can contain forward looking and backward looking terms, but not both
% at the same time. If both are desired then artificial variables should
% have been created in the original code (for this example: MonFrictions)
% within the y_t vector of endogenous variables, and additional
% corresponding F constraints. 
% e.g.:
% Want to use both future and past levels of inflation. Proceed as follows:
% 1. In the MonFrictions code add to the list of variables 'y', the
%    following artificial variable: LPi (guess value should be the same as
%    for Pi)
% 2. In the MonFrictions code add additional lagrange multiplier (in FLM 
%    vector): FLMLPI (also add the guess value to be, e.g. 0)
% 3. In MonFrictions insert constraint in the F set of constraints:
%       F(pos) = LPi_t - Pi_tL;
%    where pos should be the position in the F(:) vector of the new
%    constraint, e.g. previous nF+1
% 4. In this code, for the alternative rule, write:
%       Rule = R_t - Pi_tF^phi_1 * Pi_t^phi_2 * LPi_t^phi_3;
% 5. This code will then analyze the expression and proceed accordingly.
%
% NOTE: This rule should be formulated so that it is consistent with the 
%       optimal steady state 
%
% Convention for refering to variables in equations:
%   x_t   refers to x{t}
%   x_tF  refers to x{t+1}
%   x_tL refers to x{t-1}
%   x_ss  refers to the steady state level of 'x'
%   (where x refers to some variable with name 'x')
%
% NOTE: the rule can't use any variables not defined in the base code,
%       including exogenous/shock variables. If so desired, insert them in
%       the initial setup right from the beginning.
% 
% Required m-files:
%   - symbolic toolbox
%   - gensys, available in Chris Sims's website
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
% Updated: July 27, 2009

% -------------------------------------------------------------------------

% The previous information above can be accessed issuing the following
% command:
%    help LQAltRule
% or
%    doc LQAltRule

%% ------------------------------------------------------------------------

%% Define some variables
nR = length(Rule);
nz = length(LQz_t);
ny = size(LQmat.A0,1);
nF = size(LQmat.C0,1);
nG = size(LQmat.D0,1);
ncsi = size(LQmat.B0,2);
if nargin<6, yLogIdx = true(ny); end % check if yLogIdx was not specified


%% check number of rules
if (nR-(ny-nF-nG))~=0
    error('Number of rules is %.0f, but number of rules needed is %.0f!',nR,ny-nF-nG)
end

%% ------------------------------------------------------------------------

%% Define z_tF and z_tL in case they are needed
LQz_tF = LQz_t; LQz_tL = LQz_t;
for j=1:nz
    LQz_tF(j) = subs(LQz_t(j),LQz_t(j),[char(LQz_t(j)),'F']);
    LQz_tL(j) = subs(LQz_t(j),LQz_t(j),[char(LQz_t(j)),'L']);
end
hy_t = LQz_t(1:ny);
hy_tF = LQz_tF(1:ny);
hy_tL = LQz_tL(1:ny);

%% Define original y vectors
y_t = hy_t; y_tF = hy_tF; y_tL = hy_tL;
for j=1:ny
    yname = char(y_t(j)); 
    yname = yname(2:end); sym(yname); y_t(j) = yname;
    yname = [yname,'t']; sym(yname); y_tF(j) = yname;
    yname = [yname(1:end-1),'_1']; sym(yname); y_tL(j) = yname;
end

%% change variables to exp(log_variables)
logRule = Rule;
hy_ss = y_ss;
for j=1:ny
    if yLogIdx(j)
        hy_ss(j) = log(y_ss(j));
        logRule = subs(logRule, [y_t(j),y_tF(j),y_tL(j)],exp([hy_t(j),hy_tF(j),hy_tL(j)]));
    end
end

%% define csi_t vector
csi_t = LQz_t(ny+(1:ncsi));
csi_ss = zeros(1,ncsi);

%% define LQzz_t vectors
% zz refers to the (y,csi) vector
nzz = ny+ncsi;
LQzz_t = LQz_t(1:nzz);
LQzz_tF = LQz_tF(1:nzz);
LQzz_tL = LQz_tL(1:nzz);

%% ------------------------------------------------------------------------

%% Compute derivatives of logRule
derivs = {'LQzz_t','LQzz_tF','LQzz_tL','csi_t'};
for j=1:length(derivs)
    jd = derivs{j};
    Namejd = ['logRule_',jd];
    eval([Namejd,' = jacobian(logRule,',jd,');'])
    eval([Namejd,' = double(subs(',Namejd,',[hy_t,hy_tF,hy_tL,csi_t],[hy_ss,hy_ss,hy_ss,csi_ss]));'])
end

%% generate state space with rule

% introduce the constraints
%   0 = LQmat.D0*E_t[hy_tF]+LQmat.D1*hy_t-LQmat.D2*csi_t
%   0 = LQmat.C0*hy_t+LQmat.C1*hy_tL-LQmat.C2*csi_t
G0 = [-LQmat.D0 zeros(nG,ncsi); -LQmat.C0 LQmat.C2; zeros(ncsi,ny) eye(ncsi)];
G1 = [LQmat.D1 -LQmat.D2; LQmat.C1 zeros(nF,ncsi); zeros(ncsi,ny) LQS];
G2 = [zeros(nG+nF,ncsi); eye(ncsi);zeros(nR,ncsi)];
G3 = [eye(nG);zeros(nF+ncsi,nG)];

% in this framework the constant is zero
const = zeros(nzz,1);

for jr=1:nR
    % check rule type
    if ~all(logRule_LQzz_tF(jr,:)==0)
        if ~all(logRule_LQzz_tL(jr,:)==0)
            error('Cannot have both forward and backward looking variables in same rule!')
        end
        G0 = [G0; -logRule_LQzz_tF(jr,:)];
        G1 = [G1; logRule_LQzz_t(jr,:)];
        G3 = [G3, zeros(size(G3,1),1); zeros(1,size(G3,2)) 1];
    else
        G0 = [G0; -logRule_LQzz_t(jr,:)];
        G1 = [G1; logRule_LQzz_tL(jr,:)];
        G3 = [G3; zeros(1,size(G3,2))];        
    end
end

%% ------------------------------------------------------------------------

%% Solve for REE solution 
[Phi1,Const,Phi2,fmat,fwt,ywt,gev,eu] = gensys(G0,G1,const,G2,G3);
if any(eu~=1)
if any(eu~=1)
    fprintf('Warning: eu = (%.0f,%.0f)',eu)
    if all(eu)==-2
        fprintf(' - Coincident zeros!!!')
    elseif eu(1)~=1
        fprintf(' - Solution does not exist!!!')
    elseif eu(2)~=1
        fprintf(' - Solution is not unique!!!')
    end
    fprintf('\n')
end
end

%% ------------------------------------------------------------------------

%% Adjust matrices to bigger z vector
% if lagrange multipliers are part of z vector
if nz>nzz
    Phi1 = [Phi1 zeros(nzz,nG); zeros(nG,nz)];
    Phi2 = [Phi2; zeros(nG,ncsi)];
end    

%% ------------------------------------------------------------------------

%% Prepare output
REE.Phi1 = Phi1;
REE.Phi2 = Phi2;
REE.z_t = LQz_t;
REE.eu = eu;
REE.Sy = [eye(ny),zeros(ny,nzz-ny)];
REE.Scsi = [zeros(ncsi,ny),eye(ncsi),zeros(ncsi,nzz-ny-ncsi)];

%% ------------------------------------------------------------------------

