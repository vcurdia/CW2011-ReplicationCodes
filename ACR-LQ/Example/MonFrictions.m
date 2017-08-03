% MonFrictions
%
% Model with monetary frictions into the setup of Benigno and Woodford (2004)
%
% This is one example of how to use the LQ codes. Sections are identified
% with lines starting with '%%'. Pay attention to the comments in each
% section, especially those mentioning that you shouldn't change that 
% section unless you know what you are doing.
%
% Convention for refering to variables in equations:
%   x_t  refers to x{t}
%   x_tF refers to x{t+1}
%   x_tL refers to x{t-1}
%   x_ss  refers to the steady state level of 'x'
%   (where x refers to some variable with name 'x')
%
% Required m-files:
%   - symbolic toolbox
%   - LQGenSymVar
%   - LQ
%   - LQSolveREE
%   - LQCheckSOC
%   - gensys.m, available in Chris Sims's website
%   - csolve.m, available in Chris Sims's website
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
% Created: December 15, 2004
% Updated: March 30, 2010

% -------------------------------------------------------------------------

% The previous information above can be accessed issuing the following
% command:
%    help MonFrictions
% or
%    doc MonFrictions

%% ------------------------------------------------------------------------

%% preamble

clear all
tic

NumPrecision = 1e-10;
SolveIterations = 10000;

%% ------------------------------------------------------------------------

%% define parameters
BETTA=.99;
ALPHA=0.75;
LAMBDA=0.7;
ISIGMA_T=0.157;
OMEGA=0.473;
PHI=1;
NU = (OMEGA+1)/PHI-1;
THETA=10;
a=0.0111; % according to Schmitt-Grohe and Uribe (2004)
b=0.07524; % according to Schmitt-Grohe and Uribe (2004)
n_min = sqrt(b/a); 
mu_ss=1;
A_ss=1;
BARH_ss=1;
BARC_ss=1;
s_g=0.3;
RHO_CSI=0.7;
tau_ss=-1/(THETA-1);%0.2;%-1/(THETA-1); 
PSI=0;

%% define endogenous variables
% The first column includes the names
% The second column includes the guess values for the steady state
y = {'R', 1.010101;
     'Y', 1.142937;
     'n', 2.603532;
     'm', 0.307296;
     'Pi', 1;
     'F', 3.677436;
     'K', 3.677436;
     'lambda', 1.035643;
     'Delta', 1};
 
%% Identify which variables are to be log-linearized
% OPTIONAL -- comment this section if all variables are to be log-linearized. 
% Define a binary vector with the same length as y to identify which
% variables are to be loglinearized.

% yLogIdx = true(1,size(y,1));

%% define exogenous variables
% The list includes the names
% NOTE: Exogenous variables must be defined such that they are zero in 
%       steady state
csi = {'htau','hmu','hG','hBARC','hBARH','hA'};
% S is the matrix such that:
%   csi{t} = S*csi{t-1} + epsilon{t)
S = RHO_CSI*eye(length(csi)); 

%% define lagrange multipliers for F constraints
% The first column includes the names
% The second column includes the guess values for the steady state
% NOTE: Make sure to include one and only one lagrange multiplier for each
%       equation in F
FLM = {'FLM1', 0.458954;
       'FLM2', 0.182565;
       'FLM3', 0;
       'FLM4', -0.461469;
       'FLM5', -2.246906};

%% define lagrange multipliers for G constraints
% The first column includes the names
% The second column includes the guess values for the steady state
% NOTE: Make sure to include one and only one lagrange multiplier for each
%       equation in G
GLM = {'GLM1', 0;
       'GLM2', 0.501946;
       'GLM3', -0.501946};

%% Generate symbolic variables
% DO NOT CHANGE THIS SECTION
% Generate all symbolic variables required to setup the model
% Do not work with symbolic variables and expressions before executing this
% section!
LQGenSymVar

%% ------------------------------------------------------------------------

%% Auxiliary definitions, used in later sections
% this section is specific to this example, and doesn't need to be defined
% in other applications
s_t = a*n_t+b/n_t-2*(a*b)^(1/2);
ds_t = diff(s_t,n_t);
tau_t = 1-(1-tau_ss)*exp(-htau_t);
mu_t = mu_ss*exp(hmu_t);
if s_g==0
    G_t = hG_t;
else
    G_t = s_g*Y_ss*exp(hG_t);
end
BARC_t = BARC_ss*exp(hBARC_t);
BARH_t = BARH_ss*exp(hBARH_t);
A_t = A_ss*exp(hA_t);
U1 = BARC_t^ISIGMA_T/(1-ISIGMA_T)*((Y_t-G_t)/(1+s_t))^(1-ISIGMA_T);
V = LAMBDA/(1+NU)*Y_t^(1+OMEGA)/BARH_t^(-NU)/A_t^(-PHI*(1+NU))*Delta_t;
U1_Y = diff(U1,Y_t);
V_Y = diff(V,Y_t);

%% define utility function
% specify symbolic expression that must be in U, the period utility
U = U1-V;

%% F constraints
% specify symbolic column vector with expressions for the F constraints
F = sym(zeros(min(2,nF),1));
F(1) = (1+s_t)/(1+s_t+ds_t*n_t)*U1_Y-lambda_t;
F(2) = ds_t*n_t^2-(1-PSI)*(R_t-1)/R_t;
F(3) = (Y_t-G_t)/(1+s_t)/n_t-m_t;
F(4) = ((1-ALPHA*Pi_t^(THETA-1))/(1-ALPHA))^((1+OMEGA*THETA)/(THETA-1))-F_t/K_t;
F(5) = ALPHA*Delta_tL*Pi_t^(THETA*(1+OMEGA))...
       +(1-ALPHA)*((1-ALPHA*Pi_t^(THETA-1))/(1-ALPHA))^(THETA*(1+OMEGA)/(THETA-1))...
       -Delta_t;
if length(F)~=nF
    error('number of F constraints different from number of FLM multipliers')
end
    
%% G constraints
% specify symbolic column vector with expressions for the G constraints
G = sym(zeros(min(2,nG),1));
G(1) = BETTA*R_t*lambda_tF/Pi_tF-lambda_t;
G(2) = (1-tau_t)*lambda_t*Y_t+ALPHA*BETTA*Pi_tF^(THETA-1)*F_tF-F_t;
G(3) = THETA/(THETA-1)*mu_t*V_Y*Y_t/Delta_t+ALPHA*BETTA*Pi_tF^(THETA*(1+OMEGA))*K_tF-K_t;
if length(G)~=nG
    error('number of G constraints different from number of GLM multipliers')
end

%% ------------------------------------------------------------------------

%% Compute steady state
% CHANGE THIS SECTION ONLY IF YOU KNOW WHAT YOU ARE DOING AND YOU HAVE AN
% ALTERNATIVE WAY TO COMPUTE THE OPTIMAL STEADY STATE!
% In this section we compute the optimal steady state levels for endogenous
% variables and multipliers.

% generate the steady state system
ssSys = jacobian(U,y_t).'+...
        (FLM_ss*jacobian(F,y_t)).'+BETTA*(FLM_ss*jacobian(F,y_tL)).'+...
        (GLM_ss*jacobian(G,y_t)).'+BETTA^(-1)*(GLM_ss*jacobian(G,y_tF)).';
ssSys = [ssSys; F];
ssSys = [ssSys; G];

% plug in the steady state values in symbolic form
ssSys = subs(ssSys,[y_t,y_tF,y_tL,csi_t],[y_ss,y_ss,y_ss,csi_ss]);

% prepare variables
nSys = length(ssSys);
x = y;
x(ny+1:ny+nF) = FLM;
x(ny+nF+1:ny+nF+nG) = GLM;
nx = length(x);

% generate function for csolve
fid=fopen('ssSysF.m','w');
fprintf(fid,'function f=ssSysF(x) \n');
fprintf(fid,'for j=1:size(x,2)  \n');
for j=1:nSys
    fprintf(fid,[x(j).name,'_ss = x(',int2str(j),',j);  \n']);
end
for j=1:nSys
    fprintf(fid,['f(',int2str(j),',j) = ',char(ssSys(j)),';  \n']);
end
fprintf(fid,'end  \n');
fclose(fid);

% Solve
[x1,rc] = csolve('ssSysF',[x.guess]',[],NumPrecision,SolveIterations);
if rc~=0, error(['Solution of steady state system is not normal, rc = ', int2str(rc)]), end
delete ssSysF.m

% evaluate variables (need to keep different notation, so as not to erase
% symbolic expressions at this stage)
y_ss_1 = x1(1:ny)';
FLM_ss_1 = x1(ny+1:ny+nF)';
GLM_ss_1 = x1(ny+nF+1:ny+nF+nG)';

%% Present steady state
disp(sprintf('\n '))
disp('Steady state results:')
disp('=====================')
disp(' ')
maxnmlength = length(x(1).name);
for j=2:nx
    maxnmlength = max(maxnmlength,length(x(j).name));
end
for j=1:nx
    nmtext = x(j).name;
    ldiff = maxnmlength-length(nmtext);
    if ldiff>0, nmtext = [nmtext, repmat(' ',1,ldiff)]; end
    if sign(x1(j))~=-1, vnntxt = ' '; else vnntxt = ''; end
    disp(sprintf([nmtext,' = ',vnntxt,'%f'],x1(j)))
end
disp(' ')

%% ------------------------------------------------------------------------

%% Run LQ function to get symbolic matrices
% DO NOT CHANGE THIS SECTION
if exist('yLogIdx','var')
    [LQmat,LQk_t]=LQ(y_t,csi_t,U,F,G,BETTA,S,y_ss,FLM_ss,GLM_ss,yLogIdx);
else
    [LQmat,LQk_t]=LQ(y_t,csi_t,U,F,G,BETTA,S,y_ss,FLM_ss,GLM_ss);
end

%% Fill LQ matrices with steady state and parameter vectors
% DO NOT CHANGE THIS SECTION
for j=1:nx
    eval([x(j).name,'_ss = x1(',int2str(j),');'])
end
y_ss = eval(y_ss);
FLM_ss = eval(FLM_ss);
GLM_ss = eval(GLM_ss);
matnames = fieldnames(LQmat);
for j=1:length(matnames)
    LQmat.(matnames{j}) = reshape(eval(LQmat.(matnames{j})(:)),size(LQmat.(matnames{j})));
end

%% Solve the REE
% DO NOT CHANGE THIS SECTION
[REE.LQ,LQz_t]=LQSolveREE(LQmat,LQk_t);

%% Check for the SOC
% DO NOT CHANGE THIS SECTION
SOC = LQCheckSOC(BETTA,S,LQmat);

%% Generate IRFs
% this is optional
nsteps = 25; % can be set to a different number
IRF = zeros([size(REE.LQ.Phi2),nsteps]);
IRF(:,:,1) = REE.LQ.Phi2;
for j=2:nsteps
    IRF(:,:,j) = REE.LQ.Phi1*IRF(:,:,j-1);
end
IRF = permute(IRF,[1,3,2]); % (nz x nstep x ncsi) 

%% ------------------------------------------------------------------------

%% Consider alternative rule
% this is optional

% define rule
phi_pi = 1.5;
Rule = R_t/R_ss-Pi_t^phi_pi; 

% solve REE for alternative rule
if exist('yLogIdx','var')
    REE.Rule=LQAltRule(Rule,LQmat,S,LQz_t,y_ss,yLogIdx);
else
    REE.Rule=LQAltRule(Rule,LQmat,S,LQz_t,y_ss);
end

% get IRF for alternative rule
IRF_Rule = zeros([size(REE.Rule.Phi2),nsteps]);
IRF_Rule(:,:,1) = REE.Rule.Phi2;
for j=2:nsteps
    IRF_Rule(:,:,j) = REE.Rule.Phi1*IRF_Rule(:,:,j-1);
end
IRF_Rule = permute(IRF_Rule,[1,3,2]); % (nz x nstep x ncsi) 

%% Evaluate welfare for each shock
Policies = fieldnames(REE);
nPol = length(Policies);
VarShocks = cell(1,ncsi);
for j=1:ncsi
    VarShocks{j} = diag(ismember(csi,csi{j}));
    for jPol=1:nPol
        W.(csi{j}).(Policies{jPol}) = LQWEval(LQmat,REE.LQ,REE.(Policies{jPol}),VarShocks{j});
    end
    fprintf('\nWelfare comparison for %s:\n',csi{j})
    disp(W.(csi{j}))
end

%% ------------------------------------------------------------------------
%% ------------------------------------------------------------------------

%% Save environment
% Saves the entire workspace after running LQ and generating thed IRFs
% this section is specific to this example. Change at will or even delete/comment
% the following commands are intended to save the workspace in different
% mat files, depending on the parameters. Can be simplified if desired.
if tau_ss~=0.2
    Suffix1 = 'NoDistortions';
else
    Suffix1 = 'Distortions';
end
if PSI
    Suffix2 = 'Cashless';
else
    Suffix2 = 'Frictions';
end
eval(['save MonFrictions_LQ_',Suffix1,Suffix2])

%% Elapsed time
disp(' '), toc, disp(' ')

%% ------------------------------------------------------------------------

