function Exercise3(varargin)

% Exercise 3
%
% Compute the equilibrium under no credit policy and a simple interest rate
% rule. The point is to extract the lagrange multiplier for the zero CB
% lending.
%
% Notation:
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
% Created: February 23, 2010 by Vasco Curdia
% Updated: April 21, 2014 by Vasco Curdia
%
% Copyright 2010-2014 by Vasco Curdia

%% ------------------------------------------------------------------------

%% preamble
% clear all
% tic
ttic=toc();
format short g
NumPrecision = 1e-10;
nsteps = 32;

%% Options
isNoSpread = 0;
isNoDist = 0;
isIgnoreZLB = 0;

BindBp = 0;
dSP = 12;
PersSP = 90;

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


%% Calculation needed
nMax.ZLB = ~isIgnoreZLB*10;

%% Exercise type
ExerciseName = ['Exercise3_dSP_',int2str(dSP),'_Pers_',num2str(PersSP)];
isBind = (BindBp~=0);
if isBind
    ExerciseName = [ExerciseName,'_Bind_',num2str(BindBp),'bp'];
end
if isNoSpread
    ExerciseName = [ExerciseName,'_NoSpread'];
end
if isNoDist
    ExerciseName = [ExerciseName,'_NoDist'];
end
if isIgnoreZLB
    ExerciseName = [ExerciseName,'_IgnoreZLB'];
end
fprintf('\n****%s****',repmat('*',1,length(ExerciseName)))
fprintf('\n*   %s   *',ExerciseName)
fprintf('\n****%s****\n\n',repmat('*',1,length(ExerciseName)))

%% ------------------------------------------------------------------------

%% define parameters
rd_ss = (1.03)^(1/4)-1;
sigmabari = 0.16; %0.16
sigmabar = 1/sigmabari;
phi = 1/0.75;
alpha = 0.66;
omega_y = 0.473;
nu = (omega_y+1)/phi-1;
mu_p = 1.15;
theta = 1+1/(mu_p-1);
spread_ss = (1.02)^(1/4)-1;
omega_ss = ~isNoSpread*spread_ss; 
elast_omega_b = 1/4; 
eta = 1+elast_omega_b*(1+spread_ss)/spread_ss;
varkappa = elast_omega_b*(1+spread_ss)/spread_ss;
delta = 0.975;
pi_b = 0.5;
rho_b = 3.2;
eta_cb = 1;
Yss = 1;
s_c = 0.7;
sigma_bs = 5;
psi = 1;
mu_w_ss = 1;
tau_ss = ~isNoDist*0.2+isNoDist*(1-mu_p*mu_w_ss);
Hbar_ss = 1;
Z_ss = 1;
phi_pi = 2; %1.5; 2
phi_y = 1; %0.5; 1
rho_csi = PersSP/100;

%% Endogenous variables
y = {'Rd','Pi','Y','lambda_b','lambda_s','b','Delta','K','F','omegatil','lcb'};
ny = length(y);
y_ss = []; y_t = []; y_tF = []; y_tL = [];
for j=1:ny
    eval(['syms ',y{j},'_ss ',y{j},'_t ',y{j},'_tF ',y{j},'_tL'])
    y_ss  = [y_ss, eval([y{j},'_ss'])];
    y_t   = [y_t,  eval([y{j},'_t'])];
    y_tF  = [y_tF, eval([y{j},'_tF'])];
    y_tL  = [y_tL, eval([y{j},'_tL'])];
end

%% Identify which variables are to be log-linearized
% yLogIdx = ~ismember(y(:,1),'b');
% yLogIdx = true(1,size(y,1));
yLogIdx = ~ismember(y,{'lcb'});

%% Exogenous variables
csi = {'hchitil','hXitil','hchitiladd','hXitiladd'};
ncsi = length(csi);
S = rho_csi*eye(ncsi); 
csi_t = []; csi_tF = []; csi_tL = []; eps_t = [];
for j=1:ncsi
    eval(['syms ',csi{j},'_t ',csi{j},'_tF ',csi{j},'_tL'])
    csi_t  = [csi_t,  eval([csi{j},'_t'])];
    csi_tF = [csi_tF, eval([csi{j},'_tF'])];
    csi_tL = [csi_tL, eval([csi{j},'_tL'])];
    eval(['syms eps_',csi{j},'_t'])
    eps_t = [eps_t, eval(['eps_',csi{j},'_t'])];
end
csi_ss = zeros(1,ncsi);

%% Lagrange multipliers for F constraints
nF = 5; for jF=1:nF,FLM{jF}=sprintf('FLM%.0f',jF);end
nF = nF+1; FLM{nF}='Upsilon';
nF = nF+1; FLM{nF}='zeta';
FLM_ss = []; FLM_t = []; FLM_tF = [];
for j=1:nF
    eval(['syms ',FLM{j},'_ss ',FLM{j},'_t ',FLM{j},'_tF'])
    FLM_ss  = [FLM_ss, eval([FLM{j},'_ss'])];
    FLM_t   = [FLM_t,  eval([FLM{j},'_t'])];
    FLM_tF  = [FLM_tF, eval([FLM{j},'_tF'])];
end

%% Lagrange multipliers for G constraints
nG = 4; for jG=1:nG,GLM{jG}=sprintf('GLM%.0f',jG);end
GLM_ss = []; GLM_t = []; GLM_tL = [];
for j=1:nG
    eval(['syms ',GLM{j},'_ss ',GLM{j},'_t ',GLM{j},'_tL'])
    GLM_ss  = [GLM_ss, eval([GLM{j},'_ss'])];
    GLM_t   = [GLM_t,  eval([GLM{j},'_t'])];
    GLM_tL  = [GLM_tL, eval([GLM{j},'_tL'])];
end

%% Extended vector z
nz = ny+ncsi;
z = [y,csi];
z_t = [y_t,csi_t];
z_tF = [y_tF,csi_tF];
z_tL = [y_tL,csi_tL];

%% ------------------------------------------------------------------------

%% Auxiliary definitions, used in later sections

%% parameter values and steady state values
Rd_ss = 1+rd_ss;
Pi_ss = 1;
Delta_ss = 1;
Y_ss = Yss;
b_ss = rho_b*Y_ss;
lcb_ss = 0;
if omega_ss>0
    beta = (delta+1+omega_ss*(delta+(1-delta)*pi_b)-...
        (((delta+1)+omega_ss*(delta+(1-delta)*pi_b))^2-4*delta*(1+omega_ss))^(1/2))...
        /2/delta/(1+omega_ss)/Rd_ss;
elseif omega_ss==0
	beta = 1/Rd_ss;
else
    error('omega_ss needs to be non-negative')
end
Omega_ss = (1-Rd_ss*beta*(delta+(1-delta)*(1-pi_b)))/Rd_ss/beta/(1-delta)/pi_b;
psi_bs = Omega_ss;
psi_s = psi*(pi_b*psi_bs^(-1/nu)+(1-pi_b))^nu;
psi_b = psi_bs*psi_s;
lambda_s_ss = mu_p*(1+omega_y)*mu_w_ss*Hbar_ss^(-nu)*(Y_ss/Z_ss)^(1+omega_y)/...
    ((1-tau_ss)*Y_ss*(pi_b*Omega_ss^(1/nu)*psi_b^(-1/nu)+(1-pi_b)*psi_s^(-1/nu))^nu);
lambda_b_ss = Omega_ss*lambda_s_ss;
Lambda_ss = pi_b*lambda_b_ss+(1-pi_b)*lambda_s_ss;
lambdatil_ss = psi*(pi_b*(lambda_b_ss/psi_b)^(1/nu)+(1-pi_b)*(lambda_s_ss/psi_s)^(1/nu))^nu;
gamma_b = pi_b*(psi/psi_b*lambda_b_ss/lambdatil_ss)^(1/nu);
sb_ss = (1+pi_b*omega_ss-delta*(1+omega_ss)*Rd_ss)*rho_b/pi_b/(1-pi_b)...
    +(gamma_b-pi_b)/pi_b/(1-pi_b)*psi/lambdatil_ss*(1-tau_ss)/mu_p/(1+omega_y);
s_s = s_c-pi_b*sb_ss;
s_b = s_c+(1-pi_b)*sb_ss;
s_bs = s_b/s_s;
s_Xi = omega_ss/eta*rho_b;
s_Xi_cb = 0;
s_g = 1-s_c-s_Xi-s_Xi_cb;
G_ss = s_g*Y_ss;
sigma_s = sigmabar/(pi_b*s_b*sigma_bs+(1-pi_b)*s_s);
sigma_b = sigma_bs*sigma_s;
sigma = sigmabar/s_c;
Cbar_b_ss = s_b*Y_ss*lambda_b_ss^sigma_b;
Cbar_s_ss = s_s*Y_ss*lambda_s_ss^sigma_s;
ctil_b_ss = s_b*Y_ss;
ctil_s_ss = s_s*Y_ss;
omegatil_ss = 1+omega_ss;
K_ss = Lambda_ss*mu_p*(1+omega_y)*psi*mu_w_ss/lambdatil_ss*Hbar_ss^(-nu)*(Y_ss/Z_ss)^(1+omega_y)/(1-alpha*beta);
F_ss = Lambda_ss*(1-tau_ss)*Y_ss/(1-alpha*beta);
chitil_ss = 0;
Xitil_ss = omega_ss/eta/b_ss^(eta-1);

omega_b = (eta-1)*omega_ss/b_ss/(1+omega_ss);
omega_chi = (1+varkappa)/b_ss/(1+omega_ss);
omega_Xi = eta/b_ss/(1+omega_ss);
omega_chiadd = 1/(1+omega_ss);
omega_Xiadd = 1/(1+omega_ss);

if isBind
    Xitil_cb_ss = 0.0086905+BindBp/40000; % in ~isBind it is 0.0086905
else
    Xitil_cb_ss = 1/eta_cb*(eta*Xitil_ss*(b_ss-lcb_ss)^(eta-1)+...
        FLM5_ss/FLM2_ss*eta*(eta-1)*Xitil_ss*(b_ss-lcb_ss)^(eta-2));
end

LMNN_ss = eta*Xitil_ss*(b_ss-lcb_ss)^(eta-1)+...
    FLM5_ss/FLM2_ss*eta*(eta-1)*Xitil_ss*(b_ss-lcb_ss)^(eta-2);

%% a couple more ratios
s_Omega = pi_b*(1-pi_b)*(s_b*sigma_b-s_s*sigma_s)/sigmabar;
s_hb = (psi/psi_b*lambda_b_ss/lambdatil_ss)^(1/nu);
s_hs = (psi/psi_s*lambda_s_ss/lambdatil_ss)^(1/nu);
s_hbs = s_hb/s_hs;

%% shocks
chitil_t = chitil_ss+1/b_ss^(1+varkappa)*hchitil_t;
Xitil_t = Xitil_ss+1/b_ss^eta*hXitil_t;
chitiladd_t = hchitiladd_t;
Xitiladd_t = hXitiladd_t;

%% Monetary Policy rules
ZLBRule.NoZLB = Rd_t-(Rd_ss*Pi_t^phi_pi*(Y_t/Y_ss)^(phi_y/4));
ZLBRule.ZLB = Rd_t-1;
ZLBList = fieldnames(ZLBRule);
nZLB = length(ZLBList);

%% auxiliary definitions
Lambda_t = pi_b*lambda_b_t+(1-pi_b)*lambda_s_t;
lambdatil_t = psi*(pi_b*(lambda_b_t/psi_b)^(1/nu)+(1-pi_b)*(lambda_s_t/psi_s)^(1/nu))^nu;
Lambdatil_t = psi^(1/(1+nu))*(pi_b*psi_b^(-1/nu)*lambda_b_t^((1+nu)/nu)+...
    (1-pi_b)*psi_s^(-1/nu)*lambda_s_t^((1+nu)/nu))^(nu/(1+nu));
ctil_b_t = Cbar_b_ss*lambda_b_t^(-sigma_b);
ctil_s_t = Cbar_s_ss*lambda_s_t^(-sigma_s);
B_t = ctil_b_t-ctil_s_t-((lambda_b_t/psi_b)^(1/nu)-(lambda_s_t/psi_s)^(1/nu))...
    *(lambdatil_t/psi)^(-(1+nu)/nu)*mu_w_ss*Hbar_ss^(-nu)*(Y_t/Z_ss)^(1+omega_y)*Delta_t;

%% Utility function
U = vpa(pi_b*ctil_b_t^(1-1/sigma_b)*Cbar_b_ss^(1/sigma_b)/(1-1/sigma_b)...
    +(1-pi_b)*ctil_s_t^(1-1/sigma_s)*Cbar_s_ss^(1/sigma_s)/(1-1/sigma_s)...
    -psi/(1+nu)*(lambdatil_t/Lambdatil_t)^(-(1+nu)/nu)*Hbar_ss^(-nu)...
    *(Y_t/Z_ss)^(1+omega_y)*Delta_t);

%% G constraints
G = vpa([...
	Rd_t*omegatil_t*beta*...
        ((delta+(1-delta)*pi_b)*lambda_b_tF/Pi_tF+(1-delta)*(1-pi_b)*lambda_s_tF/Pi_tF)-lambda_b_t;
	Rd_t*beta*((1-delta)*pi_b*lambda_b_tF/Pi_tF+...
        (delta+(1-delta)*(1-pi_b))*lambda_s_tF/Pi_tF)-lambda_s_t;
	Lambda_t/lambdatil_t*mu_p*(1+omega_y)*psi*mu_w_ss*Hbar_ss^(-nu)*(Y_t/Z_ss)^(1+omega_y)+...
        alpha*beta*Pi_tF^(theta*(1+omega_y))*K_tF-K_t;
	Lambda_t*(1-tau_ss)*Y_t+alpha*beta*Pi_tF^(theta-1)*F_tF-F_t;
    ]);
nG = length(G);

%% F constraints
for jZLB=1:nZLB,ZLBj = ZLBList{jZLB};
    F.(ZLBj) = vpa([...
        pi_b*(1-pi_b)*B_t+delta*b_tL*omegatil_tL*Rd_tL/Pi_t-(1+pi_b*(omegatil_t-1))*b_t;
        pi_b*Cbar_b_ss*lambda_b_t^(-sigma_b)+(1-pi_b)*Cbar_s_ss*lambda_s_t^(-sigma_s)+G_ss+...
            Xitil_t*(b_t-lcb_t)^eta+Xitiladd_t*(b_t-lcb_t)+...
            Xitil_cb_ss*lcb_t^eta_cb-Y_t;
        alpha*Delta_tL*Pi_t^(theta*(1+omega_y))+...
            (1-alpha)*((1-alpha*Pi_t^(theta-1))/(1-alpha))^(theta*(1+omega_y)/(theta-1))-Delta_t;
        (F_t/K_t)^((theta-1)/(1+omega_y*theta))-(1-alpha*Pi_t^(theta-1))/(1-alpha);
        1+chitil_t*(1+varkappa)*(b_t-lcb_t)^varkappa+chitiladd_t+...
            Xitil_t*eta*(b_t-lcb_t)^(eta-1)+Xitiladd_t-omegatil_t;
        ZLBRule.(ZLBj);
        lcb_t-lcb_ss;
        ]);
end
nF = length(F.NoZLB);

%% ------------------------------------------------------------------------

%% Compute steady state
fprintf('\nSolving for steady state...\n')

%% generate the steady state system
ssSys = vpa(jacobian(U,y_t).'+...
        (FLM_ss*jacobian(F.NoZLB,y_t)).'+beta*(FLM_ss*jacobian(F.NoZLB,y_tL)).'+...
        (GLM_ss*jacobian(G,y_t)).'+beta^(-1)*(GLM_ss*jacobian(G,y_tF)).');
ssSys = [ssSys; G];
ssSys = [ssSys; F.NoZLB];
% plug in the steady state values in symbolic form
ssSys = subs(ssSys,[y_t,y_tF,y_tL,csi_t],[y_ss,y_ss,y_ss,csi_ss]);

%% prepare variables
ssSys1 = eval(ssSys(1:ny));
nSys = length(ssSys1);
x = [FLM,GLM];
nx = length(x);

%% generate function for csolve
SolveFileName = sprintf('ssSys%s',ExerciseName);
fid=fopen([SolveFileName,'.m'],'w');
fprintf(fid,'function f=%s(x) \n',SolveFileName);
fprintf(fid,'f = ones(size(x));\n');
fprintf(fid,'for j=1:size(x,2)  \n');
for j=1:nSys
    fprintf(fid,'%s_ss = x(%.0f,j);  \n',x{j},j);
end
for j=1:nSys
    fprintf(fid,['f(',int2str(j),',j) = ',char(ssSys1(j)),';  \n']);
end
fprintf(fid,'end  \n');
fclose(fid);

%% Solve
[x1,rc] = csolve(SolveFileName,zeros(nx,1),[],NumPrecision,1000);
if rc~=0, error(['Solution of steady state system is not normal, rc = ', int2str(rc)]), end
% [(1:nx)' feval(SolverFileName,x1)] % check system solution
delete([SolveFileName,'.m'])
% x1 = round(x1/NumPrecision)*NumPrecision;

%% evaluate variables
for j=1:nx
    eval(sprintf('%s_ss = %.16f;',x{j},x1(j)))
end

%% Check steady state
ssSys = eval(ssSys);
if ~all(abs(ssSys)<1e-6)
    fprintf('\nWARNING: system solution is not precise\n')
    [(1:length(ssSys))' ssSys]
end

%% fill in a couple of variables
Xitil_cb_ss = eval(Xitil_cb_ss);
LMNN_ss = eta*Xitil_ss*(b_ss-lcb_ss)^(eta-1)+...
    FLM5_ss/FLM2_ss*eta*(eta-1)*Xitil_ss*(b_ss-lcb_ss)^(eta-2);

%% Present steady state
fprintf('\nSteady state results:')
fprintf('\n=====================\n\n')
for j=1:ny
    fprintf('%15s_ss = %12.6f\n',y{j},eval([y{j},'_ss']))
end
for j=1:nF
    fprintf('%15s_ss = %12.6f\n',FLM{j},eval([FLM{j},'_ss']))
end
for j=1:nG
    fprintf('%15s_ss = %12.6f\n',GLM{j},eval([GLM{j},'_ss']))
end
disp(' ')
y_ss = eval(y_ss);
FLM_ss = eval(FLM_ss);
GLM_ss = eval(GLM_ss);

%% present some ratios
fprintf('\nSome ratios:')
fprintf('\n============\n\n')
RatioList = {...
    'omega_ss','delta','pi_b','beta',...
    'rho_b','s_g','s_c','s_b','s_s','s_bs','s_hb','s_hs','s_hbs',...
    'sigmabar','sigma','sigma_b','sigma_s','sigma_bs',...
    'Cbar_b_ss','Cbar_s_ss','ctil_b_ss','ctil_s_ss',...
    'Omega_ss','psi','psi_b','psi_s','psi_bs',...
    'lambdatil_ss','Lambda_ss',...
    'eta','varkappa','chitil_ss','Xitil_ss','s_Xi',...
    'eta_cb','Xitil_cb_ss',...
    'omega_chi','omega_Xi','omega_chiadd','omega_Xiadd','omega_b','elast_omega_b',...
    'zeta_ss','LMNN_ss'...
    };
xdisp(RatioList)
disp(' ')

%% Check that steady state satisfies the constraints imposed in model

% check that delta<beta
fprintf('Checking if delta<beta...')
if delta<beta
    fprintf(' Passed\n')
else
    fprintf(' Failed\n')
    fprintf('Warning: Current parametrization violates restrictions imposed by the model!\n')
end

% Check that lambda_b_ss>lambda_s_ss
fprintf('Checking if lambda_b_ss>lambda_s_ss...')
if lambda_b_ss>lambda_s_ss
    fprintf(' Passed\n')
else
    fprintf(' Failed\n')
    fprintf('Warning: Current parametrization violates restrictions imposed by the model!\n')
end

% Check that one of the following holds:
%   beta*(delta+(1-delta)*pi_b)-delta>=0
% or
%   omega_ss*(delta*(1-beta*(1-pi_b))-beta*pi_b)<=(1-beta)*(beta-delta)
fprintf('Checking if beta*(delta+(1-delta)*pi_b)-delta>=0...')
Check1 = (beta*(delta+(1-delta)*pi_b)-delta>=0);
if Check1, fprintf(' Passed\n'), else fprintf(' Failed\n'), end
fprintf('Checking if omega_ss*(delta*(1-beta*(1-pi_b))-beta*pi_b)<=(1-beta)*(beta-delta)...')
Check2 = (omega_ss*(delta*(1-beta*(1-pi_b))-beta*pi_b)<=(1-beta)*(beta-delta));
if Check2, fprintf(' Passed\n\n'), else fprintf(' Failed\n'), end
if ~Check1&&~Check2
    fprintf('Warning: Current parametrization violates restrictions imposed by the model!\n')
end

%% ------------------------------------------------------------------------

%% Solve for FOC
fprintf('Solving for FOC...\n')

%% FOC
for jZLB=1:nZLB,ZLBj = ZLBList{jZLB};
    FOC.(ZLBj) = vpa(...
        jacobian(U,y_t)+...
        FLM_t*jacobian(F.(ZLBj),y_t)+beta*FLM_tF*subs(jacobian(F.(ZLBj),y_tL),[z_t,z_tL],[z_tF,z_t],0)+...
        GLM_t*jacobian(G,y_t)+beta^(-1)*GLM_tL*subs(jacobian(G,y_tF),[z_t,z_tF],[z_tL,z_t],0)...
        ).';
end

%% ------------------------------------------------------------------------

%% Log Linearize everything
fprintf('Log-linearizing...\n')

%% Substitute variables with logs
hy_t = y_t;
hy_tF = y_tF;
hy_tL = y_tL;
hy_ss = y_ss;
for j=1:ny
    if yLogIdx(j)
        hy_t = subs(hy_t,y_t(j),['h' char(y_t(j))],0);
        hy_tF = subs(hy_tF,y_tF(j),['h' char(y_tF(j))],0);
        hy_tL = subs(hy_tL,y_tL(j),['h' char(y_tL(j))],0);
        hy_ss(j) = log(y_ss(j));
        % Plug in transformed variables into the model
        G = subs(G, [y_t(j),y_tF(j)], exp([hy_t(j),hy_tF(j)]),0);
        for jZLB=1:nZLB,ZLBj = ZLBList{jZLB};
            F.(ZLBj) = subs(F.(ZLBj),[y_t(j),y_tL(j)],exp([hy_t(j),hy_tL(j)]),0);
            FOC.(ZLBj) = subs(FOC.(ZLBj),[y_t(j),y_tL(j),y_tF(j)],exp([hy_t(j),hy_tL(j),hy_tF(j)]),0);
        end
    end
end
hz_t = [hy_t,csi_t];
hz_tF = [hy_tF,csi_tF];
hz_tL = [hy_tL,csi_tL];
hz_ss = [hy_ss,csi_ss];

%% Log-linearize expressions
FnD = {...
    'G_z',G,hz_t;
    'G_zF',G,hz_tF;
    };
for jZLB=1:nZLB,ZLBj = ZLBList{jZLB};
	FnD(end+1,:) = {['F_z.',ZLBj],F.(ZLBj),hz_t};
	FnD(end+1,:) = {['F_zL.',ZLBj],F.(ZLBj),hz_tL};
	FnD(end+1,:) = {['FOC_z.',ZLBj],FOC.(ZLBj),hz_t};
	FnD(end+1,:) = {['FOC_zL.',ZLBj],FOC.(ZLBj),hz_tL};
	FnD(end+1,:) = {['FOC_zF.',ZLBj],FOC.(ZLBj),hz_tF};
	FnD(end+1,:) = {['FOC_FLM.',ZLBj],FOC.(ZLBj),FLM_t};
	FnD(end+1,:) = {['FOC_FLMF.',ZLBj],FOC.(ZLBj),FLM_tF};
	FnD(end+1,:) = {['FOC_GLM.',ZLBj],FOC.(ZLBj),GLM_t};
	FnD(end+1,:) = {['FOC_GLML.',ZLBj],FOC.(ZLBj),GLM_tL};
end
for j=1:size(FnD,1)
    FnDj = vpa(jacobian(FnD{j,2},FnD{j,3}));
    idxSubs = find(FnDj~=0);
    FnDj(idxSubs) = subs(FnDj(idxSubs),...
        [hz_t,hz_tL,hz_tF,FLM_t,FLM_tF,GLM_t,GLM_tL],...
        [hz_ss,hz_ss,hz_ss,FLM_ss,FLM_ss,GLM_ss,GLM_ss],0);
    eval([FnD{j,1},' = eval(FnDj);'])
end
clear FnD

%% Find constants
FnD = {...
    'C_G',G;
    'C_F.NoZLB',F.NoZLB;
    'C_F.ZLB',F.ZLB;
    'C_FOC.NoZLB',FOC.NoZLB;
    'C_FOC.ZLB',FOC.ZLB;
    };
for j=1:size(FnD,1)
    FnDj = FnD{j,2};
    idxSubs = find(FnDj~=0);
    FnDj(idxSubs) = subs(FnDj(idxSubs),...
        [hz_t,hz_tL,hz_tF,FLM_t,FLM_tF,GLM_t,GLM_tL],...
        [hz_ss,hz_ss,hz_ss,FLM_ss,FLM_ss,GLM_ss,GLM_ss],0);
    eval([FnD{j,1},' = eval(FnDj);'])
end
clear FnD

%% ------------------------------------------------------------------------

%% Optimal policy
fprintf('Solving for Optimal policy state space matrices...\n')

%% some variables needed
Lhz_t = hz_t;
Lhz_tL = hz_tL;
LGLM_t = GLM_t;
LGLM_tL = GLM_tL;
for j=1:nz
    Lhz_t = subs(Lhz_t,hz_t(j),['L' char(hz_t(j))],0);
end
for j=1:nG
    LGLM_t = subs(LGLM_t,GLM_t(j),['L' char(GLM_t(j))],0);
end
k_t = [hz_t,FLM_t,GLM_t,Lhz_t,LGLM_t];
nk = length(k_t);

%% Mats for each case
for jZLB=1:nZLB,ZLBj = ZLBList{jZLB};
    G0j = [...
        -FOC_zF.(ZLBj),-FOC_FLMF.(ZLBj),zeros(ny,nz+2*nG);
        -G_zF,zeros(nG,nz+nF+2*nG);
        -F_z.(ZLBj),zeros(nF,nz+nF+2*nG);
        zeros(ncsi,ny),eye(ncsi),zeros(ncsi,nz+nF+2*nG);
        zeros(nz+nG,nz+nF+nG),eye(nz+nG);
        ];
    G1j = [...
        FOC_z.(ZLBj),FOC_FLM.(ZLBj),FOC_GLM.(ZLBj),FOC_zL.(ZLBj),FOC_GLML.(ZLBj);
        G_z,zeros(nG,nz+nF+2*nG);
        F_zL.(ZLBj),zeros(nF,nz+nF+2*nG);
        zeros(ncsi,ny),S,zeros(ncsi,nz+nF+2*nG);
        eye(nz),zeros(nz,nz+nF+2*nG);
        zeros(nG,nz+nF),eye(nG),zeros(nG,nz+nG);
        ];
    Cj = [...
        C_FOC.(ZLBj);
        C_G;
        C_F.(ZLBj);
        zeros(ncsi+nz+nG,1)
        ];
    G2j = [...
        zeros(ny+nG+nF,ncsi);
        eye(ncsi);
        zeros(nz+nG,ncsi);
        ];
    G3j = eye(nk,ny+nG);
    cv = find(all(G0j(1:ny+nG,:)==0,2)==1);
    G0j(cv,:) = -G1j(cv,:);
    G1j(cv,:) = 0;
    G3j(:,cv) = [];
    cv = find(~all(G3j==0,2));
    if ~all(all(G2j(cv,:)==0))
        error('Elements of G2j in forward looking equations are non-zero!')
    end
    Mat.(ZLBj).G0 = G0j;
    Mat.(ZLBj).G1 = G1j;
    Mat.(ZLBj).C = Cj;
    Mat.(ZLBj).G2 = G2j;
    Mat.(ZLBj).G3 = G3j;
    if strcmp(ZLBj,'NoZLB')
        [Phi1j,Cj,Phi2j,fmat,fwt,ywt,gev,eu] = gensys(G0j,G1j,Cj,G2j,G3j);
        if any(eu~=1),fprintf('WARNING: eu = (%.0f,%.0f)\n',eu),end
        REE.NoZLB.Phi1 = Phi1j;
        REE.NoZLB.C = Cj;
        REE.NoZLB.Phi2 = Phi2j;
    end
end

%% ------------------------------------------------------------------------

%% Generate IRFs
fprintf('\nGenerating IRFs...\n')
zz = {z{:},FLM{:},'RdLevel','Rrd','cb','cs','w','wb','ws','Omega',...
    'gammacb','gammacbLevel','lTot','lCB','lPriv','zetalevel','L_LHS'};
nzz = length(zz);
for j=1:nzz
    eval(sprintf('[tf,idx%1$s] = ismember(''%1$s'',zz);',zz{j}))
end
ShockSize = ones(ncsi,1);
ShockSize(ismember(csi,'hchitil')) = dSP/400/omega_chi;
ShockSize(ismember(csi,'hXitil')) = dSP/400/omega_Xi;
ShockSize(ismember(csi,'hchitiladd')) = dSP/400/omega_chiadd;
ShockSize(ismember(csi,'hXitiladd')) = dSP/400/omega_Xiadd;
ShockSize = diag(ShockSize);
if eta_cb==1
    LHScorrection = 1;
else
    LHScorrection = eta_cb*lcb_ss^(eta_cb-1);
end
for jS=1:ncsi,Sj = csi{jS};
    for TZLBj=0:nMax.ZLB
        % Generate REE
        for t=nsteps:-1:TZLBj+1
            REEj(t) = REE.NoZLB;
        end
        for t=TZLBj:-1:1
            G0j = Mat.ZLB.G0;
            G1j = Mat.ZLB.G1;
            Cj  = Mat.ZLB.C;
            G2j = Mat.ZLB.G2;
            G3j = Mat.ZLB.G3;
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
        isFoundSequence = CheckZLBj;
        if isFoundSequence, break, end
    end
    if ~isFoundSequence
        fprintf('WARNING: Did not find sequence for %s!\n',Sj)
    end
    IRF.(Sj) = irfj(1:nz+nF,:);
    TZLB.(Sj) = TZLBj;
    CheckZLB.(Sj) = CheckZLBj;
    clear irfj
end

%% Add more variables
for jS=1:ncsi,Sj = csi{jS};
    irf = IRF.(Sj);
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
    % level, LHS
    irf(idxL_LHS,:) = Xitil_cb_ss+(irf(idxzeta,:)+zeta_ss)./(irf(idxFLM2,:)/100+FLM2_ss);
    % save the irf
    IRF.(Sj) = 100*irf;
    clear irf
end

%% display checks and number of periods
fprintf('\nCheckZLB:\n')
disp(CheckZLB)
fprintf('TZLB:\n')
disp(TZLB)

%% ------------------------------------------------------------------------

%% save information
ExerciseName = ['Output_',ExerciseName];
fprintf('\nMAT file: %s\n',ExerciseName)
save(ExerciseName)

%% Elapsed time
% disp(' '), vctoc(ttic), disp(' ')
% diary off

%% ------------------------------------------------------------------------

