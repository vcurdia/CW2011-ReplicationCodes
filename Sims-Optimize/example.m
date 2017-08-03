%example script for csminwel.m, minimizing Rosenbrock function
% Get the function's own help comments
input('Hit return to get the function''s own help comments')
help csminwel
% First with numberical derivatives (grad=[])
input('Hit return to minimize with numerical gradient')
[fhat,xhat,ghat,Hhat,itct,fcount,retcodehat] = csminwel('rsnbrck',[10,-9],eye(2)*.5,[] ,1e-14,100)
% Then with analytic derivatives
input('Hit return to minimize with analytic gradient')
[fhat,xhat,ghat,Hhat,itct,fcount,retcodehat] = csminwel('rsnbrck',[10,-9],eye(2)*.5,'drsnbrck' ,1e-14,100)
