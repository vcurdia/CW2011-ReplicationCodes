function  [F,eno]=lyapcsd(R,S)
%function  [F,eno]=lyapcsd(R,S)
% solves RFR'+S=F
% Will fail if R has a root that is one in absolute value.
IR=eye(size(R))+R;
[F,eno]=lyapcs(-inv(IR),R'/IR',(IR\S)/IR');
