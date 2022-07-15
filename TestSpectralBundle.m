%Test SpectralBundle
%clc;clear; 
dx = 1;
n = 200;m=101;
load('n100m101.mat');
At_sdp = full(At_sdp); b_sdp = full(b_sdp); c_sdp = full(c_sdp);

% [OBJ_Inner,X] = InnerApproximation(At_sdp,b_sdp,c_sdp,K_sdp,1,10)

[x_sdp,~,~]=sedumi(At_sdp,b_sdp,c_sdp,K_sdp);
TrueObj = c_sdp'*x_sdp;
%TrueObj = 16.904135906958306 ;

opts.Maxiter =100; opts.dx = dx;  %size of a block
opts.n = n; opts.m = m; opts.epislon = 10^-3;
opts.ml = 0.2; opts.u = 0.25;
opts.MaxCols = 25; %for paper version
tic;
    Obj = SpectralBundle(At_sdp,b_sdp,c_sdp,K_sdp,opts);
    %Obj = SpectralBundle_Paper(At_sdp,b_sdp,c_sdp,K_sdp,opts);
toc;
%plot(1:length(Obj),TrueObj*ones(1,length(Obj)),'r-',1:length(Obj),-Obj,'b-o');
plot(1:length(Obj),TrueObj+Obj,'b-o');
xlabel('iteration'); ylabel('error');