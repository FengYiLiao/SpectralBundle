function [Wstar,y_next,xstar,Vstar] = Direction_Paper(xk,Paras,W1,P)   
    %Author: Feng-Yi Liao
    %W1 is a fixed atoms
    %P is the transformation matrix
    
    G = Paras.At_sdp'*xk+Paras.c_sdp;
    cnew = [G'*W1;0;1/Paras.u/2;-1/Paras.u*Paras.b_sdp;reshape(P'*mat(G)*P,[],1)];
    Anew = [];
    Al = [Paras.At_sdp*W1];
    APSD = zeros(Paras.m,Paras.MaxCols^2); %PSD part
    for i = 1:Paras.m
          APSD(i,:) = reshape(P'*mat(Paras.At_sdp(i,:))*P,1,[]);
    end
    Aq = [zeros(Paras.m,2),-eye(Paras.m)];
    Anew = [Al,Aq,APSD];
    bnew = zeros(Paras.m,1);
    Knew.s = Paras.MaxCols;
    Knew.r = [Paras.m+2];
    Knew.l = 1;
    
    %add trace constraint for feasible set
    [hei,wei] = size(Anew);
    Atrace = zeros(1,wei);
    Atrace(1) = 1; 
    Atrace((Knew.l+Knew.r+1):end) = reshape(eye(Paras.MaxCols),[],1);
    
    Anew = [Anew;Atrace];
    bnew = [bnew;1];
    
    %add constraint on rotated cone %force the first element to be 1/2
    Ar = zeros(1,wei);
    Ar(Knew.l+1) = 1;
    
    Anew = [Anew;Ar];
    bnew = [bnew;1/2];
    
    prob = SedumiToMosek(Anew,bnew,cnew,Knew);
    [~, res1] = mosekopt('minimize info', prob);
    status = res1.sol.itr.prosta;
    if ~strcmp(status,'PRIMAL_AND_DUAL_FEASIBLE')
       warning('Mosek:  PRIMAL_AND_DUAL_FEASIBLE') ;
    end
    xpsd = res1.sol.itr.barx;% psd
    xl = res1.sol.itr.xx; %non negative and rotated
    %[x,y,info] = sedumi(Anew,bnew,cnew,Knew);
    
    %recover W
    
    V = zeros(Paras.MaxCols);
    V(Paras.IndicesPSD) = xpsd;
    V = V + tril(V,-1)';
    Vstar = V ;%for return

    W = P*V*P'; 
    
    Wstar = xl(1)*W1 + reshape(W,[],1);
    y_next = xk + (Paras.At_sdp*Wstar-Paras.b_sdp)/Paras.u;
%     [eigvec,eigval] = eig(mat(Wstar));
%     U_next = eigvec';
    xstar = xl(1);

end