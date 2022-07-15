function [Wstar,y_next,U_next,xstar] = Direction(xk,Paras,W1,W2,U)   
    %Author: Feng-Yi Liao    
    %W1 and W2 are fixed atoms
    %U is the transformation matrix for basis pursuit
    
    G = Paras.At_sdp'*xk+Paras.c_sdp;
    cnew = [G'*W1;G'*W2;0;1/Paras.u/2;-1/Paras.u*Paras.b_sdp];
    Anew = [];
    Al = [Paras.At_sdp*W1,Paras.At_sdp*W2];
    AFW2 = zeros(Paras.m,Paras.NumOfVar_new*Paras.NumOfP);%factor width part
    IndicesAll = Paras.IndicesAll;
    Indices = Paras.Indices;
    G_update = reshape(U*mat(G)*U',[],1);
    At_sdp_update = zeros(Paras.m,Paras.n^2);
    for i = 1:Paras.m
        At_sdp_update(i,:) = reshape(U*mat(Paras.At_sdp(i,:))*U',[],1);
    end
    for i = 1:Paras.NumOfP
        idx = IndicesAll(i,:);
        cnew = [cnew;G_update(idx)];
        AFW2(:,Paras.NumOfVar_new*(i-1)+1:Paras.NumOfVar_new*i) = At_sdp_update(:,idx);
    end
    
    Aq = [zeros(Paras.m,2),-eye(Paras.m)];
    Anew = [Al,Aq,AFW2];
    bnew = zeros(Paras.m,1);
    Knew.s = Paras.n_new*ones(1,Paras.NumOfP);
    Knew.r = [Paras.m+2];
    Knew.l = 2;
    
    %add trace constraint for feasible set
    [hei,wei] = size(Anew);
    Atrace = zeros(1,wei);
    Atrace(1) = 1; Atrace(2) = 1;
    Atrace((Knew.l+Knew.r+1):end) = repmat(reshape(eye(Paras.n_new),1,[]),1,Paras.NumOfP);
    
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
    W = zeros(Paras.n);
    start = 1;
    for i = 1:Paras.NumOfP
        idx = Indices(i,:);
        W(idx) = W(idx) + xpsd(start:(start+Paras.NumOfVar_new_sym-1))';
        start = start + Paras.NumOfVar_new_sym;
    end
    W = W+tril(W,-1)';
    W = reshape(U'*W*U,[],1);
    Wstar = xl(1)*W1+xl(2)*W2+W;
    y_next = xk + (Paras.At_sdp*Wstar-Paras.b_sdp)/Paras.u;
    [eigvec,eigval] = eig(mat(Wstar));
    U_next = eigvec';
    xstar = xl(1:2);

end
