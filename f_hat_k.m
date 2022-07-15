function [obj]=f_hat_k(y,Paras,W1,W2,U)
    G = Paras.c_sdp+Paras.At_sdp'*y;
    cnew = [G'*W1;G'*W2];
    IndicesAll = Paras.IndicesAll;
    G_update = reshape(U*mat(G)*U',[],1);
    for i = 1:Paras.NumOfP
        idx = IndicesAll(i,:);
        cnew = [cnew;G_update(idx)];
    end
    Knew.l = 2;
    Knew.s = Paras.n_new*ones(1,Paras.NumOfP);
    Anew = zeros(1,Knew.l+Paras.NumOfP*Paras.NumOfVar_new);
    
    %add trace constraint
    Anew(1) = 1; Anew(2) = 1;
    Anew((Knew.l+1):end) = repmat(reshape(eye(Paras.n_new),1,[]),1,Paras.NumOfP);
    bnew = [1];
    
    prob = SedumiToMosek(Anew,bnew,cnew,Knew);
    [~, res1] = mosekopt('minimize info', prob);
    status = res1.sol.itr.prosta;
    if ~strcmp(status,'PRIMAL_AND_DUAL_FEASIBLE')
       warning('Mosek:  PRIMAL_AND_DUAL_FEASIBLE') ;
    end
    
    %[x,~,info]=sedumi(Anew,bnew,cnew,Knew);
    %obj = -cnew'*x+Paras.b_sdp'*y;
    obj = -res1.sol.itr.pobjval+Paras.b_sdp'*y;
end