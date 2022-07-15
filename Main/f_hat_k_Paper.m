function [obj]=f_hat_k_Paper(y,Paras,W1,P)
    %Author: Feng-Yi Liao
    
    G = Paras.c_sdp+Paras.At_sdp'*y;
    cnew = [G'*W1];
    cnew = [cnew;reshape(P'*mat(G)*P,[],1 )];
    Knew.l = 1;
    %Knew.s = Paras.n_new*ones(1,Paras.NumOfP);
    Knew.s = Paras.MaxCols;
    
    Anew = zeros(1,Knew.l+Paras.MaxCols^2);
    %add trace constraint
    Anew(1) = 1;
    Anew((Knew.l+1):end) = reshape(eye(Paras.MaxCols),[],1);
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
