function Obj = SpectralBundle_Paper(A_sdp,b_sdp,c_sdp,K_sdp,opts)
    %Author: Feng-Yi Liao
    %Spectral Bundle method (the same algorithm as the original paper)
    Paras = Initialize(A_sdp,b_sdp,c_sdp,opts);
    %W1 = reshape(eye(opts.n),[],1);
    x = ones(opts.m,1); %initial point
    [eig_vec,eig_val] = eig(mat(-c_sdp-A_sdp'*x));
    [~,I] = sort(diag(eig_val),'descend');
    eig_vec = eig_vec(:,I);
    P = eig_vec(:,1:opts.MaxCols);  %transformation matrix
    W1 = reshape(eig_vec(:,1)*eig_vec(:,1)',[],1);
    %W2 = reshape(V*V',[],1);
    Obj = [];
    for iter = 1:Paras.Maxiter

        [~,y_next,xstar,Vstar] = Direction_Paper(x,Paras,W1,P);
        [eig_vec,eig_val] = eig(Vstar);
        [~,I] = sort(diag(eig_val),'descend');
        eig_vec = eig_vec(:,I);
        eig_val = eig_val(:,I);
        Q1 = eig_vec(:,1:Paras.MaxCols-1);
        Q2 = eig_vec(:,Paras.MaxCols:end);
        %Sigma1 = eig_val(1:Paras.MaxCols-1,1:Paras.MaxCols-1);
        Sigma2 = eig_val(Paras.MaxCols:end,Paras.MaxCols:end);
        
        %Stopping criteria
        f1 = max(eig(mat(-Paras.c_sdp-Paras.At_sdp'*x)))+Paras.b_sdp'*x;
        f2 = f_hat_k_Paper(y_next,Paras,W1,P);
        if f1-f2<0
            warning('something wrong');
        elseif f1-f2 < Paras.epislon
            Obj = [Obj,f2];
            break;
        end
        
        %Serious Step or Null Step
        f3 = max(eig(mat(-Paras.c_sdp-Paras.At_sdp'*y_next)))+Paras.b_sdp'*y_next;
        Threshold = f1-Paras.ml*(f1-f2);
        if f3 <= Threshold
            %serious step
            x = y_next;
        end
        
        %Update W1 and P
        W1 = (xstar(1)*W1+ reshape(P*Q2*Sigma2*Q2'*P',[],1))/(xstar(1)+trace(Sigma2));
        [eig_vec,eig_val] = eig(mat(-c_sdp-A_sdp'*y_next));
        [~,I] = sort(diag(eig_val),'descend');
        eig_vec = eig_vec(:,I);
        v = eig_vec(:,1);    
        [Q,~] = gson([P*Q1,v]); %Orthogonalization the function is written from another group
        P = Q;
        Obj = [Obj,f2];
    end
end

function Paras = Initialize(At_sdp,b_sdp,c_sdp,opts)
    Paras.epislon = opts.epislon;
    Paras.ml = opts.ml;
    Paras.u = opts.u;
    IndicesAll = BIGPSDpositionAll(opts.n,opts.dx);%The nonzero indices in a n x n matrix (including lower and upper)
    Indices = BIGPSDposition(opts.n,opts.dx); %The nonzero indices in a n x n matrix (only symmetric part)
    Paras.IndicesAll = IndicesAll;
    Paras.Indices = Indices;
    Paras.Maxiter = opts.Maxiter;
    Paras.dx = opts.dx;
    Paras.n = opts.n;
    Paras.m = opts.m;
    Paras.n_new = 2*Paras.dx;
    Paras.NumOfVar_new = Paras.n_new^2;
    Paras.NumOfVar_new_sym = Paras.n_new*(Paras.n_new+1)/2;
    Paras.NumOfP = nchoosek(Paras.n/Paras.dx,2);%Number of Blocks
    Paras.b_sdp = b_sdp;
    Paras.At_sdp = At_sdp;
    Paras.c_sdp = c_sdp;
    Paras.MaxCols = opts.MaxCols;
    [xIndSym,~,~,~,~,~] = SymmetricIndices(Paras.MaxCols);
    Paras.IndicesPSD = xIndSym;
end



