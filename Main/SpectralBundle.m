function Obj = SpectralBundle(A_sdp,b_sdp,c_sdp,K_sdp,opts)
    %Author: Feng-Yi Liao
    %Spectral Bundle method with Basis Pursuit using factor width matrices
    Paras = Initialize(A_sdp,b_sdp,c_sdp,opts);

    W1 = reshape(eye(opts.n),[],1);
    W2 = reshape(eye(opts.n),[],1);
    U = eye(opts.n);
    x = ones(opts.m,1); %initial point
    [eig_vec,eig_val] = eig(mat(-c_sdp-A_sdp'*x));
    [~,I] = sort(diag(eig_val),'descend');
    eig_vec = eig_vec(:,I);
    V = eig_vec(:,1);    
    W2 = reshape(V*V',[],1);
    Obj = [];
    for iter = 1:Paras.Maxiter
        [W,y_next,U_next,xstar] = Direction(x,Paras,W1,W2,U);
        
        %Stopping criteria
        f1 = max(eig(mat(-Paras.c_sdp-Paras.At_sdp'*x)))+Paras.b_sdp'*x;
        f2 = f_hat_k(y_next,Paras,W1,W2,U);
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
        
        %Update W1, W2, and U
        U = U_next;
        %it's important to update W1 first
        W1 = (xstar(1)*W1 + xstar(2)*W2)/(xstar(1)+xstar(2));
        [eig_vec,eig_val] = eig(mat(-c_sdp-A_sdp'*y_next));
        [~,I] = sort(diag(eig_val),'descend');
        eig_vec = eig_vec(:,I);
        V = eig_vec(:,1);    
        W2 = reshape(V*V',[],1); 
        
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
end



