function [x_out,llh_out,n_iter,conv,b]  = Momentum_convergence_2(I,K)
tic;
i_err = (min(K,[],2)>0);
K = K(i_err,:);
I = I(:,i_err);


lI = log(I)-1;



n  = size(I,1);
k  = size(K,2);

hd = repmat(sum(K,1),[n 1]);

% X = K\I';
% X = mKx(cKt(3,zeros(n,k), X'),[],3);

X = repmat(sum(I,2),[1 k])./k/5;
Exitflag = false;

update_rate = 1.1; % 110 percent of the difference from previous update
mu = 0.1;
residual  =  [];
v = zeros(size(X));



me = zeros(1,20);
conv=[];
n_iter=0;


while ~Exitflag
    
    
    
    Z  = X*K';
    x  = (X./hd).*((I./Z)*K);
    dx = (x - X);  % Gradient at each update step
    
    %% Momentum update
    v_prev   = v; % back this up
    v        = mu * v_prev + update_rate * dx; % velocity update stays the same
    X        = X +  (update_rate - 1) * v; % position update changes form
    %         residual = [residual 10*max(sum(abs(dx),2))]; %#ok<AGROW>
    residual = 100*max(sum(abs(dx),2));
    
    
    n_iter  = n_iter  + 1;
    
    if mod(n_iter,10) == 0 && update_rate < 1.8
        
        update_rate  = update_rate + 0.1;
        mu           = mu + 0.1;
    end
    
    if residual < 1
        Exitflag = true;
    else
        Exitflag = false;
    end
    
end

m   = Z.*(log(Z)-lI);
llh = sum(m,2);
x_out   = x;
llh_out = llh;
b=toc;
end