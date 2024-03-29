function [x_out,llh_out,n_iter,conv,b]  = Convergence(I,K)
tic;
i_err = (min(K,[],2)>0);
K = K(i_err,:);
I = I(:,i_err);

lI = log(I)-1;

n  = size(I,1);
k  = size(K,2);

hd = repmat(sum(K,1),[n 1]);

X = repmat(sum(I,2),[1 k])./k/5;
Exitflag = false;

me = zeros(1,20);
conv=[];
n_iter=0;
update = 1.2;
while ~Exitflag
    
    rate = update;
    
    for sl = 1:20
        
        Z  = X*K';
        
        
        x  = (X./hd).*((I./Z)*K);
        dx = (x - X);
        
        X  = X + rate.*dx;
        
        me(sl) = 25*max(sum(abs(dx),2));
    end
    
    if min(me) < 50
        update = 1.6;
    end
    
    if min(me) < 1
        Exitflag = true;
    else
        Exitflag = false;
    end
    
    n_iter=n_iter+20;
    conv=[conv me];  %#ok<AGROW>
    
end

m   = Z.*(log(Z)-lI);
llh = sum(m,2);
x_out   = x;
llh_out = llh;
b=toc;
end