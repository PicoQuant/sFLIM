function Z = smooth1D(Y,lambda)

if size(Y,1)>1
    E  = eye(size(Y,1));
    D1 = diff(E,1);
    D2 = diff(D1,1);
    P  = lambda.^2 .* D2'*D2 + 2.*lambda .* D1'*D1;
    Z  = (E + P) \ Y;
else
    Z = Y;
end;

