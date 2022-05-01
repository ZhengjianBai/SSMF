function X = SimplexProj(Y)

[N,D] = size(Y);
X = sort(Y,2,'descend');
if abs(sum(X))>1e-36
Xtmp = (cumsum(X,2)-1)*diag(sparse(1./(1:D)));
X = max(bsxfun(@minus,Y,Xtmp(sub2ind([N,D],(1:N)',sum(X>Xtmp,2)))),0);
else
    X=ones(1,D)/D;
end