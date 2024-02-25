function [beta B h] = Mimd1(X,y,nd)

[n,p] = size(X);
X = X - repmat(mean(X),n,1);
S = inv(X'*X/n+eye(p)/n/n)^0.5;
X = X*S;

X1 = [X ones(n,1)];
beta = inv(X1'*X1 + eye(p+1)/n/n )*(X1'*y);
%y = y-X1*beta;

B = eye(p);
m = p;
for iter = 1:p+5;
	ab = ones(p,n);
    if iter<3
        h = n^(-1/(m+4));
    end
    if iter/3 == floor(iter/3)
        h = n^(-1/(m+4));
%        h = cvm(X*B, y);
    end
	for i = 1:n;
       	xi = X - repmat(X(i,:),n,1);
        kernel = exp(-sum((xi*B).^2,2)/(2*h*h));
        onexi = [xi ones(n,1)];
        xk = onexi.*repmat(kernel, 1, p+1);
        abi = inv(xk'*onexi+eye(p+1)*0.0001)*xk'*y;
        ab(:,i) = abi(1:p);
	end;
    ab = ab - repmat(mean(ab,2), 1, n);
	[B0 D] = eig(ab*ab');
	[D I] = sort(diag(D));
	B = B0(:,I);
    B = B(:, p+1-(1:p));
    m = max(nd, m-1);
    B = B(:,1:m);
end;

beta = S*(eye(p)-B*B')*beta(1:p);

B = S*B;
for i = 1:size(B,2);
    B(:,i) = B(:,i)/sqrt(B(:,i)'*B(:,i));
end