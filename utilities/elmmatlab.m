function  [logelr, lambda, grad, hess, wts, nits] ...
= elmmatlab(z, lam, elmmaxit, gradtol, svdtol, itertrace )

%
%  elmmatlab        computes the empirical likelihood for the mean
%
%  elmmatlab(x, mu, lam, maxit, gradtol, svdtol, itertrace )
%-----------------------------------------------------------------
%  Input      Meaning                            Default
%  
%  z          estimating equation, size is n*p 
%             ================================  
%             hypothesized value of parameters
%             ================================
%  lam        guess for Lagrange multiplier      zeros(1,p)
%  maxit      maximum iterations                 25
%  gradtol    solve grad=0 to this tolerance     1e-7
%  svdtol     tolerance in SVD                   1e-9
%  itertrace  1 to trace iterations              0
%------------------------------------------------------------------
%  Output     Meaning                            Default
%  
%  logelr     log empirical likelihood
%  lambda     Lagrange multiplier
%  grad       gradient vector
%  hess       Hessian matrix
%  wts        output weights
%  nits       number of iterations


%  Check input

if   nargin <1    error('No data provided.');      end

[n p] = size(z);

if  nargin <= 1    lam = zeros(p,1);               end

%if size(lam') ~= size(mu)
%  error('Lam must be sized like mu'' (like a row of x'')');
%end

if  nargin  <= 2      elmmaxit  = 25;     end
if  nargin  <= 3      gradtol   = 1e-7;   end
if  gradtol <  1e-16  gradtol   = 1e-16;  end
if  nargin  <= 4      svdtol    = 1e-9;   end
if  svdtol  <  1e-8;  gradtol   = 1e-8;   end
if  nargin  <= 5      itertrace = 0;      end





% Step weights: inner search starts with
% Newton step, then if necessary gets smaller
% and more parallel to the gradient

newton_wts     =   [ 3.^-(0:3), zeros(1,12)]';
gradient_wts   =   [ 2.^-(0:15) ]';
gradient_wts   =   (gradient_wts.^2-newton_wts.^2).^(0.5);

gradient_wts(12:16,1) = gradient_wts(12:16,1).*[ 10.^-(1:5) ]';


%
%    Outer loop of iteration.  When all goes well
% this is simply a Newton iteration. 
%

nits = 0;
gsize = gradtol + 1.0;
while nits < elmmaxit & gsize > gradtol
  arg  = 1 + z*lam;
  wts1 = plog( arg,1/n,1 );
  wts2 = ( - plog( arg,1/n,2 )).^(0.5);

  grad = -z.*repmat(wts1,1,p);
  grad = sum(grad)';
  gsize = mean(abs(grad));

  hess = z.*repmat(wts2,1,p); % matrix sqrt of hessian

  [hu,hs,hv] = svd(hess,0);
  dhs = diag(hs);
  if  min(dhs) < max(dhs)*svdtol + 1e-128; %guarantee the matrix is not ill
    dhs = dhs + svdtol*max(dhs) + 1e-128;
  end
  dhs = 1./dhs;
  hs = diag(dhs);

  nstep = hv*hs*hu'*(wts1./wts2);
  gstep = -grad;
  if  sum(nstep.^2) < sum(gstep.^2)
    gstep = gstep*sum(nstep.^2)^.5/sum(gstep.^2)^.5;
  end
  
  ologelr = - sum( plog(arg,1/n) );

  ninner = 0;

  for i = 1:length(newton_wts);

    nlam    = lam + newton_wts(i)*nstep + gradient_wts(i)*gstep;
    nlogelr = -sum( plog(  1+z*nlam, 1/n ));

    if  nlogelr < ologelr
      lam = nlam;
      ninner = i;
      break;
    end
  end    

  nits = nits+1;
  if  ninner == 0
     nits = elmmaxit;
  end
  if  itertrace
    disp( [lam,nlogelr,gsize,ninner,nits]);
  end
end


logelr = nlogelr;

if  nargout > 1
  lambda = lam;
end

if  nargout > 2
  grad = grad;
end

if  nargout > 3
  hess = hess'*hess;
end

if nargout > 4
  wts = wts1;
end

if nargout > 5
  nits = nits;
end


%--------------------------------------------------------------------------
function llogz = plog(z,eps,d)
%  Pseudo logarithm
%  The pseudo logarithm agrees with log for arguments larger
%  than eps.  Below eps it is a quadratic.  It has two continuous
%  derivatives.

%    Input
%    
%    z    =  matrix of pseudolog arguments
%    eps  =  threshold
%    d    =  0,1,2, for function, 1st, 2nd deriv
%

if  nargin < 3
  d = 0;
end

zsize = size(z);
out = z(:);
low = out<eps;

if  d == 0
  out(~low) = log( out(~low));
  out( low) = log(eps) -1.5 +2*out(low)/eps -0.5*(out(low)/eps).^2;
elseif d == 1
  out(~low) = 1./out(~low);
  out( low) = 2/eps -out(low)/eps^2;
elseif d == 2
  out(~low) = -1./out(~low).^2;
  out( low) =  -1/eps^2;
else
  error('Unknown option d for pseudologarithm');
end

llogz = reshape(out,zsize); % guarantee that the size of the reaults is the
                            % same as the input 'z'.
