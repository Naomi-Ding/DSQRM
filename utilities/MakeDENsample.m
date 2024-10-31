% #' Convenience function for converting log quantile densities to densities                                                                                                     
% #' 
% #' See 'lqd2dens' and 'DeregulariseByAlpha' for more details.
% #' This function transforms the log quantile densities in 'qmatrix' to density functions, optionally followed by deregularisation.
% #' 
% #' @param qmatrix Matrix holding the log quantile density values on [0,1]
% #' @param lqdSup Support grid for input log quantile densities (default = seq(0, 1, length.out = ncol(qmatrix)))
% #' @param dSup Support grid for output densities (default = seq(0, 1, length.out = ncol(qmatrix)))
% #' @param useAlpha Logical indicator to deregularise the densities (default = FALSE)
% #' @param alpha Scalar to deregularise the density - where possible, this will be the minimum value for the deregularised densities (default=0)

% #' @return list with the 'DEN' transformed data, and 'dSup' that matches the input argument.

function [DEN, dSup] = MakeDENsample(qmatrix, lqdSup, dSup)
    [n, m] = size(qmatrix);
    if nargin < 3
        if nargin < 2
            lqdSup = linspace(0,1, m);
        end
        dSup = linspace(0,1, m);
    end

    % A. Get the densities from the log-quantile-density projections
    qfun = @(u) lqd2dens(u, lqdSup, dSup, true);
    N = length(dSup); 
    DEN = zeros(n, N);
    for i = 1:n
        DEN(i, :) = qfun(qmatrix(i, :));
    end

end 
