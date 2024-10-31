% #' Function for converting densities to log quantile density functions 
% #' 
% #' @param dens density values on dSup - must be strictly positive and integrate to 1
% #' @param dSup support (grid) for Density domain
% #' @param lqdSup support for lqd domain - must begin at 0 and end at 1; default [0,1] with N-equidistant support points
% #' @param N desired number of points on a [0,1] grid for lqd function; default length(dSup)
% #' 
% #' @return lqd log quantile density on lqdSup
% #' 

function [lqd] = dens2lqd(dens, dSup, N, lqdSup)
    if nargin == 3
        lqdSup = linspace(0,1,N);
    end
    if any(dens <= 0)
        error('Please correct negative or zero probability density estimates.')
    end
    % Check density requirements
    if abs(trapz(dSup, dens) - 1) > 1e-5
        warning('Density does not integrate to 1 with tolerance of 1e-5 - renormalizing now.')
        dens = dens / trapz(dSup, dens);
    end

    % Get CDF
    qtemp = cumtrapz(dSup, dens);
    lqd_temp = -log(dens);
    lqd = spline(qtemp, lqd_temp, lqdSup);
end 
