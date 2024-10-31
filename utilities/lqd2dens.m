% #' Function for converting log quantile densities to densities
% #' 
% #' @param  lqd log quantile density on lqdSup
% #' @param  lqdSup support forlqd domain - must begin at 0 and end at 1
% #' @param  dSup support for Density domain - max and min values mark the boundary of the support.
% #' @param  useSplines fit spline to the lqd when doing the numerical integration (default: TRUE)
% #' 
% #' @return dens density values on dSup

function [dens] = lqd2dens(lqd, lqdSup, dSup, useSplines)
    if nargin == 3 
        useSplines = true;
    end

    if ~isequal([min(lqdSup), max(lqdSup)], [0,1] )
        error("Please check the support of the LQD domain's boundaries.")
    end

    if useSplines
        % Could fit spline if this yields more accurate numerical integration
        lqd_sp = @(t) spline(lqdSup, lqd, t); 
        lqd_exp = @(t) exp(lqd_sp(t));
        % Get grid and function for density space
        integrate = @(i) integral(lqd_exp, lqdSup(i-1), lqdSup(i));
        dtemp_cumsum = zeros(1, length(lqdSup));
        for i = 2:length(lqdSup)
            dtemp_cumsum(i) = dtemp_cumsum(i-1) + integrate(i); 
        end
        dtemp = dSup(1) + dtemp_cumsum; 
    else
        % Get grid and function for density space
        dtemp = dSup(1) + cumtrapz(lqdSup, exp(lqd));
    end

    dens_temp = exp(-lqd);

    % Fix the support of the density to match dSup
    % r1 = range(dtemp);
    % r2 = range(dSup);
    r1 = max(dtemp) - min(dtemp); 
    r2 = max(dSup) - min(dSup);
    dtemp = (dtemp - dtemp(1)) * r2 / r1 + dtemp(1); 

    % Remove duplicates
    [dtemp, ia, ~] = unique(dtemp);

    % Interpolate to dSup and normalize
    dens = interp1(dtemp, dens_temp(ia), dSup, "linear", "extrap");
    dens = dens / trapz(dSup, dens); % Normalize 
end