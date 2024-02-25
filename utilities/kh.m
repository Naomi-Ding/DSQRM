%% Kernel Function
function kernel = kh(w, h)

tmp = w/h;
m = length(w);
kernel = zeros(m,1);

for i = 1:m
    if (abs(tmp(i)) <= 1)
        kernel(i) = (1 - tmp(i)^2) / h * 0.75;
    end
end

end
