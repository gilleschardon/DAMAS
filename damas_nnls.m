function [x, T, Titer] = damas_nnls(D, Bf, epsilon)

%% fast NNLS for CMF-NNLS, diagonal removal

n = size(D, 2);

tic;

R = true(n, 1);
N = 1:n;

x = zeros(n, 1);

tic
%% w = At (y - Ax)
Ay = proddamas(D, Bf);
w = Ay - proddamas(D, proddamas(D, x));

iter = 1;
%%
while (any(R)) && (max(w(R)) > epsilon)
    
    % logs
    %obj = norm(D(:, ~R)*spdiags(x(~R),0,sum(~R), sum(~R))*D(:, ~R)' - Data, 'fro');
    fprintf("T %.2fs Iter %u    Maxw %.4e\n", toc, iter, max(w(R)));

    [~, idx] = max(w(R));
    Ridx = N(R);
    idx = Ridx(idx);
    R(idx) = false;
    s = zeros(size(x));
    
    %% small LS problem
    AR = abs(D' * D(:, ~R)).^2;

    s(~R) = AR\ Bf;
       
    
    
    while min(s(~R)) <= 0
        Q = (s <= 0) & (~R);
        alpha = min(x(Q)./(x(Q)-s(Q)));
        
        x = x + alpha * (s - x);
        R = ((x <= eps) & ~R) | R;

        s(:) = 0;
        
        %% small LS problem
    AR = abs(D' * D(:, ~R)).^2;
    
       s(~R) = AR\ Bf;       
    end
    x = s;
    %% w = At (y - Ax)



    w = Ay - proddamas(D, proddamas_support(D, x, ~R));
    iter = iter + 1;

end

T = toc;
Titer = T/iter;
end
    