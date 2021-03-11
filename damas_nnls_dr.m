function x = damas_nnls_dr(D, Bf, epsilon)

%% fast NNLS for CMF-NNLS, diagonal removal

n = size(D, 2);

tic;

R = true(n, 1);
N = 1:n;

x = zeros(n, 1);


%% w = At (y - Ax)
Ay = proddamasDR(D, Bf);
w = Ay - proddamasDR(D, proddamasDR(D, x));

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
    DR = D(:, ~R);
    aDR2 = (abs(D).^2);
    aDR2r = (abs(DR).^2);
    
    AR = abs(D'*DR).^2 - aDR2' * aDR2r;% (abs(DR').^2) *  (abs(DR).^2);
   
    s(~R) = AR\ Bf;
       
    
    
    while min(s(~R)) <= 0
        Q = (s <= 0) & (~R);
        alpha = min(x(Q)./(x(Q)-s(Q)));
        
        x = x + alpha * (s - x);
        R = ((x <= eps) & ~R) | R;

        s(:) = 0;
        
        %% small LS problem
       DR = D(:, ~R);
    aDR2 = (abs(D).^2);
    aDR2r = (abs(DR).^2);

    AR = abs(D'*DR).^2 - aDR2' * aDR2r;% (abs(DR').^2) *  (abs(DR).^2);
   
    s(~R) = AR\ Bf;     
    end
    x = s;
    %% w = At (y - Ax)

    w = Ay - proddamasDR(D, proddamasDR_support(D, x, ~R));
    iter = iter + 1;


end
end
    