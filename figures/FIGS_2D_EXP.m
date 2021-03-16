%% Compares DAMAS, quadprog, lsqnonneg and Lawson-Hanson on experimental data
%
% 2D grid of 10 800 points
% plots the objective function in function of time
% and power maps for
% cyclic and random DAMAS
% quadprog
% Lanson-Hawnson, with and without factorization

% Approx 1h



% loads the data (cov. matrix, frequency, array coordinates)

load damasdata2D
addpath('..')

% we select the 64 inner microphones
N = Pmic(:, 2).^2 + Pmic(:, 1).^2;
[~, order] = sort(N);
Nm = 64;
Z = order(1:Nm);
% positions of the microphones
Pmic = Pmic(Z, :);
% spatial covariance matrix
Data = Data(Z, Z);

close all

% source grid
Lx = 180;
Ly = 60;
Lz = 1;
xx = linspace(-2, 1, Lx)';
yy = linspace(-1, 0, Ly)';
zz = 4.4;

[Xg, Yg, Zg] = meshgrid(xx, yy, zz);

% dictionary of sources
D = dictionary(Pmic, [Xg(:) Yg(:) Zg(:)], k);

%% Beamforming (no normalization)
% used as input for DAMAS
Cbf = sum(conj(D) .* (Data*D), 1);
Cbf = real(Cbf');

%% Beamforming (normalized)
% used to plot the power map
Dbf = D ./ sum(abs(D).^2, 1);
Cbfn = sum(conj(Dbf) .* (Data*Dbf), 1);
Cbfn = real(Cbfn');

%% Large matrix for DAMAS (A in the paper)
DD = abs(D'*D).^2;

%% Quadprog, interior point
% an iteration is a new barrier

niters = 5:5:25;
for uu = 1:length(niters)
options = optimoptions('quadprog', 'MaxIterations', niters(uu), 'OptimalityTolerance', 1e-12);

tic
[xquad,fval,exitflag,output] = quadprog(DD, -real(Cbf), [], [], [], [], zeros(Lx*Ly, 1), [], [], options);
% Time
Tquad(uu) = toc;
% Objective
Cquad(uu) = norm(D*diag(xquad)*D' - Data, 'fro')^2;
end


%% cyclic DAMAS
Niter = [5 10 20 50 100 200 500]*1e5;
Tdamas = zeros(size(Niter));
Cdamas = zeros(size(Niter));

xdamas = zeros(size(DD, 1));

    tic;
    xdamas = damas(DD, real(Cbf), Niter(1), xdamas);
    Tdamas(1) = toc;
    Cdamas(1) = norm(D*diag(xdamas)*D' - Data, 'fro')^2;
    
    % use the last output as a starting point
for u = 2:length(Niter)
    tic;
    xdamas = damas(DD, real(Cbf), Niter(u)- Niter(u-1), xdamas);
    Tdamas(u) = toc + Tdamas(u-1);
    Cdamas(u) = norm(D*diag(xdamas)*D' - Data, 'fro')^2;
end


%% random DAMAS
% this one is shown to attain the minimum, and to converge towards the
% argmin if it is unique
Niter = [5 10 20 50 100 200 500]*1e5;
Tdamas_rand = zeros(size(Niter));
Cdamas_rand = zeros(size(Niter));

xdamas_rand = zeros(size(DD, 1));

    tic;
    xdamas_rand = damas_rand(DD, real(Cbf), Niter(1), xdamas_rand);
    Tdamas_rand(1) = toc;
    Cdamas_rand(1) = norm(D*diag(xdamas_rand)*D' - Data, 'fro')^2;
    
for u = 2:length(Niter)
    tic;
    xdamas_rand = damas_rand(DD, real(Cbf), Niter(u)- Niter(u-1), xdamas_rand);
    Tdamas_rand(u) = toc + Tdamas_rand(u-1);
    Cdamas_rand(u) = norm(D*diag(xdamas_rand)*D' - Data, 'fro')^2;
end

%% lsqnonneg
% using the MATLAB function
% (convergence is not plotted, for legal reasons, see results in the paper)

M = size(D, 2);
N = size(D, 1);

% We prepare the matrix Dtilde (lsqnonneg does not take complex matrices)
DDDri = zeros(2*N^2, M);

for u = 1:M
    d = D(:, u)*D(:, u)';

    DDDri(:, u) = [real(d(:)) ; imag(d(:))];
end

Xri = [real(Data(:)) ; imag(Data(:))];

tic
xlsqnonneg = lsqnonneg(DDDri, Xri);
Tlsqnonneg = toc;
Clsqnonneg = norm(D*diag(xlsqnonneg)*D' - Data, 'fro')^2;

%% optimized LH
% the fastest, using the factorization of A

[xlh, Clh, Tlh, w] = cmf_nnls(D, Data, 0);

%% DAMAS-NNLS
% /!\ this does not solve the same optimization problem

tic
[xnnls] = damas_nnls(D, Cbf, 0);
TDAMASNNLS = toc;

%%
clear DD DDDri
save 2Dexp


%% Plots

MMM = Clh(end);

% Objective function in function of computational time
figure
set(gcf, 'Position',  [100, 100, 500, 400])

%scatter(Tlh, Clh, 100, 'black', 'x', 'linewidth', 2)
loglog(Tlh, (Clh-MMM)/MMM, 'k--', 'linewidth', 1);
hold on
%plscatter(Tlh, Clh, 100, 'black', 'x', 'linewidth', 2)
plot(Tdamas, (Cdamas-MMM)/MMM, 'k-x', 'linewidth', 2, 'Markersize', 15)
plot(Tdamas_rand(1:end-1), (Cdamas_rand(1:end-1)-MMM)/MMM, 'k--o', 'linewidth', 2, 'Markersize', 10)

plot(Tquad, (Cquad-MMM)/MMM, 'k-.+', 'linewidth', 2, 'Markersize', 15)

step = 25;

scatter(Tlsqnonneg(step:step:end), (Clsqnonneg(step:step:end)-MMM)/MMM, 'k', 'filled')
scatter(Tlh(step:step:end), (Clh(step:step:end)-MMM)/MMM, 'k', 'filled')

for u = 1:length(Niter)
    text(Tdamas(u)*1.1, 3*(Cdamas(u)-MMM)/MMM, sprintf("%.0e", Niter(u)))
end
for u = 1:length(niters)
    text(Tquad(u)*0.9, 0.6*(Cquad(u)-MMM)/MMM, sprintf("%.u", niters(u)), 'HorizontalAlignment', 'right')
end
for u = step:step:length(Tlh)
    text(Tlh(u)*0.9, 0.6*(Clh(u)-MMM)/MMM, sprintf("%u", u), 'HorizontalAlignment', 'right')
end

xlabel('Time (s)')
ylabel('$\delta_k$', 'Interpreter', 'LaTeX')

legend('Lawson-Hanson (optimized)', 'cyclic DAMAS', 'random DAMAS', 'quadprog')
ylim([1e-7 1e1])
xlim([8e-3 1000])

grid on



%% Maps: DAMAS

% actual positions
Xgt = [-1.53, -0.8, 0.15, -0.66];
Ygt = [-0.43, -0.45, -0.47, -0.68];
m = 20;
M = max(10*log10(Cbfn));

figure
subplot(4, 1, 2)

set(gcf, 'Position',  [100, 100, 500, 200])

imagesc(xx,yy,reshape(10*log10(xdamas_rand), Ly, Lx))
axis xy
axis image
colormap((hot))

hold on
scatter(Xgt, Ygt, 300, 'w', 'LineWidth', 2)

xlabel("X (m)")
ylabel("Y (m)")
ax = gca;
ax.CLim =[m M];

colorbar
title("(b) DAMAS")
% Maps: CMF-NNLS


subplot(4, 1, 3)

set(gcf, 'Position',  [100, 100, 500, 200])

imagesc(xx,yy,reshape(10*log10(xlh), Ly, Lx))
axis xy
axis image
colormap((hot))
hold on
scatter(Xgt, Ygt, 300, 'w', 'LineWidth', 2)
xlabel("X (m)")
ylabel("Y (m)")
ax = gca;
ax.CLim =[m M];
colorbar
title("(c) CMF-NNLS (opt. LH)")


subplot(4, 1, 4)

set(gcf, 'Position',  [100, 100, 500, 200])

imagesc(xx,yy,reshape(10*log10(xnnls), Ly, Lx))
axis xy
axis image
colormap((hot))
hold on
scatter(Xgt, Ygt, 300, 'w', 'LineWidth', 2)
xlabel("X (m)")
ylabel("Y (m)")
ax = gca;
ax.CLim =[m M];
colorbar
title("(d) DAMAS-NNLS")

% Maps: Beamforming


subplot(4, 1, 1)

set(gcf, 'Position',  [100, 100, 500, 200])

%imagesc(xx,yy,reshape(10*log10(Cbfn), Ly, Lx))

contourf(xx, yy, reshape(10*log10(Cbfn), Ly, Lx), [30:0.5:57],'linewidth', 1)
colormap(hot)
hold on
axis xy
axis image
scatter(Xgt, Ygt, 300, 'kx', 'LineWidth', 2)
colorbar
ax = gca;
ax.CLim =[m+20 M];
xlabel("X (m)")
ylabel("Y (m)")

title("(a) Beamforming")


%%

MMM = Clh(end);

% Objective function in function of computational time
figure
set(gcf, 'Position',  [100, 100, 500, 400])

loglog((Clsqnonneg-MMM)/MMM, 'k-', 'linewidth', 1)
hold on
%scatter(Tlh, Clh, 100, 'black', 'x', 'linewidth', 2)
plot((Clh-MMM)/MMM, 'k--', 'linewidth', 1);

%plscatter(Tlh, Clh, 100, 'black', 'x', 'linewidth', 2)
plot(Niter, (Cdamas-MMM)/MMM, 'k-x', 'linewidth', 2, 'Markersize', 15)
plot(Niter, (Cdamas_rand-MMM)/MMM, 'k--o', 'linewidth', 2, 'Markersize', 10)

plot(niters, (Cquad-MMM)/MMM, 'k-.+', 'linewidth', 2, 'Markersize', 15)
legend('lsqnonneg', 'Lawson-Hanson (optimized)', 'cyclic DAMAS', 'random DAMAS', 'quadprog')
grid on

xlabel('Iteration k')
ylabel('$(J(\mathbf p^k) - J(\mathbf p^\star)) / J(\mathbf p^\star)$', 'Interpreter', 'LaTeX')
