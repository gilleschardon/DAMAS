
%% Lawson-Hanson on experimental data
%
% 3D grid of 3.75e5 points

% Approx 45min

% load the data
load damas3D
addpath('..')

% we use all microphones
N = Pmic(:, 2).^2 + Pmic(:, 1).^2;
[~, order] = sort(N);
Nm = 128;
Z = order(1:Nm);
Pmic = Pmic(Z, :);
Data = Data(Z, Z);

close all

% 3D grid
Lx = 150;
Ly = 50;
Lz = 50;
xx = linspace(-2, 1, Lx)';
yy = linspace(-1, 0, Ly)';
zz = linspace(4, 5, Lz)';
[Xg, Yg, Zg] = meshgrid(xx, yy, zz);

D = dictionary(Pmic, [Xg(:) Yg(:) Zg(:)], k);
%% CMF-NNLS, Lawson-Hanson algorithm

tic;[xlh, Tlh, Clh] = cmf_nnls(D, Data, 0); Tnd = toc;
%%
%% DAMAS-NNLS, Lawson-Hanson algorithm

% we compute the beamforming first (no normalization)
Cbf = sum(conj(D) .* (Data*D), 1);
Cbf = real(Cbf');

tic;[xsn, Tsn, Csn] = damas_nnls(D, Cbf, 0); Tdn = toc;

clear D
save 3Dexp


%% 3D plot

CAM = [8.2276   -7.7357    6.5083];
figure


% actual positions and powers
Xgt = [-1.54, -0.8, 0.15, -0.66];
Ygt = [-0.43, -0.45, -0.48, -0.68];
Zgt = [4.39 4.41 4.43 4.43];

pp = [63.8 62.5 61.1 63];

subplot(3, 1, 1)

% color
C = [0.85,0.85,0.85]; % gray for the projections

scatter3(Xgt, Zgt, Ygt, 10.^(pp/10)/10000, pp, 'filled')
colormap(flipud(summer))
colorbar
ax = gca;
ax.CLim =[40 62];

campos(CAM)

title('(a) Sources')
hold on

% projections (shaded)
scatter3(Xgt, Zgt*0+5, Ygt, 10.^(pp/10)/10000, C, 'filled')
scatter3(Xgt, Zgt, Ygt*0-1, 10.^(pp/10)/10000, C, 'filled')



axis equal
xlim([-2, 1])
ylim([4, 5])
zlim([-1, 0])

xlabel("X")
ylabel("Z")
zlabel("Y")


subplot(3, 1, 2)
% CMF-NNLS
% we threshold the powers for clarity
support = xlh(:) > max(xlh(:))/100;

scatter3(Xg(support), Zg(support), Yg(support), (xlh(support)+eps)/10000, 10*log10(xlh(support)), 'filled')
colormap(flipud(summer))
colorbar
ax = gca;
ax.CLim =[40 62];
campos(CAM)

title('(b) CMF-NNLS')
hold on

% projections (shaded)
scatter3(Xg(support), Zg(support)*0+5, Yg(support), (xlh(support)+eps)/10000, C, 'filled')
scatter3(Xg(support), Zg(support), Yg(support)*0-1, (xlh(support)+eps)/10000, C, 'filled')

axis equal
xlim([-2, 1])
ylim([4, 5])
zlim([-1, 0])

xlabel("X")
ylabel("Z")
zlabel("Y")



subplot(3, 1, 3)
% DAMAS-NNLS
support = xsn(:) > max(xsn(:))/100;


scatter3(Xg(support), Zg(support), Yg(support), (xsn(support)+eps)/10000, 10*log10(xsn(support)), 'filled')
colormap(flipud(summer))
colorbar
ax = gca;
ax.CLim =[40 62];
campos(CAM)
hold on

% projections (shaded)
scatter3(Xg(support), Zg(support)*0+5, Yg(support), (xsn(support)+eps)/10000, C, 'filled')
scatter3(Xg(support), Zg(support), Yg(support)*0-1, (xsn(support)+eps)/10000, C, 'filled')

title('(c) DAMAS-NNLS')

axis equal
xlim([-2, 1])
ylim([4, 5])
zlim([-1, 0])

xlabel("X")
ylabel("Z")
zlabel("Y")

set(gcf,'Renderer','Painter')

%%

figure
supportlh = xlh(:) > 0;
supportdn = xsn(:) > 0;


% actual positions and powers
Xgt = [-1.54, -0.8, 0.15, -0.66];
Ygt = [-0.43, -0.45, -0.48, -0.68];
Zgt = [4.39 4.41 4.43 4.43];

pp = [63.8 62.5 61.1 63];

subplot(3, 2, 1)

scatter(Xgt, Ygt, 50,'xk', 'linewidth', 2)

ax = gca;
ax.CLim =[40 62];

axis equal
xlim([-2, 1])
ylim([-1, 0])

xlabel("X")
ylabel("Y")

title('(a) Sources')

colorbar

subplot(3, 2, 2)

scatter(Xgt, Zgt, 50,'xk', 'linewidth', 2)

ax = gca;
ax.CLim =[40 62];

axis equal
xlim([-2, 1])
ylim([4, 5])

xlabel("X")
ylabel("Y")
colorbar

subplot(3, 2, 3)

scatter(Xg(supportlh), Yg(supportlh), (xlh(supportlh)+eps)/10000, 10*log10(xlh(supportlh)), 'filled')
hold on
scatter(Xgt, Ygt, 'xw')

colorbar

ax = gca;
ax.CLim =[40 62];

axis equal
xlim([-2, 1])
ylim([-1, 0])

xlabel("X")
ylabel("Y")

title('(a) CMF-NNLS')

subplot(3, 2, 4)

scatter(Xg(supportlh), Zg(supportlh), (xlh(supportlh)+eps)/10000, 10*log10(xlh(supportlh)), 'filled')
hold on
scatter(Xgt, Zgt, 'xw')

colorbar

ax = gca;
ax.CLim =[40 62];

axis equal
xlim([-2, 1])
ylim([4, 5])

xlabel("X")
ylabel("Z")


subplot(3, 2, 5)

scatter(Xg(supportdn), Yg(supportdn), (xsn(supportdn)+eps)/10000, 10*log10(xsn(supportdn)), 'filled')
hold on
scatter(Xgt, Ygt, 'xw')

colorbar

ax = gca;
ax.CLim =[40 62];

axis equal
xlim([-2, 1])
ylim([-1, 0])

xlabel("X")
ylabel("Y")
title('(b) DAMAS-NNLS')

subplot(3, 2, 6)

scatter(Xg(supportdn), Zg(supportdn), (xsn(supportdn)+eps)/10000, 10*log10(xsn(supportdn)), 'filled')
hold on
scatter(Xgt, Zgt, 'xw')

colorbar

ax = gca;
ax.CLim =[40 62];

axis equal
xlim([-2, 1])
ylim([4, 5])

xlabel("X")
ylabel("Z")
colormap(flipud(gray))

