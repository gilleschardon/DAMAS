

%% Compares DAMAS, CMF-NNLS, DAMAS-NNLS
%% on simulated data
%

addpath('..')
% we do not use the data, just the array geometry
load damas3D

N = Pmic(:, 2).^2 + Pmic(:, 1).^2;
[~, order] = sort(N);
Nm = 64;
Z = order(1:Nm);
Pmic = Pmic(Z, :);
Data = Data(Z, Z);

close all

Lx = 100;
Ly = 100;
Lz = 1;
xx = linspace(0, 1, Lx+1)';
yy = linspace(0, 1, Ly+1)';
xx = xx(1:end-1);
yy = yy(1:end-1);
zz = 4.4;
[Xg, Yg, Zg] = meshgrid(xx, yy, zz);


D = dictionary(Pmic, [Xg(:) Yg(:) Zg(:)], k);

% indexes of the sources

idx = [6050 6767 2545];
% index XXYY = coordinates (0.XX, 0.YY, 4.4)

% sources power 
% /!\ variances of real and imaginary part, power is multiplied by 2
sigma_source1 = 1e8;
sigma_source2 = 4e8;
sigma_source3 = 2e7;

% noise variances (first and last 32 sensors)
% same remark
sigma_noise1 = 4e6;
sigma_noise2 = 1.6e7;

% number of snapshots
N = 500;

% large DAMAS matrix
DD = abs(D'*D).^2;

% normalized beamforming matrix
Dbf = D ./ sum(abs(D).^2, 1);

% masks for estimating the powers, we sum on the source and its 8
% neighbouring pixels
[r1, c1] = ind2sub([Ly, Lx], idx(1));
[r2, c2] = ind2sub([Ly, Lx], idx(2));
[r3, c3] = ind2sub([Ly, Lx], idx(3));

mask1 = zeros(Ly, Lx);
mask1(r1-1:r1+1, c1-1:c1+1) = 1;
mask1 = logical(mask1);

mask2 = zeros(Ly, Lx);
mask2(r2-1:r2+1, c2-1:c2+1) = 1;
mask2 = logical(mask2);

mask3 = zeros(Ly, Lx);
mask3(r3-1:r3+1, c3-1:c3+1) = 1;
mask3 = logical(mask3);


% sample size
Nt = 200;


% powers of the sources

% beamforing with diagonal removal
PBF = zeros(4, Nt);
% DAMAS without/with diagonal removal
PD = zeros(4, Nt);
% LH without/with noise estimation
PLH = zeros(4, Nt);
% DAMAS-NNLS
PDN = zeros(4, Nt);

JJ = ones(64)-eye(64);

% LH with noise estimation, estimated noise variances
LHnoiseest = zeros(64, Nt);
%%
% %Uncomment to compute the large matrix for BF and DAMAS with diagonal removal
% DDdr = zeros(size(DD));
% 
% for u = 1:size(DD, 1)
% 
%     dd = diag(diag(D(:, u)*D(:, u)'));
%     
%     DDdr(u, :) = DD(u, :) - sum(conj(D).* (dd*D), 1);
% 
% end
% 
% save DDdr DDdr
%% otherwise load from a precomputed matrix
load DDdr

% DAMAS iterations
Niter = 1e7;

% beamforming steering vectors, diagonal removal, normalized
C1 = D(:, idx(1)) * D(:, idx(1))';
C1 = C1 - diag(diag(C1));
C1 = C1 / sum(abs(C1(:)).^2);

C2 = D(:, idx(2)) * D(:, idx(2))';
C2 = C2 - diag(diag(C2));
C2 = C2 / sum(abs(C2(:)).^2);

C3 = D(:, idx(3)) * D(:, idx(3))';
C3 = C3 - diag(diag(C3));
C3 = C3 / sum(abs(C3(:)).^2);

for u = 1:Nt
    
    waitbar(u/Nt)

    %% simulation of the data
    sigs =  D(:, idx) * (sqrt([sigma_source1 sigma_source2 sigma_source3]') .* (randn(length(idx), N) + 1i * randn(length(idx), N))); 
    noise = sqrt([ones(32, 1) * sigma_noise1 ; ones(32, 1) * sigma_noise2]) .* (randn(size(D, 1), N) + 1i * randn(size(D, 1), N));
    Data = (sigs + noise) * (sigs + noise)' /N;

    SNR = 20*log10(norm(sigs, 'fro')/norm(noise, 'fro'));

%% Beamforming
    Cbf = sum(conj(D) .* (Data*D), 1);
    Cbf = real(Cbf');

% normalized
    Cbfp = sum(conj(Dbf) .* (Data*Dbf), 1);
    Cbfp = real(Cbfp');

    
    
     PBF(1, u) = real(trace(Data*C1));
     PBF(2, u) = real(trace(Data*C2));
     PBF(3, u) = real(trace(Data*C3));



% diagonal removal, input for DAMAS
    DataDR = Data - diag(diag(Data));
    CbfDR = sum(conj(D) .* (DataDR*D), 1);
    CbfDRdamas = real(CbfDR');


%% optimized NNLS DR

tic
xlh = cmf_nnls_dr(D, DataDR, 0);

Tlh(u) = toc;
Clh(u) = norm(JJ.*(D*diag(xlh)*D' - Data), 'fro')^2;

% powers of the three sources and total power
PLH(1, u) = sum(xlh(mask1));
PLH(2, u) = sum(xlh(mask2));
PLH(3, u) = sum(xlh(mask3));
PLH(4, u) = sum(xlh);


%% DAMAS, diagonal removal

tic
xd = damas_rand(real(DDdr), real(CbfDRdamas), Niter, zeros(size(DDdr, 1), 1));

Cdamas(u) = norm(JJ.*(D*diag(xd)*D' - Data), 'fro')^2;
Td(u) = toc;

PD(1, u) = sum(xd(mask1));
PD(2, u) = sum(xd(mask2));
PD(3, u) = sum(xd(mask3));
PD(4, u) = sum(xd);

%% DAMAS-NNLS
tic
xdn = damas_nnls_dr(D, real(CbfDRdamas), 0);
Tdn(u) = toc;

PDN(1, u) = sum(xdn(mask1));
PDN(2, u) = sum(xdn(mask2));
PDN(3, u) = sum(xdn(mask3));
PDN(4, u) = sum(xdn);


end

clear DDdr DD D Dbf

save 2Dsim

%%
close all

bins = (-15:0.5:25);

color1 = [0.4 0.4 0.4];
color2 = color1; 
color3 = color2;

%%
m = 60;
M = 90;

BBB = [40 90];

figure

subplot(2, 2, 1)

contourf(xx, yy, 10*log10(reshape(Cbfp, Lx, Ly)), 60:1:90)
hold on
scatter(Xg(idx), Yg(idx), 100,'kx', 'LineWidth', 2)

axis xy
axis square
colormap((hot))
colorbar
ax = gca;
ax.CLim =[m M];
title("(a) Beamforming")

subplot(2, 2, 2)

imagesc(xx, yy, 10*log10(reshape(xd, Lx, Ly)))
hold on
scatter(Xg(idx), Yg(idx), 300,'wo', 'LineWidth', 1)

axis xy
axis square
colormap((hot))
colorbar
ax = gca;
ax.CLim =BBB;
title("(b) DAMAS")

subplot(2, 2, 3)

imagesc(xx, yy, 10*log10(reshape(xlh, Lx, Ly)))
hold on
scatter(Xg(idx), Yg(idx), 300,'wo', 'LineWidth', 1)

axis xy
axis square
colormap((hot))
colorbar
ax = gca;
ax.CLim =BBB;
title("(c) CMF-NNLS")

subplot(2, 2, 4)

imagesc(xx, yy, 10*log10(reshape(xdn, Lx, Ly)))
hold on
scatter(Xg(idx), Yg(idx), 300,'wo', 'LineWidth', 1)

axis xy
axis square
colormap((hot))
colorbar
ax = gca;
ax.CLim =BBB;
title("(d) DAMAS-NNLS")

%%
figure

B1 = [-3 1];
B2 = [-4 1];
B3 = [-8 1];
B4 = [-3 0];

subplot(3, 4, 1)
plothist(PD(1, :), sigma_source1*2, B1, bins, color2, '', 'DAMAS')
subplot(3, 4, 2)
plothist(PD(2, :), sigma_source2*2, B2, bins, color2, '', '')
subplot(3, 4, 3)
plothist(PD(3, :), sigma_source3*2, B3, bins, color2, '', '')
subplot(3, 4, 4)
plothist(sum(PD(1:3, :)),PD(4, :),  B4, bins, color2, '', '')


subplot(3, 4, 5)
plothist(PLH(1, :), sigma_source1*2, B1, bins, color2, '', 'CMF-NNLS')
subplot(3, 4, 6)
plothist(PLH(2, :), sigma_source2*2, B2, bins, color2, '', '')
subplot(3, 4, 7)
plothist(PLH(3, :), sigma_source3*2, B3, bins, color2, '', '')
subplot(3, 4, 8)
plothist(sum(PLH(1:3, :)),PLH(4, :),  B4, bins, color2, '', '')


subplot(3, 4, 9)
plothist(PDN(1, :), sigma_source1*2, B1, bins, color2, '\Delta_1', 'DAMAS-NNLS')
subplot(3, 4, 10)
plothist(PDN(2, :), sigma_source2*2, B2, bins, color2, '\Delta_2', '')
subplot(3, 4, 11)
plothist(PDN(3, :), sigma_source3*2, B3, bins, color2, '\Delta_3', '')
subplot(3, 4, 12)
plothist(sum(PDN(1:3, :)), PDN(4, :), B4, bins, color2, '\Delta_0', '')

set(gcf,'Renderer','Painter')
