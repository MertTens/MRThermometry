clear all; close all;
load('recon_for_synaptive/recon.mat');
addpath("Simulation\")
TR = 85.653; TE = 23.75;
alpha = -0.008;
gamma = 2*pi*42.577/1000;
dt = TR*size(data,1)/1000;
meanIm = mean(abs(data),3);
%figure(); imagesc(meanIm); axis equal; axis tight;
mask1 = meanIm > 0.05;
mask1 = imerode(mask1,strel('diamond',1));
cc = bwconncomp(mask1);
m1 = zeros(size(meanIm));
m1(cc.PixelIdxList{1}) = 1;
m2 = zeros(size(meanIm));
m2(45:50,110:115) = 1;
m3 = zeros(size(meanIm));
m3(cc.PixelIdxList{2}) = 1;

mask = m3;

u = data(:,:,1).*m1;
imagesc(angle(u))
figure
imagesc(abs(u))
u = reshape(u,[],1);

% t = data(:,:,1).*mask;
% t = reshape(t,[],1);
% u = [u t];

[rows, cols, pts] = size(data);
d = data.*mask;
ub = reshape(d,[],237);
u = [u ub(:,1)];
% for row = 1:rows
%     for col = 1:cols
%         if(mask(row,col) == 1)
%             temp = zeros(size(mask));
%             temp(row,col) = 1;
%             temp = reshape(temp,[],1);
%             u = [u temp];
%         end
%     end
% end

num_basis = size(u,2);
div = 237;
[dimx, dimy] = size(mask);
kdata = fft2(data);
[dimx, dimy, time_pts] = size(kdata);
nspokes = 20;
[usmat] = cartesian_undersampling(dimx, dimy, nspokes, time_pts);
usmat(find(usmat > 1)) = 1;
usmat = pagectranspose(usmat);
% usmat = zeros(size(data));
% usmat(35:55,:,:) = 1;
% usmat = circshift(usmat, dimx/2, 1);
% usmat = circshift(usmat, dimy/2, 2);
% [~, usmat] = radial_undersampling(img_dim, num_spokes, start_angle, increment)
masku = usmat;
kdatau = kdata.* masku;
kdatag = kdatau;
maskg = masku;
a_factor = dimx*dimy / sum(sum(usmat(:,:,1)));
fprintf("The acceleration factor is %f\n", a_factor)


svd_recon = zeros(dimx, dimy, div);
ku = fft2(reshape(u, dimx, dimy,1, []));
ku = reshape(ku, [], num_basis);
vec_usmat = reshape(maskg, [], div);
vec_observations = reshape(kdatag, [], div);
for m = 1:div
    testu = zeros(size(ku));
    for n = 1:num_basis
        testu(:,n) = ku(:,n).*vec_usmat(:,m);
    end
    testy = vec_observations(:,m);

    testu = reshape(nonzeros(testu),[],num_basis);
    testy = nonzeros(testy);
    coeffs = inv(testu'*testu)*testu'*testy;

    img = u(:,1:num_basis)*coeffs;
    img = reshape(img,dimx,dimy);
    svd_recon(:,:,m) = (img);
end

proposed = svd_recon;
TR = 85.653; TE = 23.75;
alpha = -0.008;
gamma = 2*pi*42.577/1000;
dt = TR*size(data,1)/1000;
meanIm = mean(abs(data),3);
%figure(); imagesc(meanIm); axis equal; axis tight;
mask1 = meanIm > 0.05;
mask1 = imerode(mask1,strel('diamond',1));
cc = bwconncomp(mask1);
m1 = zeros(size(meanIm));
m1(cc.PixelIdxList{1}) = 1;
m2 = zeros(size(meanIm));
m2(45:50,110:115) = 1;
m3 = zeros(size(meanIm));
m3(cc.PixelIdxList{2}) = 1;

x = linspace(-1,1,size(data,1));
y = linspace(-1,1,size(data,2));
[X,Y] = meshgrid(x,y);
N = numel(X);
model = [ones(N,1), X(:), Y(:), X(:).*Y(:)];
model_fit = model(cc.PixelIdxList{1},:);

imagesc(imfuse(meanIm,m1));

pdiff_data = angle(conj(data(:,:,2:end)).*data(:,:,1));
pdiff_prop = angle(conj(proposed(:,:,2:end)).*proposed(:,:,1));

nTime = size(pdiff_data,3);
t = dt:dt:nTime*dt;
vals1 = zeros(1,nTime);
vals2 = zeros(1,nTime);
temp = zeros(size(vals1));
temp1 = zeros(size(pdiff_data));
temp2 = zeros(size(temp1));
for i = 1:nTime
    pdiff_data2D = CorrectPhase(pdiff_data(:,:,i),model,cc.PixelIdxList{1}, pdiff_data(:,:,i));
    vals1(i) = mean(pdiff_data2D(m2>0))./(alpha*gamma*TE);
    
    pdiff_prop2D = CorrectPhase(pdiff_prop(:,:,i),model,cc.PixelIdxList{1}, pdiff_prop(:,:,i));
    vals2(i) = mean(pdiff_prop2D(m2>0))./(alpha*gamma*TE);
    a = pdiff_data(:,:,i);
    temp(i) = mean(a(m2>0));
    q = 0.5;
    subplot(2,2,1)
    imagesc(pdiff_data(:,:,i).*mask1, [-q q])
    subplot(2,2,2)
    imagesc(pdiff_prop(:,:,i).*mask1, [-q q])
    subplot(2,2,3)
    imagesc(pdiff_data2D.*mask1, [-q q])
    subplot(2,2,4)
    imagesc(pdiff_prop2D.*mask1, [-q q])
    pause(0.01)
end
figure(); plot(t,vals1-vals1(1),'linewidth',2); hold on; plot(t,vals2-vals2(1),'linewidth',2);
xlabel('Time (s)'); ylabel('Temperature ({\circ}C)'); legend('Original','Proposed');
title('Mean temperature over ROI');
figure();
clims = [-10,10];
subplot(1,2,1); imagesc(pdiff_data2D(30:65,90:130).*m3(30:65,90:130)./(alpha*gamma*TE),clims); axis equal; axis tight; 
hold on; rectangle('Position',[20,15,5,5],'EdgeColor','r','linewidth',2'); hold off;
subplot(1,2,2); imagesc(pdiff_prop2D(30:65,90:130).*m3(30:65,90:130)./(alpha*gamma*TE),clims); axis equal; axis tight;
p = get(gca,'Position');
a = colorbar;
ylabel(a,'Temperature ({\circ}C)');
set(gca,'Position',p);
set(gcf,'Position',[488,342,927,420]);
sgtitle('Temperature at last time point');


obsv1 = vals1 - vals1(1);
kf1 = zeros(size(obsv1));
param.a = 1;
param.pred = vals1(1);
param.sig = 1;
param.m = param.sig;
param.sigs = 0.001;

for n = 2:max(size(obsv1))
    [kf1(n), param] = kalman_filter(obsv1(n), param);
end

obsv2 = vals2 - vals2(1);
kf2 = zeros(size(obsv2));
param.a = 1;
param.pred = vals1(1);
param.sig = 2;
param.m = param.sig;
param.sigs = 0.001;

for n = 2:max(size(obsv2))
    [kf2(n), param] = kalman_filter(obsv2(n), param);
end


figure
plot(obsv1,'linewidth',2)
hold all
plot(kf2,'linewidth',2)
legend("original", "9x accelerated, proposed")


function corr = CorrectPhase(input,model,maskInds, input2)
model_fit = model(maskInds,:);
data = input2(maskInds);
fit_params = mldivide(model_fit,data);
model_apply = reshape(model*fit_params,size(input));

corr = input-model_apply;
% corr = input;
% corr = model_apply;
end