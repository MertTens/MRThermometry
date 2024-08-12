function k = osu_times(x,mask,b1);
dimx = size(x,1);
dimy = size(x,2);
x = reshape(x,size(x,1), size(x,2),1,size(x,3));
bx = b1.*x;
kbx = fft2c(bx);

k = mask.*kbx;
k = squeeze(k)./sqrt(dimx*dimy);
end