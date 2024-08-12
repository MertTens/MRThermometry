function x = osu_ctimes(k, mask, b1);
dimx = size(k,1);
dimy = size(k,2);
mk = mask.*k;
xmk = ifft2c(mk);
bxmk = conj(b1).*xmk;
x = sum(bxmk,3)*sqrt(dimx*dimy);
x = squeeze(x);
end