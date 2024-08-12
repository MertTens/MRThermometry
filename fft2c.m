function [x] = fft2c(x)
x = ifftshift(fft(fftshift(x,1),[],1),1);
x = ifftshift(fft(fftshift(x,2),[],2),2);
end