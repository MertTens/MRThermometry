function [b1] = get_b1(kdata)

% assumes 4th dimension is time
kdata = mean(kdata,4);
coil_imgs = ifft2c(kdata);
combined_img = sqrt(sum(abs(coil_imgs.^2),3));
combined_img = reshape(combined_img,size(combined_img,1),size(combined_img,2),1);
b1 = coil_imgs./combined_img;

end