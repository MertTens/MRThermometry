function [mask] = gaussian_mask(dimx, dimy, factor);

    varx = dimx*1.5;
    vary = dimy*1.5;

    pts = round(dim*dim/factor);
    mask = zeros(dim,dim);
    R = mvnrnd([dim/2 dim/2],[var 0; 0 var],pts*100);

    R = round(R);

    m = 0;
    for n = 1:pts*100
        x = R(n,1);
        y = R(n,2);
        if x > dim || x < 1
            continue
        end

        if y > dim || y < 1
            continue
        end

        if mask(x,y) == 1
            continue;
        else
            mask(x,y) = 1;
            m = m + 1;
        end

        if m == pts
            break;
        end
    end
end