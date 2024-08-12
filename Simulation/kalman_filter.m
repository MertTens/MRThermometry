function [out, outparam] = kalman_filter(obsv, param);

    % Get the guess for the current time point
    pred = param.a.*param.pred;

    % The prediction MSE
    predmse = param.a.^2 .* param.m + param.sigs;
    K = predmse./(predmse + param.sig);
    % The correction
    out = pred + K.*(obsv - pred);
    outparam.a = param.a;
    outparam.pred = out;
    outparam.m = (1-K).*predmse;
    outparam.sig = param.sig;
    outparam.sigs = param.sigs;

end