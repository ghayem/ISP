function y=SL0Thr(x, tau)

% SL0 thresholding
y = x.*(1-exp(-x.^2/tau^2)); 