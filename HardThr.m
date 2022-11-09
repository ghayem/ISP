function y=HardThr(x, tau)

% Hard thresholding
y = x.*(abs(x) >= tau);