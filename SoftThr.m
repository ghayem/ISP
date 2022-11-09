function y=SoftThr(x, tau)

% Soft thresholding
y = sign(x).*(abs(x) >= tau).*(abs(x) - tau); 