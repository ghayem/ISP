
% A simple demo of how to use the "ISP" function for recovery of sparse signals from
%  their compressed (possibly noisy) measurements.

n =200;         % sparse signal length
m =40;          % number of measurements

S=2:2:20;      % number of non-zeros

numiter=100;        % number of Monte-Carlo simulations

% ISP parameters
maxiter=300;             % maximum number of outer-loop iterations
tauf=5e-6;                   % minimum value of the threshold
opts.c=0.9;                  % Threshold decaying factor (default: c=0.9)
opts.L =3;                     %  number of inner-loop iterations (default: L=3)
opts.gam=0.4;            %  step-size of the projection step (default: gam=0.4)
opts.proj_mode=0;    % set to "1" for using ADMM projection and "0" for applying Eftekhari's projection algorithm

noise_std=0;               % noise standard deviation

Ns=length(S);
mse_hard=zeros(Ns, 1);
mse_soft=zeros(Ns, 1);
mse_sl0=zeros(Ns, 1);

p=1;

% For your experiments, you can use "matlabpool" and "parfor" structures to speed up
% the Monte-Carlo simulations

for s=S   % number of non-zeros of x
    
    mse1=0;
    mse2=0;
    mse3=0;
    
    disp(['Monte-Carlo for s= ', num2str(s),' ...']);
    
    for iter=1:numiter
        
        A = randn(m, n);         % measurement matrix
        A_pinv=pinv(A);
        
        nn=noise_std*randn(m,1);
        xo = zeros(n,1);
        xo(randsample(n,s))= randn(s,1);
        y=A*xo+nn;
        
        e=noise_std*sqrt(m);        % noise power
        
        %% ISP-Hard
        T='HardThr';                  % thresholding function
        x_hard=ISP(A, y, T, e, maxiter, tauf, opts);
        
        %% ISP-Soft
        T='SoftThr';                  % thresholding function
        x_soft=ISP(A, y, T, e, maxiter, tauf, opts);
        
        %% ISP-SL0
        T='SL0Thr';                   % thresholding function
        x_sl0=ISP(A, y, T, e, maxiter, tauf, opts);
        

        mse1=mse1+20*log10(norm(x_hard-xo)/norm(xo));
        mse2=mse2+20*log10(norm(x_soft-xo)/norm(xo));
        mse3=mse3+20*log10(norm(x_sl0-xo)/norm(xo));
        
    end
    
    mse_hard(p)= mse1/numiter;
    mse_soft(p)= mse2/numiter;
    mse_sl0(p)= mse3/numiter;
    
    p=p+1;
    
end

%% Plot the results

plt=plot(S, mse_hard,'b',S, mse_soft,'--m',S, mse_sl0,'r');set(gca,'fontsize',15);grid on
set(plt,'linewidth',3);xlabel('s','fontsize',18);ylabel('MSE (dB)','fontsize',18);
leg=legend('ISP-Hard','ISP-Soft','ISP-SL0');
set(leg,'fontsize',18);
