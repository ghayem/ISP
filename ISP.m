function x=ISP(A, y, T, e, maxiter, tauf, opts)
%
%  This function implements the ISP algorithms proposed in [1] which aim at recovery
%   of sparse signals from their compressed linear measurements by solving the
%   following problem:
%
%    min_x  f(x)   s.t.    ||y-Ax||_2 <=  e          (P)
%
%   where, "f" is a smooth or non-smooth sparsity-inducing function, "y" is an (m*1)
%   vector of linear measurements, "A" is an (m*n) measurement matrix, and "e >=0" 
%   is the noise power.
%
%  As its name suggests, the family of ISP algorithms follow an iterative
%  sparsification-projection approach to get an approximate solution of (P).
%  The sparsification step can be realized by a thresholding function
%  resulting from either the proximal mapping of "f" when it is easy-to-compute,
%  or one-step gradient descent on "f" when it is smooth.
%
%  Syntax:
%
%               x = ISP(A, y, e, maxiter, tauf);                         default version
%               x = ISP(A, y, e, maxiter, tauf, opts);               full version
%
% %%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%
%
%  Required:
%
%                       y:                                     m*1 measurement vector
%                       A:                                    m*n measurement matrix
%                       e:                                    noise power (set e=0 for noiseless case)
%                       maxiter:                        maximum number of outer-loop iterations
%                       tauf:                               minimum value of the threshold
%                       T:                                    thresholding function, for instance, hard/soft
%                                                               sharinkage or one-step gradient descent of a
%                                                               smooth sparsity promoting function
%                                                               **** See "ReadMe" ****
% Optional:
%
%                        opts.c:                           Threshold decaying factor (default: c=0.9)
%                        opts.L :                           number of inner-loop iterations (default: L=5)
%                        opts.A_pinv :                pseudo inverse of the measurment matrix
%                        opts.gam :                    step-size of the projection step (default: gam=0.4)
%                        opts.proj_mode :        set to "1" to use ADMM projection and "0"
%                                                                to apply the projection algorithm
%                                                                proposed in [2] (default: proj_mode=0)
%                                                                **** See "ReadMe" ****
%
%
%
% %%%%%%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%%%%%%%
%
% The algorithm returns "x", an estimate of the solution of (P).
%
%
%
% Reference:
%
%   [1] M. Sadeghi and M. Babaie-Zadeh, “Iterative Sparsification-Projection: Fast
%         and Robust Sparse Signal Approximation”,
%         IEEE Trans. on Signal Proc., vol. 64, no. 21, pp. 5536-5548, November, 2016.
%
%   [2] A. Eftekhari, M. Babaie-Zadeh, C. Jutten, and H. Abrishami-Moghaddam,
%        “Robust-SL0 for stable sparse representation in noisy settings,” in Proc.
%        ICASSP, 2009, pp. 3433–3436.
%


%  Mostafa Sadeghi
%  EE Department, Sharif University of Technology, Tehran, Iran.
%  Email: m.saadeghii@gmail.com

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 6
    error('Not enough input parameters!');
end

% parse input parameters

if (isfield(opts,'c'))
    c=opts.c;
else
    c=0.9;
end

if (isfield(opts,'L'))
    L=opts.L;
else
    L=5;
end

if (isfield(opts,'A_pinv'))
    A_pinv=opts.A_pinv;
else
    A_pinv=pinv(A);
end

if (isfield(opts,'gam'))
    gam=opts.gam;
else
    gam=0.4;
end

if (isfield(opts,'proj_mode'))
    proj_mode=opts.proj_mode;
else
    proj_mode=0;
end

x = A_pinv*y;            % signal initialization by minimum l_2 norm solution of y=Ax

tau=3*max(abs(x));  % initial value of the threshold
e2=e^2;

if proj_mode
    K=size(A,2);
    M=inv(eye(K)+gam*(A'*A));
    lam=zeros(size(y));
end

if e < 1e-10
    noise_mode=0;     % noiseless recovery
else
    noise_mode=1;
end

% Main Loop

iter=1;

if noise_mode==0
    
    while iter <= maxiter || tau >= tauf
        
        for l=1:L
            x = feval(T, x, tau);                             % sparsification step
            x = x - A_pinv*(A*x-y);                      % projection step
        end
        
        tau=c*tau;
        iter=iter+1;
        
    end
    
else
    
    if proj_mode==0
        
        while iter <= maxiter || tau >= tauf
            
            for l=1:L
                x = feval(T, x, tau);                             % sparsification step
                
                zx=y-A*x;
                if (zx'*zx) > e2                                     % Eftekhari's projection step
                    x = x - A_pinv*(A*x-y);                      
                end
                
            end
            
            tau=c*tau;
            iter=iter+1;
            
        end
        
    elseif proj_mode==1
        
        while iter <= maxiter || tau >= tauf
            
            for l=1:L
                x0 = feval(T, x, tau);                             % sparsification step
                
                zx=y-A*x0;
				x = x0;
				
                while (zx'*zx) > e2                            % proposed projection step
                    
                    % z-update
                    z=y-A*x+1/gam*lam;
                    zn=sqrt(z'*z);
                    
                    if zn > e
                        z=z/zn *e;                         % z projection to A_z set
                    end
                    
                    % x-update
                    G=y-z+1/gam*lam;
                    x=M*(x0+gam*(A'*G));
                    
                    zdif=z-y+A*x;
                    
                    % lam-update
                    lam=lam-gam*zdif;
                    zx=z-zdif;
                    
                end
                
            end
            
            tau=c*tau;
            iter=iter+1;
            
        end
        
    end
    
end