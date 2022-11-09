ISP version 1
October 11, 2016

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

This MATLAB collection includes the implementation of the Iterative Sparsification-Projection (ISP) 
algorithms [1] for sparse signal recovery from compressed linear measurements. 

ISP.m			     : Main function
HardThr.m		: Hard thresholding
SoftThr.m		 : Soft thresholding	
SL0Thr.m		 : SL0 shrinkage	
demo.m			 : A demo of how to use the main code

%%%%%%%%%%%%%%%%%%%%%% IMPORTANT NOTES %%%%%%%%%%%%%%%%%%%%%%%

1- The default projection is the one proposed in [2]. However, when the measurement matrix is 
     problematic, e.g., ill-conditioned or low-rank, it is recommended to run
     the algorithm with ADMM projection, though it is more complex. Another case is when the measurement 
     matrix is not so overcomplete, i.e., the number of its rows is close to the number of its columns.
2- For harder problems, i.e., in the case of very incomplete measurements or when the underlying
     signal is assumed to have many non-zeros, it is recommended to set "L"
    (number of inner-loop iterations) to a higher value, say 6.
3- There are lots of shrinkage functions in the literature to be used in the sparsification step of 
     the ISP algorithms. We have included three samples: hard, soft, and SL0 shrinkage. Nevertheless, you can
     wrtie your own thresholding function as an "m" file, like "HardThr.m", and pass it to "ISP.m".
    The choice of the thresholding function highly affects the performance of the algorithm.

References:

   [1] M. Sadeghi and M. Babaie-Zadeh, “Iterative Sparsification-Projection: Fast
         and Robust Sparse Signal Approximation”,
         IEEE Trans. on Signal Proc., vol. 64, no. 21, pp. 5536-5548, November, 2016.

   [2] A. Eftekhari, M. Babaie-Zadeh, C. Jutten, and H. Abrishami-Moghaddam,
        “Robust-SL0 for stable sparse representation in noisy settings,” in Proc.
        ICASSP, 2009, pp. 3433–3436.


  Mostafa Sadeghi
  EE Department, Sharif University of Technology, Tehran, Iran.
  Email: m.saadeghii@gmail.com
