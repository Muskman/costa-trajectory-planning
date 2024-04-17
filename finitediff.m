% Code for generating a finite difference matrix using forward and
% backwards differentiation rules


function [A, B, D] = finitediff(n,opts)
    
    A = zeros(n,n); B = zeros(n,n); D = zeros(n,n);
    p = opts.goal_weight;

    A(1,1) = 1; A(n,n) = 1-p*(n/n);
    A(1,2) = -1; A(n,n-1) = -(1-p*(n/n));
    
    B(1,1) = 1; 
    B(n,n-1) = -(1-p*(n/n));
    
    D(1,1) = 1; D(n,n) = 1;
    D(1,2) = -1;

    for j = 2:(n / 2)
        A(j,j-1) = -(1-p*j/n);
        A(j,j) = 2-p*(j+j+1)/n;
        A(j,j+1) = -(1-p*(j+1)/n);
        
        A(n-j+1,n-j+2) = -(1-p*(n-j+2)/n);
        A(n-j+1,n-j+1) = 2-p*(n-j+1+n-j+2)/n;
        A(n-j+1,n-j) = -(1-p*(n-j+1)/n);

        B(j,j-1) = -(1-p*j/n);
        B(j,j) = 1-p*(j+1)/n;
        
        B(n-j+1,n-j+1) = 1-p*(n-j+2)/n;
        B(n-j+1,n-j) = -(1-p*(n-j+1)/n);

        D(j,j) = 1;
        D(j,j+1) = -1;
        
        D(n-j+1,n-j+2) = -1;
        D(n-j+1,n-j+1) = 1;
    end
    
    if mod(n, 2) == 1
        index = int8(n / 2);
        A(index,index-1) = -(1-p*index/n);
        A(index,index) = 2-p*(index+index+1)/n;
        A(index,index+1) = -(1-p*(index+1)/n);

        B(index,index-1) = -(1-p*index/n);
        B(index,index) = (1-p*(index+1)/n);

        D(index,index) = 1;
        D(index,index+1) = -1;
    end
    
    A = 2*kron(A,eye(2)); B = 2*kron(B,eye(2));
    D = kron(D,eye(2));

end