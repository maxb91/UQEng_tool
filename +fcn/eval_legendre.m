function P = eval_legendre( X, deg )
% Evaluate the Legendre polynomial of degree deg
% for x-values given in a vector X

% Make sure that the degree is valid
if deg < 0
    error('The degree should be at least 0!')
end

%% Compute the unnormalized Legendre polynomials (Part 1.1)
% case distinction
if deg == 0
    P = ones(size(X));
elseif deg == 1
    P = X;
else
    % Recursion rule:
    % (n+1)*P_nplus1 = (2*n+1)*X*P_n - n*P_nminus1;
    %P_nplus1 = (2*n+1)/(n+1)*X*P_n  - n/(n+1)*P_nminus1;
    
    % initialize recursion
    P_nminus1 = ones(size(X));
    P_n = X;
    
    for n = 1:(deg-1)
        % Recursion rule for deg >= 2
        % ".*" for element-wise multiplication
        P_nplus1 = (2*n+1)/(n+1)*X.*P_n  - n/(n+1)*P_nminus1;
        % assign values for the next iteration
        P_nminus1 = P_n; % P_nminus1 must be assigned first!
        P_n = P_nplus1;
    end

    
    % P contains the evaluations of the unnormalized Legendre polynomial of
    % degree deg
    P = P_nplus1; 
end

%% Compute the normalized Legendre polynomials (Part 1.2)
P = sqrt(2*deg+1)*P;
