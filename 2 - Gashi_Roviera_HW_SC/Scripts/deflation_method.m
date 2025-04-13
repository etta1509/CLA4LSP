function [eigvalues, eigvects, residualnorms] = deflation_method(A, v, eigvects, eigvalues, M, maxIter, relTol)

    % Function  that computes the M smallest eigenvalues and eigenvectors of a
    % symmetrix matrix 
    % Inputs:
    % A: matrix whose values will be computed
    % v: starting vector for the inverse power method
    % eigvects, eigvalues : contain the smallest eigenvalue and eigenvector
    % already computed
    % M: number of values to be computed 
    % maxIter : max number of iterations inside the inverse power method
    % relTol : minimum relative tolreance to stop inverse power method
    % iterations
    % Outputs:
    % eigvects, eigvalues : smallest M eigenvalues and corresponding
    % eigenvectors
    % residualnorms: vector containing the norm of A*x - lambda*x
    
    % Number of points
    N = size(v,1);
    
    % Vector of all zeros except the first entry
    e1 = zeros(N,1);
    e1(1,1) = 1;
    
    % Construct the j-th Householder matrix
    Pj = eye(N,N) - 2*((eigvects(:, 1) + e1)*((eigvects(:, 1) + e1)')) /norm(eigvects(:, 1) + e1)^2;
    % Construct the similar matrix Bj to A
    Bj = Pj * A * Pj;
    
    % inizialize the vector of residual norms
    residualnorms = zeros(M,1);
    residualnorms(1) = norm(A * eigvects(:,1) - eigvalues(1) * eigvects(:,1));
    
    %Inizialize the matrix S, tha will be at each j step: P1*P2*...*Pj_1
    S = eye(N);
    
    % Loop to caclulate the eigenvalues and eigenvectors
    for j = 2: M
        S = S * Pj;
    
        % Submatrix containg the remaining eigenvalues
        Aj =  Bj(j:end, j:end);
        
        % Find the next smallest eigenvalue and the eigenvector of Bj by using the inverse power method
        [eigvalues(j), x_bar] = invpower_method(Aj, v(j:end), maxIter, relTol); 
        
        % Compute the sumbmatrix for the next matrix Pj
        Pbarj = eye(N -(j - 1),N -(j - 1)) - 2*(x_bar + e1(1: N - (j - 1)))*((x_bar + e1(1:N - (j - 1)))') /norm(x_bar + e1(1: N - (j - 1)))^2;
        Pj = zeros(N);
        
        % Position Pjbar below, in the right corner of Pj
        Pj(j : end, j : end) = Pbarj; 
        
        % Fill the remaining elements on the diagonal with 1
        Pj(1:j-1, 1:j-1) = eye(j-1); 
    
        % Construct the new matrix Bj
        Bj = Pj * Bj * Pj;
        
        % Using the property of symmetric matrices: the eigenvectors of
        % distinct eigenvalues are orthogonal (so we are assuming to have
        % distinct eigenvalues)
    
        % Matrix containing the eigenvectors computed previously
        X = eigvects(:, 1 :j - 1);
        
        % Solve the linear system to find the first y elements that along with
        % x_bar, will be used the jth eigenvector of A
        
        
        D =  X' * S(:,1: j - 1);
        c = - X' * S(:,j : end) * x_bar;
        y = D \ c;
   
        % Compute the jth eigenvector and the residual norm
        eigvects(:, j) = S *[y ; x_bar];
        eigvects(:, j) = eigvects(:, j) / norm(eigvects(:, j));
        residualnorms(j) = norm(A * eigvects(:,j) - eigvalues(j) * eigvects(:,j));
    end

end