function [eigvalue, eigvect] = invpower_method(A, v, maxIter, relTol)
    
    % Function that computes the smallest eigenvalue and corresponding
    % eigenvector of the given matrix
    % Inputs:
    % A: Symmetric matrix (could also not be symm)
    % v: initial vector
    % maxIter: max number of iterations
    % relTol: relative tolerance for the stopping criteria
    % Outputs:
    % eigvalue, eigenvect: smallest eigenvalue and corresponding eigenvector of the given matrix

    % Size of v
    v_size = size(A, 1); 
    
    % Inizialize values for the method

    % Will contain the approximation of the eigenvalue computed ad each step
    lambdas = zeros(maxIter, 1);
    % Inizialize at a very large number
    lambdas(1) = 1e8; 

    % Will contain the approximation of the eigenvector computed ad each step
    v_vectors = zeros(v_size, maxIter + 1);

    % Normalize v
    v_vectors(:, 1) = v / norm(v);
    
    % Start the loop
    for k = 1:maxIter
        % Solve the linear system to find v_tilde = A^-1*v_vectors(:, k)
        v_tilde = A \ v_vectors(:, k);

        % Compute the new approximation of the eigenvalue
        lambdas(k + 1) = dot(v_vectors(:, k), v_tilde);
        
        % Normalize the new approximation of the eigenvector
        v_vectors(:, k + 1) = v_tilde / norm(v_tilde);
        
        % Check if the stopping criteria is met
        if abs(lambdas(k + 1) - lambdas(k)) < relTol * abs(lambdas(k + 1))
            break;
        end
    end
    
    % Save the computed values and vector
    eigvalue = 1 / lambdas(k + 1)  ;
    eigvect = v_vectors(:, k + 1);

end