function W = knn_graph(X, k)

    % Function  that computes the adjacency matrix W based on the knn
    % Inputs:
    % X : Matrix containing the points as rows, their coordinates as columns
    % k : number of neighbors that will be used
    % Output:
    % W: knn adjacency matrix
    
    % Number of points
    N = size(X, 1); 
    
    % Lower part of the matrix, since it is symmetric
    W_lower = (zeros(N));  
    
    % Parameter sigma  for the similarity function
    sigma = 1;
    
    % Similarity function
    similarity = @(xi, xj, sigma) exp(-sum((xi - xj).^2) / (2 * sigma^2));

    % Create row and column indices for the lower triangle (excluding diagonal)
    [row, col] = find(tril(true(N), -1));
    
    % Use arrayfun to apply the similarity function
    W_lower(sub2ind([N, N], row, col)) = arrayfun(@(i, j) similarity(X(i, :), X(j, :), sigma), row, col);

    % The graph is undirected, so the matrix is symmetric
    W = W_lower + W_lower';
    
    % Copy of the matrix W
    W_copy = W; 
    
    % Computing the k-nearest neighborhood adjacency matrix
    for i = 1:N
        % Order the points in descending order wrt the similarity with i
        [~, idx] = sort(W_copy(i, :), 'descend'); 
        % knn = idx(1:k);
    
        % Indices of points that are not in the k-nearest neighborhood
        not_knn = idx(k + 1 : end);
    
        % The entry for these points will be 0
        W(i,not_knn) = 0; 
    end
    
    % Tolerance that will be used to express if an entry of W is 0
    % It is importanto to have the certainty that the matrix will be
    % symmetric
    eps = 1e-20; 
    
    % Symmetrizing method for the matrix W: if j is in the k-nn of i, then also i becomes a
    % neighbor of j and viceversa, so now each row could have more than k values != 0
    
    % Find indices where either W(i,j) or W(j,i) is greater than eps
    idx = (W > eps) | (W' > eps);

    % Update W using logical indexing
    W(idx) = W_copy(idx);
    
    % Check if it is actually symmetric
    dist = norm(W-W');
end