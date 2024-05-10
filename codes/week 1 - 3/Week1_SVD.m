A = [1 0 2; 
     1 1 2; 
     0 0 0]; % 3x3 rank 2 matrix

[U, S, V] = svd(A);

B = [1 0 1; 1 0 -1];
[U, S, V] = svd(B)

[u, v] = eig(transpose(B)*B)
%To check if first two columns of U are lin indep.
transpose(U(:,1))*U(:,2);
%To check if first two columns of V are lin indep.
transpose(V(:,1))*V(:,2);
%To check that R(A) and N(A*) are orthogonal
transpose(U(:,1))*U(:,3);
transpose(U(:,2))*U(:,3);
%To check that R(A*) and N(A) are orthogonal
transpose(V(:,1))*V(:,3);
transpose(V(:,2))*V(:,3);