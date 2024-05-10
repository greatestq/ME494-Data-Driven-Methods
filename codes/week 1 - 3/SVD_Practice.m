%% Example 1
X = rand(5,3);
[U, S, V] = svd(X)
[Uhat, Shat, V] = svd(X,'econ')
%Question: How do I prove that U^orth spans the orthogonal space

%% Example 2
A = imread('KDB.jpg');
X = double(rgb2gray(A));
nx = size(X,1); ny = size(X,2);
imagesc(X), axis off, colormap gray;

[U, S, V] = svd(X);
for r = [10 20 30 50 600]
    Xapprox = U(:,1:r)*S(1:r, 1:r)*V(:,1:r)';
    figure, imagesc(Xapprox), axis off, colormap gray;
    title(['r=', num2str(r,'%d'),'']);
end

%subplot(1,2,1), semilogy(diag(S),'k')
%subplot(1,2,2), plot(cumsum(diag(S))/sum(diag(S)),'k')
