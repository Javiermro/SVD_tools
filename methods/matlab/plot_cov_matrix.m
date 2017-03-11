function plot_cov_matrix(phi, W, i)

if(W==0)
    W=eye(size(phi,1));
end
%normphi = sqrt(sum((phi).*(phi)));
%phiC = phi./repmat(normphi,size(phi,1),1);
covMatrix = abs(phi'*W*phi) ;
figure(i); imagesc(covMatrix); colorbar;

end