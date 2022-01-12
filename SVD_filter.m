function Signal=SVD_filter(Data)
% SVD filter
[U,~,latent] = pca(Data);
Qlatent=latent(latent>latent(end)*10); % thresholding
QU=U(:,1:size(Qlatent,1));
P=QU*QU';
Signal=Data*P;
end