% tic;
% [z,gam,ngam,muu,Phi,nu,kappa,S]=asugs_m(xa(1:450*samplingrate,:),A,params); 
% time1 = toc

max_subplots_per_fig = 4;
num_dim = K*numCh;

numClasses=max(gam);
col=hsv(numClasses);
%% Plot non-trivial clusters

n = 0;
cluster_pairs = zeros(nchoosek(numCh,2),2);
for i = 1:4
    for j = (i+1):4
        n = n + 1;
        cluster_pairs(n,:) = [(i-1)*K+1,(j-1)*K+1];
    end
end

subCount = 0;
for n = 1:length(cluster_pairs)
    subCount = subCount + 1;
    if (subCount > max_subplots_per_fig)
        figure;
        subCount = 1;
    end
    subplot(2,2,subCount)
    hold on
    for c=1:numClasses
        plot(S(cluster_pairs(n,1),gam==c),S(cluster_pairs(n,2),gam==c),'.','Color',col(c,:),'markersize',10)
    end
    hold off
    xlabel('PCA Component 1','FontSize',16)
    ylabel('PCA Component 2','FontSize',16)
    title('Inferred y_k Values for Detected Spikes','FontSize',18)
end
