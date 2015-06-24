% what part of the dataset to use
mT=size(x,1);
xa=x(1000:mT,:);

[N,numCh]=size(xa);
samplingrate=3e4; % 30 kHz
fprintf('%d s of data\n',N/samplingrate);

%% Set paramters
P=round(3e-3*samplingrate); % window size is 3ms
maxpoint=round(1.2e-3*samplingrate); % where to align waveform peaks
K=3; % Number of PCA components to use.
sig=std(xa); % Noise standard deviation estimate
thres=3.5*sig; % Detection voltage threshold

%% Detect spike waveforms
maxtimepoints=120*samplingrate; % limit to first 120s to simulate online system
[timepoints,spikes]=detectspikes_thresh_multi(xa(1:maxtimepoints,:),thres,samplingrate,P,maxpoint);

%% Reduce dimensionality
% [A,S,V]=svds(spikes,K);
A = zeros(P*numCh, K*numCh);
% A = cell(numCh,1);
for d = 1:numCh
    [u,s,v] = svds(spikes([1:P] + (d-1)*P,:),K);
    % A{d} = u;
    A([1:P]+(d-1)*P, [1:K] + (d-1)*K) = u;
end

%% Set parameters:
params.alph=0.001;
params.alph_lamda = 1/params.alph;
params.kappa_0=0.2;
params.nu_0=K*numCh;
params.Phi_0=eye(K*numCh);
params.a_pii=1;
params.b_pii=2.01;
% params.b_pii=1e7;
params.bet=1./(30*samplingrate);
params.samplingrate = samplingrate;
params.maxtimepoints = maxtimepoints;

tic;
[z,gam,ngam,muu,Phi,nu,kappa,S]=asugs_m(xa,A,params); 
time1 = toc

%% Plot non-trivial clusters
% Plot spikes
C=max(gam);
col=hsv(C);
figure(1);clf;hold on
for c=1:C
plot(S(1,gam==c),S(2,gam==c),'.','Color',col(c,:),'markersize',20)
end
hold off
xlabel('PCA Component 1','FontSize',16)
ylabel('PCA Component 2','FontSize',16)
title('Inferred y_k Values for Detected Spikes','FontSize',18)
