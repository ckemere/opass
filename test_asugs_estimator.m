

% assume x is already loaded - tetrode data

% what part of the dataset to use
% mT=size(x,1);
mT = 30000 * 180; % 3 minutes
xa=x(1000:mT,:);

[N,numCh]=size(xa);
samplingrate=3e4; % 30 kHz
fprintf('%d s of data\n',N/samplingrate);

%% Set paramters
P=round(4e-3*samplingrate); % window size is 3ms
maxpoint=round(1.7e-3*samplingrate); % where to align waveform peaks
K=3; % Number of PCA components to use.
sig=std(xa); % Noise standard deviation estimate
thres=3*sig; % Detection voltage threshold %% was 3.5

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


CC = zeros(K*numCh,maxtimepoints);
fprintf('Forming lag matrix.');
for i=0:(size(CC,2)-1)
    CC(:,i+1) = A'*reshape(x(i+[1:P],:), P*numCh, 1);
end
fprintf('\n');
%%%%% "sig" will be the covariance matrix of a 90x4 (i.e., 360x1) observation
sig = CC*CC';
sig = sig/size(CC,2);

%% Set parameters:
params.alph=0.05;
params.alph_lamda = 1/params.alph;
params.kappa_0=0.01;
params.nu_0=K*numCh;
params.Phi_0=100*eye(K*numCh);
params.a_pii=1;
params.b_pii=2.1;
% params.b_pii=1e7;
params.bet=1./(30*samplingrate);
params.samplingrate = samplingrate;
params.maxtimepoints = maxtimepoints;
params.verbose = 1;

% Phi_options = [1, 10, 100, 1000]; % Phi_0
% options = [0.001, 0.005, 0.01, 0.1, 0.2, 0.5]; % kappa_0
options = [0.001, 0.01, 0.05, 0.1, 0.2];
for i = 1:length(options)
params.alph=options(i);
% params.Phi_0=Phi_options(i)*eye(K*numCh);
tic;
[z,gam,ngam,muu,Phi,nu,kappa,S]=asugs_estimator(xa,A,sig,params); 
time1 = toc
number_of_classes(i) = sum(ngam > 0);
end

run plot_clusters.m;
