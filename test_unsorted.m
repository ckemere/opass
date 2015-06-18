%% test all algorithms
% set which algorithms to run
run_FAKEOPASS=false;
run_OPASS=true;
run_OPASS_A=false;
run_M_OPASS=false;
run_M_OPASS_A=false;

% load example data set


% what part of the dataset to use
mT=size(X,1);
x=X(1000:mT,1); % load channel 1
xa=X;

[N,numCh]=size(x);

samplingrate=3e4; % 30 kHz

%% Set paramters
P=round(3e-3*samplingrate); % window size is 3ms
maxpoint=round(1.5e-3*samplingrate); % where to align waveform peaks
K=5; % Number of PCA components to use.
sig=std(x); % Noise standard deviation estimate
thres=3*sig; % Detection voltage threshold

%% Detect spike waveforms
[timepoints,spikes]=detectspikes_thresh(-x,thres,samplingrate,P,maxpoint);
%% Reduce dimensionality
maxtimepoints=120*samplingrate; % limit to first 120s to simulate online system
[U,S,V]=svd(spikes(:,timepoints<maxtimepoints),'econ');
A=U(:,1:K);

%% Set parameters:
params.alph=0.1;
params.kappa_0=0.01;
params.nu_0=0.1;
params.Phi_0=0.1*eye(K);
params.a_pii=1;
params.b_pii=1e5;
params.bet=1./(30*samplingrate);
params.samplingrate = samplingrate;

%% Run fake_opass
if run_FAKEOPASS
    y=A'*spikes; % convert all spikes to low-dimensional feature space
    % normalize
    y=bsxfun(@minus,y,mean(y,2));
    y=bsxfun(@rdivide,y,std(y')');

    %%
    [gam,ngam,muu,Lam,nu,kappa,Phi]=fake_opass(y,params);    
    %% Plot non-trivial clusters;
    C=sum(ngam>10);
    [~,rendx]=sort(ngam,'descend');
    colors=hsv(C);
    figure(1);clf; hold on
    set(0,'defaulttextinterpreter','latex')
    for c=1:C
        plot3(y(1,gam==rendx(c)),y(2,gam==rendx(c)),y(3,gam==rendx(c)),'.','MarkerSize',15,'color',colors(c,:))
    end
    xlabel('pc-1');ylabel('pc-2');zlabel('pc-3');title('\tt FAKE-OPASS','FontSize',20)
    hold off
end
%% run OPASS
if run_OPASS
%     clear params
tic;
    [z,gam,ngam,muu,lamclus,nu,kappa,Phi,S]=opass(x,A,params); time1 = toc;
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
    snames=cell(C+1,1);
    for c=1:C
        snames{c}=num2str(c);
    end
    a=legend(snames);
end
%% run OPASS-A
if run_OPASS_A
%     clear params
    tic;[z,gam,ngam,muu,lamclus,nu,kappa,Phi,S]=opass_a(x,A,params);time2 = toc;
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
    snames=cell(C+1,1);
    for c=1:C
        snames{c}=num2str(c);
    end
    a=legend(snames);
    
end
%%
%% run M_OPASS
if run_M_OPASS
%     clear params
    [z,gam,ngam,muu,lamclus,nu,kappa,Phi,S]=m_opass(xa,A,params);
    %% Plot non-trivial clusters
    % Plot spikes
    C=max(gam);
    col=hsv(C);
    figure(1);clf;hold on
    for c=1:C
        plot(squeeze(S(1,1,gam==c)),squeeze(S(2,1,gam==c)),'.','Color',col(c,:),'markersize',20)
    end
    hold off
    xlabel('PCA Component 1','FontSize',16)
    ylabel('PCA Component 2','FontSize',16)
    title('Inferred y_k Values for Detected Spikes','FontSize',18)
    snames=cell(C+1,1);
    for c=1:C
        snames{c}=num2str(c);
    end
    a=legend(snames);
end

%% run M_OPASS_A
if run_M_OPASS_A
%     clear params
    [z,gam,ngam,muu,lamclus,nu,kappa,Phi,S]=m_opass_a(xa,A,params);
    %% Plot non-trivial clusters
    % Plot spikes
    C=max(gam);
    col=hsv(C);
    figure(1);clf;hold on
    for c=1:C
        plot(squeeze(S(1,1,gam==c)),squeeze(S(2,1,gam==c)),'.','Color',col(c,:),'markersize',20)
    end
    hold off
    xlabel('PCA Component 1','FontSize',16)
    ylabel('PCA Component 2','FontSize',16)
    title('Inferred y_k Values for Detected Spikes','FontSize',18)
    snames=cell(C+1,1);
    for c=1:C
        snames{c}=num2str(c);
    end
    a=legend(snames);
end
