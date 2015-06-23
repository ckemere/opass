function [z,gam,ngam,muu,Phi,Nu,Kappa,S]=asugs(x,A,params);
% Runs the modified OPASS algorithm (using the ASUGS updates)
% x is of size NxK
% contains the data including spikes
% params passes in a list of parameters:
% params.alph is the parameter of the CRP
% params.kappa_0, prior precision of mean on NW distribution
% params.nu_0, prior precision of Wishart part of NW distribution
% params.Phi_0, prior cluster covariance*nu_0
% params.a_pii and params.b_pii are the hyperparameters on the probability
% of seeing a spike
% params.samplingrate is the sampling rate
% params.maxtimepoints is the size of the window for covariance matrix estimation
%%
apii=params.a_pii;
bpii=params.b_pii;
alph=params.alph;
alph_lamda=params.alph_lamda;
Phi0=params.Phi_0;
Nu_0=params.nu_0;
Kappa_0=params.kappa_0;
samplingrate=params.samplingrate;
maxtimepoints=params.maxtimepoints
%% Internal Parameters
Cmax=300;
curndx=0;
lookahead=10*samplingrate/1000;
rang=3 * samplingrate/1000;
%%
[N,D]=size(x); % N is number of time points, D is number of dimensions
if ~iscell(A)
    if D > 1
        error('The matrix A should be a cell array of length size(x,2).');
    end
    Atemp{1} = A
    A = Atemp;
    clear Atemp;
end

[P,K] = size(A{1});

for d = 1:D
    CC = zeros(P,maxtimepoints);
    fprintf('Forming lag matrix.');
    for i=0:(size(CC,2)-1)
        CC(:,i+1) = x(i+[1:P],d);
    end
    %%%%% "sig" will be the covariance matrix of a 90x4 (i.e., 360x1) observation
    %%%%% "lamda" is just sig inverse
    sig{d} = CC*CC' / size(CC,2); 
    logDetSig{d} = 2*sum(log(diag(chol(sig{d}))));
    lamda{d}=inv(sig{d});
end
fprintf('\n');


%%
thr=log(apii/(bpii-apii));

Nu=repmat(Nu_0,Cmax,1);
Phi=cell(Cmax,1);
for c=1:Cmax
    Phi{c}=Phi0;
end
muu0=zeros(D*K,1);
Kappa=Kappa_0*ones(Cmax,1);
muu=zeros(K*D,Cmax);
%%
xpad=[x;zeros(P,D)];
%%
C=0;
nz=0;
z=zeros(N,1);
gam=zeros(N,1);
lpi_c=zeros(Cmax,1);
ngam=zeros(Cmax,1);
S=zeros(K,N);
tlastspike=zeros(Cmax,1);
muuS=cell(Cmax,1);
PhiS=cell(Cmax,1);
mT=N;
sz=0;
%%

while curndx<N-P-rang
    if (curndx > 1*samplingrate)
        return;
    end
    %% set up parameters
%     pii=(apii+sz)./(bpii+curndx);
%     thr=log(pii./(1-pii));
    ndx=(curndx+1:min(mT-P-rang,curndx+lookahead));n=numel(ndx);
    ndxwind=bsxfun(@plus,ndx,[0:P-1]');

    xwind = zeros(P,n,D);
    for d = 1:D
        for i = 1:n
            xwind(:,i,d) = reshape(xpad(ndxwind(:,i),d),P,1);
        end
    end

    %% calc llk
    lnone = 0;
    for d = 1:D
        lnone = lnone + -P/2*log(2*pi) - 0.5*logDetSig{d} + ...
            -.5*sum((xwind(:,:,d).*((lamda{d})*xwind(:,:,d))));
    end
    lon=zeros(C+1,n);
    for c=1:C+1
        % Calculate the log likelihood of a spike from cluster c being in the data
        %         lon(c,:)=getllk(xwind,muu(:,c),A,lamclus,sig,Kappa(c));
        % (includes C+1 as a special case)
        
        % Q = covariance of observation with neuron prior
        %   = sig + r*A*Phi{c}*A'
        % lamda = inv(sig) % Phi{c}=inv(lamclus{c})
        r = Kappa(c)/(1+Kappa(c));
        for d = 1:D
            phi_d = Phi{c}((d-1)*K + [1:K],(d-1)*K + [1:K]);
            Qupdate = 1/r*inv(phi_d) + A{d}'*lamda{d}*A{d};
            cholQupdate = chol(Qupdate);
            Qinv = lamda{d} - lamda{d} * A{d} * inv(Qupdate) * A{d}' * lamda{d};
            % Qinv = lamda{d} - lamda{d} * A{d} * inv(cholQupdate)*inv(cholQupdate)' * A{d}' * lamda{d};
            logDetQ = 2*log(det(cholQupdate)) + K*log(r) + log(det(phi_d)) + logDetSig{d};

            xwindm=bsxfun(@minus,xwind(:,:,d),A{d}*muu([1:K] + (d-1)*K,c));

            lon(c,:) = lon(c,:) + -P/2*log(2*pi) - 0.5 * logDetQ  - 0.5*sum(xwindm .* (Qinv * xwindm));
        end
        if (c < C+1)
            Re=(ndx-tlastspike(c)) < 5 * samplingrate/1000; % refractory period
            lon(c,:) = lon(c,:) - double(Re)*1e5;
        else
            Re = 0;
        end
    end

    lpi_c=log(ngam./(alph+nz));
    lpi_c(C+1)=log(alph./(alph+nz));
    lon=bsxfun(@plus,lpi_c(1:C+1,:),lon);

    % Sum over all neurons
    H=bsxfun(@minus,lon,max(lon));
    Hadj=log(sum(exp(H)));
    lthr=lnone-max(lon)-Hadj;
    %% Find new spike
    Q=find(lthr<thr,1,'first');
    % no spike
    if (numel(Q)==0) || Q>lookahead-rang
        curndx=curndx+lookahead-rang;
        continue
    end


    % Spike detected!!
    [~,offset]=min(lthr(Q:min(Q+rang,numel(lthr))));
    Q=Q+offset-1; % this is the "peak"
    nz=nz+1;
    Qt=Q+curndx;
    z(Qt)=1;

    [~,cTemp]=max(lon(:,Q));
    for d = 1:D
        yhat([1:K] + (d-1)*K,1) = A{d}'*xwind(:,Q,d);
    end
    for c = 1:C+1
        r = Kappa(c)/(1+Kappa(c));
        dmuu = yhat - muu(:,c);

        llClass(c) = lpi_c(c) + ...
           (K/2)*log(r/Nu(c)) + gammaln((Nu(c)+1)/2) - gammaln((Nu(c)+1 - K)/2) + ...
           -(1/2)*log(det(Phi{c})) + (-(Nu(c)+1)/2)*log(1 + r/Nu(c) * sum(dmuu.* (Phi{c} \ dmuu)));

        % evid(c) = -(1/2)*log(det(Phi{c})) + (-(Nu(c)+1)/2)*log(1 + r/Nu(c) * sum(dmuu.* (Phi{c} \ dmuu)));
        % kappa_gamma(c) = (K/2)*log(r/Nu(c)) + gammaln((Nu(c)+1)/2) - gammaln((Nu(c)+1 - K)/2);
    end
    [~,cSpike] = max(llClass);


    if cSpike>C
        C=cSpike;
        if (C == Cmax)
            keyboard
        end
    end
    tlastspike(cSpike)=Qt;

    %%% Update class statistics.
    % Phi is the inverse of the mean of the Wishart distribution (~ covariance of mean of cluster)
    Phi{cSpike} = Nu(cSpike)/(1+Nu(cSpike)) * Phi{cSpike} + ...
       1 / (1 + Nu(cSpike)) * Kappa(cSpike) / (1 + Kappa(cSpike)) * (yhat - muu(:,cSpike))*(yhat - muu(:,cSpike))';
    muu(:,cSpike) = (yhat + muu(:,cSpike)*Kappa(cSpike))./(Kappa(cSpike)+1);
    Kappa(cSpike) = Kappa(cSpike) + 1;
    Nu(cSpike) = Nu(cSpike) + 1;

    % Lamda is going to be the precision (inverse cov) of the cluster
    % lamclus{cSpike}=inv(Phi{cSpike});
    muuS{cSpike}=[muuS{cSpike},muu(:,cSpike)];
    PhiS{cSpike}{size(muuS{cSpike},2)}=PhiS{cSpike};


    gam(Qt)=cSpike;
    ngam(cSpike)=ngam(cSpike)+1;

    keyboard
    for d = 1:D
        xpad(Qt:Qt+P-1,d)=xpad(Qt:Qt+P-1,d)-A{d}*yhat([1:K] + (d-1)*K);
    end

    S(:,Qt)=yhat;
    sz=sz+1;
    curndx=Qt+samplingrate/5000; % step forward 0.2 ms

    % adapt alpha
    alph = C / (alph_lamda + log(sz));
    %     continue

end

