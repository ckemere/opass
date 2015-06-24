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
Cmax=75;
curndx=0;
lookahead=10*samplingrate/1000;
rang=3 * samplingrate/1000;
%%
[N,D]=size(x);
[PD,K] = size(A);
P = PD/D;

CC = zeros(K,maxtimepoints);
fprintf('Forming lag matrix.');
for i=0:(size(CC,2)-1)
    CC(:,i+1) = A'*reshape(x(i+[1:P],:), PD, 1);
end
fprintf('\n');

%%%%% "sig" will be the covariance matrix of a 90x4 (i.e., 360x1) observation
%%%%% "lamda" is just sig inverse
sig = CC*CC';
sig = sig/size(CC,2);
logDetSig = 2*sum(log(diag(chol(sig))));
lamda=inv(sig);
% detlamb=det(lamda);
%%
thr=log(apii/(bpii-apii));

Nu=repmat(Nu_0,Cmax,1);
Phi=cell(Cmax,1);
for c=1:Cmax
    Phi{c}=Phi0;
end
muu0=zeros(K,1);
Kappa=Kappa_0*ones(Cmax,1);
muu=zeros(K,Cmax);
%%
xpad=[x;zeros(PD,D)];
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

% generate fractional shift matrices.
%  these are PDxPD matrices which circularly shift the data
%  between -0.4 and + 0.5 samples
shift_fractions = [-0.4:0.1:0.5];
for i = 1:length(shift_fractions)

end

while curndx<N-P-rang
    %% set up parameters
%     pii=(apii+sz)./(bpii+curndx);
%     thr=log(pii./(1-pii));
    ndx=(curndx+1:min(mT-P-rang,curndx+lookahead));n=numel(ndx);
    ndxwind=bsxfun(@plus,ndx,[0:P-1]');

    xwind = zeros(K,n);
    for i = 1:n
        xwind(:,i) = A'*reshape(xpad(ndxwind(:,i),:),PD,1);
    end

    %% calc llk
    lnone=-P/2*log(2*pi) - 0.5*logDetSig-.5*sum((xwind.*((lamda)*xwind)));
    lon=zeros(C+1,n);
    for c=1:C+1
        % Calculate the log likelihood of a spike from cluster c being in the data
        %         lon(c,:)=getllk(xwind,muu(:,c),A,lamclus,sig,Kappa(c));
        % (includes C+1 as a special case)
        
        % Q = covariance of observation with neuron prior
        %   = sig + r*A*Phi{c}*A'
        % lamda = inv(sig) % Phi{c}=inv(lamclus{c})
        r = Kappa(c)/(1+Kappa(c));
        Q = sig + r*Phi{c};
        % Qupdate = r*inv(Phi{c}) + lamda;
        % % cholQupdate = chol(Qupdate);
        % Qinv = lamda - lamda * inv(cholQupdate)*inv(cholQupdate)' * lamda;
        % logDetQ = 2*log(diag(cholQupdate)) + K*log(r) + log(det(Phi{c})) + logDetSig ;

        xwindm=bsxfun(@minus,xwind,muu(:,c));
        if (c < C+1)
            Re=(ndx-tlastspike(c)) < 5 * samplingrate/1000; % refractory period
        else
            Re = 0;
        end

        lon(c,:)=-P/2*log(2*pi)-sum(log(diag(chol(Q))))-.5*sum(xwindm.*(Q\xwindm))-double(Re)*1e5;
        % lon(c,:) = -P/2*log(2*pi) - 0.5 * logDetQ  - 0.5*sum(xwindm .* (Qinv * xwindm)) + ...
            % -double(Re)*1e5;
    end

    lpi_c=log(ngam./(alph+nz));
    lpi_c(C+1)=log(alph./(alph+nz));
    lon=bsxfun(@plus,lpi_c(1:C+1,:),lon);

    % Sum over all neurons
    if (C > 0)
        H=bsxfun(@minus,lon,max(lon));
        Hadj=log(sum(exp(H)));
        lthr=lnone-max(lon)-Hadj;
    else
        lthr = lnone-lon;
    end
    %% Find new spike
    Q=find(lthr<thr,1,'first');
    % no spike
    if (numel(Q)==0) || Q>lookahead-rang
        curndx=curndx+lookahead-rang;
        continue
    end


    % Spike detected!!
    [peak_thr,offset]=min(lthr(Q:min(Q+rang,numel(lthr))));
    Q=Q+offset-1; % this is the "peak"
    nz=nz+1;
    Qt=Q+curndx;
    z(Qt)=1;

    [~,cTemp]=max(lon(:,Q));
    yhat=xwind(:,Q);

    % if (C==6)
        % keyboard
    % end


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
        fprintf('New class (%d). Thr: %d\n', C, peak_thr);
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
    % muuS{cSpike}=[muuS{cSpike},muu(:,cSpike)];
    % hiS{cSpike}{size(muuS{cSpike},2)}=PhiS{cSpike};


    gam(Qt)=cSpike;
    ngam(cSpike)=ngam(cSpike)+1;

    % for d = 1:D
        % xpad(Qt:Qt+P-1,d)=xpad(Qt:Qt+P-1,d)-A([1:P] + (d-1)*P,:)*yhat;
    % end

    S(:,Qt)=yhat;
    sz=sz+1;
    curndx=Qt+samplingrate/1000; % step forward 1 ms

    if (mod(curndx,100000) == 0)
        fprintf('*');
    end

    % adapt alpha
    alph = C / (alph_lamda + log(sz));
    %     continue

end

