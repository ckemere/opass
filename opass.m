function [z,gam,ngam,muu,lamclus,nu,kappa,Phi,S]=opass(x,A,params);
% Runs the OPASS algorithm
% x is of length N
% A is size PxK
% of spikes in the data.
% params passes in a list of parameters:
% params.alph is the parameter of the CRP
% params.kappa_0, prior precision of mean on NW distribution
% params.nu_0, prior precision of Wishart part of NW distribution
% params.Phi_0, prior cluster covariance*nu_0
% params.a_pii and params.b_pii are the hyperparameters on the probability
% of seeing a spike
%%
apii=params.a_pii;
bpii=params.b_pii;
alph=params.alph;
Phi0=params.Phi_0;
nu_0=params.nu_0;
kappa_0=params.kappa_0;
samplingrate=params.samplingrate;
%% Internal Parameters
Cmax=50;
curndx=0;
lookahead=5*samplingrate/100;
rang=40 * samplingrate/1000;
%%
N=numel(x);
[P,K]=size(A);


%%%%% "sig" will be the covariance matrix of a 90x1 observation
%%%%% "lamda" is just sig inverse
%% Calculate precision matrix
[acf] = xcorr(x,x,1, 'coeff');
%acf(3)=0;
if abs(acf(3))<1e-3
    acf(3)=0;
end
%acf(3)=0;
lambi=zeros(P);
for p=1:P
    lambi(p,:)=1-p:P-p;
end
sig=acf(3).^abs(lambi)*cov(x(1:10*samplingrate));
sig(1:(P+1):P^2)=cov(x(1:10*samplingrate));
lamda=inv(sig);
detlamb=det(lamda);




%%
thr=log(apii/(bpii-apii));

nu=repmat(nu_0,Cmax,1);
Phi=cell(Cmax,1);
lamclus=cell(Cmax,1);
for c=1:Cmax
    Phi{c}=Phi0;
    lamclus{c}=inv(Phi0)*nu_0;
end
muu0=zeros(K,1);
kappa=kappa_0*ones(Cmax,1);
muu=zeros(K,Cmax);
%%
xpad=[x;zeros(P,1)];
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
lamclusS=cell(Cmax,1);
mT=N;
sz=0;
%%



while curndx<N-P-rang
    %% set up parameters
%     pii=(apii+sz)./(bpii+curndx);
%     thr=log(pii./(1-pii));
    ndx=(curndx+1:min(mT-P-rang,curndx+lookahead));n=numel(ndx);
    ndxwind=bsxfun(@plus,ndx,[0:P-1]');

    lpi_c=log(ngam./(alph+nz));
    lpi_c(C+1)=log(alph./(alph+nz));

    xwind=xpad(ndxwind);
    %% calc llk
    lnone=0.5*log(detlamb)-.5*sum((xwind.*((lamda)*xwind)));
    lon=zeros(C+1,n);
    for c=1:C+1
        % Calculate the log likelihood of a spike from cluster c being in the data
        %         lon(c,:)=getllk(xwind,muu(:,c),A,lamclus,sig,kappa(c));
        % (includes C+1 as a special case)
        Q=sig+(1+kappa(c))./kappa(c)*A*(lamclus{c}\A');
        xwindm=bsxfun(@minus,xwind,A*muu(:,c));
        if (c < C+1)
            Re=(ndx-tlastspike(c)) < 5 * samplingrate/1000; % refractory period
            lon(c,:)=-sum(log(diag(chol(Q))))-.5*sum(xwindm.*(Q\xwindm))-double(Re)*1e5;
        else
            lon(c,:)=-sum(log(diag(chol(Q))))-.5*sum(xwindm.*(Q\xwindm))-double(Re)*1e5;
        end
    end

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
    [~,Cnew]=max(lon(:,Q));
    if Cnew>C
        C=Cnew;
        if (C == 50)
            keyboard
        end
    end
    tlastspike(Cnew)=Qt;


    Qmat=A'*lamda*A+lamclus{Cnew};
    yhat=Qmat\(A'*lamda*xwind(:,Q)+lamclus{Cnew}*muu(:,Cnew));
    %     yhat=A'*xwind(:,Q);
    gam(Qt)=Cnew;
    ngam(Cnew)=ngam(Cnew)+1;

    muuold=muu(:,Cnew);
    muu(:,Cnew)=(yhat + muu(:,Cnew)*kappa(Cnew))./(kappa(Cnew)+1);
    dmuu=muuold-muu(:,Cnew);

    Phi{Cnew}=Phi{Cnew}+(kappa(Cnew)+1)./(kappa(Cnew)+1)*(yhat-muu(:,Cnew))*(yhat-muu(:,Cnew))'+inv(Qmat)+ngam(Cnew)*dmuu*dmuu';
    kappa(Cnew)=kappa(Cnew)+1;
    nu(Cnew)=nu(Cnew)+1;
    lamclus{Cnew}=inv(Phi{Cnew})*nu(Cnew);
    %     S(:,Qt)=Bci*(A'*lamb*xwind(:,Q)+lamclus*muu(:,gam(Qt)));
    
    muuS{Cnew}=[muuS{Cnew},muu(:,Cnew)];
    lamclusS{Cnew}{size(muuS{Cnew},2)}=lamclus{Cnew};
    xpad(Qt:Qt+P-1)=xpad(Qt:Qt+P-1)-A*yhat;

    S(:,Qt)=yhat;
    sz=sz+1;
    curndx=Qt+1;
    %     continue
    % keyboard
end

