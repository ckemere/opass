function [z,gam,ngam,muu,lamclus,Nu,Kappa,Phi,S]=opass(x,A,params);
% Runs the OPASS algorithm
% x is of length N
% A is size PxK
% of spikes in the data.
% params passes in a list of parameters:
% params.alph is the parameter of the CRP
% params.Kappa_0, prior precision of mean on NW distribution
% params.Nu_0, prior precision of Wishart part of NW distribution
% params.Phi_0, prior cluster covariance*Nu_0
% params.a_pii and params.b_pii are the hyperparameters on the probability
% of seeing a spike
%%
apii=params.a_pii;
bpii=params.b_pii;
alph=params.alph;
Phi0=params.Phi_0;
Nu_0=params.Nu_0;
Kappa_0=params.Kappa_0;
%% Internal Parameters
Cmax=50;
curndx=0;
lookahead=500*3;
rang=40*3;
%%
N=numel(x);
[P,K]=size(A);

%% Calculate precision matrix
%%% This creates a PxP matrix that has the data variance on the diagonal,
%%%  and which falls off along the off diagonal as the lag-1 autocorrelation.
%%%  Why not calculate the PxP autocorrelation matrix, I don't know.
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
sig=acf(3).^abs(lambi)*cov(x(1:1e5));
sig(1:(P+1):P^2)=cov(x(1:1e5));
lamda=inv(sig);
detlamb=det(lamda);

%%
pii=apii./bpii;
Nu=repmat(Nu_0,Cmax,1);
Phi=cell(Cmax,1);
lamclus=cell(Cmax,1);
for c=1:Cmax
    Phi{c}=Phi0;
    % lamclus{c}=inv(Phi0)*Nu_0;
    lamclus{c}=inv(Phi0) / (2*Nu_0);
end
muu0=zeros(K,1);
Kappa=Kappa_0*ones(Cmax,1);
muu=zeros(K,Cmax);
%%
xpad=[x;zeros(P,1)];
%%
C=0;
nz=0;
z=zeros(N,1);
gam=zeros(N,1);
piip=zeros(N,1);
lthet=zeros(Cmax,1);
ngam=zeros(Cmax,1);
S=zeros(K,N);
lpii=log(pii);
lnpii=log(1-pii);
thr=log(pii./(1-pii));
xm=x;
tlastspike=zeros(Cmax,1);
muuS=cell(Cmax,1);
lamclusS=cell(Cmax,1);
mT=N;
%%
while curndx<N-P-rang
    ndx=(curndx+1:min(mT-P-rang,curndx+lookahead));n=numel(ndx);
    ndxwind=bsxfun(@plus,ndx,[0:P-1]');
    lthet=log(ngam./(alph+nz));
    lthet(C+1)=log(alph./(alph+nz));
    xwind=xpad(ndxwind);
    %% calc llk
    lnone=-P/2*log(2*pi)+.5*log(detlamb)-.5*sum((xwind.*((lamda)*xwind)));  % assuming Gaussian background noise
    lon=zeros(C+1,n);
    for c=1:C+1 % (assuming proper initialization for Kappa and lamclus for C+1
        %         lon(c,:)=getllk(xwind,muu(:,c),A,lamclus,sig,Kappa(c));
        Q=sig+(1+Kappa(c))./Kappa(c)*A*(lamclus{c}\A');
        xwindm=bsxfun(@minus,xwind,A*muu(:,c));
        Re=(ndx-tlastspike(c))<5 * 30000/1000;
        lon(c,:)=-P/2*log(2*pi)-sum(log(diag(chol(Q))))-.5*sum(xwindm.*(Q\xwindm))-double(Re)*1e5;
    end

    lon=bsxfun(@plus,lthet(1:C+1,:),lon);
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
    % new spike
    [~,offset]=min(lthr(Q:min(Q+rang,numel(lthr))));
    Q=Q+offset-1;
    nz=nz+1;
    Qt=Q+curndx;
    z(Qt)=1;
    [~,Cspike]=max(lon(:,Q));
    if Cspike>C
        C=Cspike; % Set up a new class!
    end
    tlastspike(Cspike)=Qt;
    Qmat=A'*lamda*A+lamclus{Cspike};
    yhat2=Qmat\(A'*lamda*xwind(:,Q)+lamclus{Cspike}*muu(:,Cspike));
    yhat = A'*xwind(:,Q);
    % xhat=A*yhat ~= xwind(:,Q);

    gam(Qt)=Cspike; % assigned class for time Qt
    ngam(Cspike)=ngam(Cspike)+1; % number of spikes assigned to class Cspike


    Phi{Cspike} = 2*Nu(Cspike)/(1+2*Nu(Cspike))*Phi{Cspike} + 1/(1+2*Nu(Cspike)) * (Kappa(Cspike)/(1+Kappa(Cspike))) * (yhat - muu(:,Cspike)) * (yhat-muu(:,Cspike))';

    % Update mean for class Cspike
    muuold=muu(:,Cspike);
    muu(:,Cspike)=(yhat + muu(:,Cspike)*Kappa(Cspike))./(Kappa(Cspike)+1);
    dmuu=muuold-muu(:,Cspike);

    % Update Kappa 
    Kappa(Cspike)=Kappa(Cspike) + 1;

    Nu(Cspike)=Nu(Cspike) + 1/2; % I think this should be 1/2??? (rather than 1)
    lamclus{Cspike}=inv(Phi{Cspike}) / (2*Nu(Cspike));
    
    curndx=Qt+1;
    muuS{Cspike}=[muuS{Cspike},muu(:,Cspike)]; % tracks the evolution of the mean spike waveform
    lamclusS{Cspike}{size(muuS{Cspike},2)}=lamclus{Cspike}; % tracks the evolution of the variability
    xpad(Qt:Qt+P-1)=xpad(Qt:Qt+P-1)-A*yhat;
    %     continue

    keyboard
end


















