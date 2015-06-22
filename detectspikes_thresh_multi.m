function [timepoints,spikes]=detectspikes_thresh_multi(x,thres,samplingrate,P,maxpoint)
% function [timepoints,spikes]=detectspikes_thresh_multi(x,thres,samplingrate,P,maxpoint)
% x is size NxK
% thres is a vector of size 1XK
% P is the window size
% maxpoint is the offset from the maximum point
% we assume spikes go upward!

timepoints=[];
t=2*P;
[N,K]=size(x);
while t<N-P
    wind=x(t-P:t+P,:);
    [val,ndx]=max(wind);
    ch = find(val>thres,1);
    if ~isempty(ch)
        if ndx(ch)<P*3/2;
            timepoints=[timepoints,t-P-1+ndx];
            t=t+round(4e-3*samplingrate);
        end
    end
    t=t+P;
end
% get array of detected spikes
spikes=zeros(P*K,numel(timepoints));
for t=1:numel(timepoints)
    for k = 0:(K-1)
        spikes([1:P] +k*P,t)=x(timepoints(t)+(-maxpoint+1:P-maxpoint),k+1);
    end
end
