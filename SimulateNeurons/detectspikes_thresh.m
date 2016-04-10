function [timepoints,spikes]=detectspikes_thresh(x,thres,samplingrate,P,maxpoint)
% x - NxnumCh matrix were N is number of observations, numCh is number of
%     channels
% thres - threshold level. When there is a threshold crossing in ANY
%         channel, a spike is extracted. 
% samplingrate - sampling rate in Hz
% P - 2P is the length of the spike snippet
% maxpoint - the point where the maximum of the spike waveform will be
%            located
timepoints=[];
t=2*P;
N=size(x, 1);
numCh = size(x, 2);
while t<N-P
    wind=x(t-P:t+P, :);
    [val,ndx]=max(wind);
    [globalmaxval, chWithMaxVal] = max(val);
    globalndx = ndx(chWithMaxVal);
    if globalmaxval>thres
        if globalndx<P*3/2;
            timepoints=[timepoints,t-P-1+globalndx];
            t=t+round(4e-3*samplingrate);
            
        end
    end
    t=t+P;
end
% get array of detected spikes
spikes=zeros(P,numel(timepoints), numCh);
for n = 1:numCh
    for t=1:numel(timepoints)
        spikes(:,t,n)=x(timepoints(t)+(-maxpoint+1:P-maxpoint), n);
    end
end