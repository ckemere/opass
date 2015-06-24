function [timepoints,spikes]=detectspikes_thresh_multi(x,thres,samplingrate,P,maxpoint)
% function [timepoints,spikes]=detectspikes_thresh_multi(x,thres,samplingrate,P,maxpoint)
% x is size NxK
% thres is a vector of size 1XK
% P is the window size
% maxpoint is the offset from the maximum point
% we assume spikes go upward!

tRefract = 2e-3 * samplingrate;
timepoints=[];
t=2*P;
[N,D]=size(x);
below_thresh = prod(bsxfun(@lt, x, thres),2);
thresh_excursion_times = find(diff(below_thresh) == -1);
thresh_return_times = find(diff(below_thresh) == 1);

if (thresh_excursion_times(1) < thresh_return_times(1))
    thresh_excursion_times(1) = [];
end

L = min(length(thresh_excursion_times), length(thresh_return_times));
thresh_excursion_times = thresh_excursion_times(1:L);
thresh_return_times = thresh_return_times(1:L);
previous_refractory_periods = thresh_excursion_times - ...
    thresh_return_times;

timepoints=thresh_excursion_times(previous_refractory_periods > tRefract);

peak_window = [-0.1e-3*samplingrate : 0.2e-3 * samplingrate];
L = length(peak_window);

% get array of detected spikes
spikes=zeros(P*D,numel(timepoints));
for t=1:numel(timepoints)
    % we are going to align to the center of mass in a window around the peak
    xx = max(x(timepoints(t)+peak_window,:),0);
    com = round(sum(abs([1:L]*xx))/sum(abs(ones(1,L)*xx)));
    offset = timepoints(t) + com - peak_window(1) - maxpoint;

    for d = 0:(D-1)
        spikes([1:P] +d*P,t) = x(offset + [1:P],d+1);
    end
end

