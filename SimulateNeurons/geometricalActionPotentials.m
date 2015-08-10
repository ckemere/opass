
% Parameters
TetrodeSpacing = 15; % microns - space between adjacent tetrode channels
MinDistanceToTetrode = 15; % microns - minimum radial distance between any neuron and a tetrode;
MinDistanceBetweenNeurons = 50; % microns - minimum radial distance between neurons
NeuronBoundsZ = [-25 125]; % microns - height of cubic box in which neurons are generated
NeuronBoundsXY = [-200 200]; % microns - width of cubic box in which neurons are generated
Scaling = 1500; % microvolts / microvolt - scaling factor for excitatory neuron
Alpha = 2; % exponent for distance function (how fast signal level drops off)


% Tetrode is centered at 0, 0
TetrodeLocation = ...
    [-0.5, -0.5, 0;
     -0.5,  0.5, 0;
      0.5,  0.5, 0;
      0.5, -0.5, 0]' * TetrodeSpacing;

load excitatory_waveform;
T = length(excitatory_waveform);
SamplingRate = 30000; % Hz, for this waveform set


% Generate neuron locations
NeuronLocations = [];
for i = 1:10000
    loc(1:2,1) = rand(2,1)* diff(NeuronBoundsXY) + NeuronBoundsXY(1);
    loc(3,1) = rand(1,1)* diff(NeuronBoundsZ) + NeuronBoundsZ(1);
    if (i > 1)
        dd = sqrt(sum(abs(bsxfun(@minus,NeuronLocations, loc))).^2);
        if (min(dd) > MinDistanceBetweenNeurons)
            NeuronLocations = cat(2, NeuronLocations, loc);
        end
    else
        NeuronLocations = loc;
    end
end

N = size(NeuronLocations,2); % number of neurons

td = zeros(N,4); % distance to each tetrode channel
for i = 1:4
    td(:,i) = sqrt(sum(abs(bsxfun(@minus,TetrodeLocation(:,i),NeuronLocations)),1).^2);
end
% sort by distance
[td_min, idx] = sort(min(td,[],2),1,'ascend');
td_sorted = td(idx,:);
% get rid of neurons that are too close to tetrodes
NeuronLocations = NeuronLocations(:,idx);
NeuronLocations(:,td_min < MinDistanceToTetrode) = [];
td_sorted(td_min < MinDistanceToTetrode,:) = [];

N = size(NeuronLocations,2); % number of neurons (modified after dropping)

UpsamplingFactor = 100;

% Simulate
NeuronWaveform = zeros(T,N,4);
UpsampledWaveform = zeros(T * UpsamplingFactor,N,4);
SquishedNeuronWaveforms = zeros(T*4,N);
for d = 1:4
    NeuronWaveform(:,:,d) = excitatory_waveform * transpose(1./(td_sorted(:,d).^Alpha)) * Scaling;
    for n = 1:N
        UpsampledWaveform(:,n,d) = interp(NeuronWaveform(:,n,d), UpsamplingFactor);
    end
end

for i = 1:N
    SquishedNeuronWaveforms(:,i) = reshape(squeeze(NeuronWaveform(:,i,:)),T*4,1);
end
Peaks = squeeze(max(NeuronWaveform,[],1));

figure
peakPlot(Peaks)

fprintf('%d neurons with peak peak more than 1000.\n', sum(max(Peaks,[],2)>1000));
fprintf('%d neurons with peak peak more than 100.\n', sum(max(Peaks,[],2)>100));
fprintf('%d neurons with peak peak more than 40.\n', sum(max(Peaks,[],2)>40));

% Now simulate neurons in time
% Assumptions (see Cheng and Frank 2007)
% 
% 
SimTime = 20 * 60; % 20 minutes
[SpikeTimes, NeuronProps] = simulateSpikeTimes(N, SimTime, @simulatePosition);

TD = SimTime * SamplingRate;
SimData = zeros(TD, 4);

TotalSpikes = 0;
for n = 1:N
    TotalSpikes = TotalSpikes + length(SpikeTimes{n});
end

JitteredPeaks = nan(TotalSpikes,4);
ClusterIDs = nan(TotalSpikes,1);
k = 1;

for n = 1:N
    SpikeIndex{n} = nan(length(SpikeTimes{n}),1);
    ti = fix(SpikeTimes{n} * SamplingRate);
    toff = ceil((SpikeTimes{n} * SamplingRate - ti) * UpsamplingFactor);
    for i = 1:length(SpikeTimes{n})
        SpikeIndex{n}(i) = ti(i);
        w =  squeeze(UpsampledWaveform(toff(i):UpsamplingFactor:end,n,:));
        JitteredPeaks(k,:) = max(w)';
        ClusterIDs(k) = n;
        k = k + 1;
        idx = (1:T);
        if (idx(1) + ti(i) < 1) | (idx(end) + ti(i) > TD)
            idx( (idx + ti(i) < 1) | (idx + ti(i) > TD) ) = []; 
        end
        SimData(ti(i) + idx, :) = SimData(ti(i) + idx, :) + w(idx,:);
    end
end


