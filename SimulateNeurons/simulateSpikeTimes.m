function [SpikeTimes, NeuronProperties] = simulateSpikeTimes(NumNeurons, TotalTime, PositionFunction)
% function [SpikeTimes] = simulateSpikeTimes(NumNeurons, TotalTime)
% NumNeurons - number of neurons to simulate
% TotalTime - total length of simulation in time

NeuronStats.PlaceFieldFiring = 20; % Spikes per pass
NeuronStats.PlaceFieldWidth = 10; % cm
NeuronStats.RippleFiringRate = 10; % Hz
NeuronStats.RippleFiringProbability = 0.02; % per ripple


% Simulated environment: Linear track
TrackLength = 200; % cm
NeuronStats.PlaceFieldCenterRange = [0, TrackLength];

RunningSpeed = 25; % cm/s
NeuronStats.ScaleFactor = NeuronStats.PlaceFieldFiring * RunningSpeed / sqrt(2*pi) / ...
    NeuronStats.PlaceFieldWidth;

EatingTime = 5; % seconds

% We assume that neurons behave according to 2 patterns - ripples and running
%  Ripples happen at 0 rad (reward site) each cycle

RippleRate = 2; % per second - Poisson process

RippleRegionFunction = @(x) ( (x == 0) | (x == TrackLength));

% Place cell firing
for n = 1:NumNeurons
    NeuronProperties(n) = generateRandomNeuron(NeuronStats, RippleRegionFunction);
end

for n = 1:NumNeurons
    poissonSpikes = poissrnd(NeuronProperties(n).ScaleFactor * TotalTime);
    spikeTimes = rand(poissonSpikes,1) * TotalTime;
    deletionFlag = rand(poissonSpikes,1) < NeuronProperties(n).RateFunction( ...
        PositionFunction(spikeTimes, RunningSpeed, TrackLength, EatingTime) );
    SpikeTimes{n} = spikeTimes(deletionFlag);
    % NeuronProperties(n).FiringRate = NeuronProperties(n).RateFunction( ...
        % PositionFunction(linspace(0,TotalTime,100*TotalTime), RunningSpeed, TrackLength, EatingTime) );
end

function [NeuronProperties] = generateRandomNeuron(NeuronStats, RippleRegionFunction)
    NeuronProperties = NeuronStats;
    NeuronProperties.PlaceFieldCenter = rand(1,1) * diff(NeuronStats.PlaceFieldCenterRange) + ...
        NeuronStats.PlaceFieldCenterRange(1);

    NeuronProperties.RateFunction = @(x) ...  % x is an angle, so we'll make this a Von Mises. 
        ~RippleRegionFunction(x) .* (NeuronStats.ScaleFactor * ...
        exp(-1/(2*NeuronStats.PlaceFieldWidth^2) * (x - NeuronProperties.PlaceFieldCenter).^2 ));

