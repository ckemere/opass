function [SpikeTimes, NeuronProperties] = simulateSpikeTimes(NumNeurons, TotalTime)
% function [SpikeTimes] = simulateSpikeTimes(NumNeurons, TotalTime)

% NumNeurons - number of neurons to simulate
% TotalTime - total length of simulation in time

NeuronStats.PlaceFieldFiringRate = 20; % Spikes per pass
NeuronStats.PlaceFieldWidth = 20; % cm
NeuronStats.RippleFiringRate = 10; % Hz
NeuronStats.RippleFiringProbability = 0.02; % per ripple


% We assume that neurons behave according to 2 patterns - ripples and running
%  Ripples happen at 0 rad (reward site) each cycle


TrackRadius = 35; % cm -> ~2 m circumference
RunningSpeed = 25; % cm/s
RadialRunningSpeed = RunningSpeed/TrackRadius; % rad/s = 2*pi*RunningSpeed/(2*pi*TrackRadius);
NeuronStats.PlaceFieldConcentration = 1 / (NeuronStats.PlaceFieldWidth/TrackRadius)^2;
NeuronStats.PeakFiring = NeuronStats.PlaceFieldFiringRate * RunningSpeed / sqrt(2*pi) / ...
    NeuronStats.PlaceFieldWidth;
NeuronStats.ScaleFactorVonMises = NeuronStats.PeakFiring * RunningSpeed / ...
    (2*pi * besselk(0, NeuronStats.PlaceFieldConcentration));
NeuronStats.MaximumRateVonMises = NeuronStats.ScaleFactorVonMises * exp(NeuronStats.PlaceFieldConcentration);

EatingTime = 5; % seconds
RippleRate = 2; % per second - Poisson process


PositionFunction = @(t) RadialRunningSpeed * t;

% Place cell firing
for n = 1:NumNeurons
    NeuronProperties(n) = generateRandomNeuron(NeuronStats);
end

for n = 1:NumNeurons
    poissonSpikes = poissrnd(NeuronProperties(n).MaximumRateVonMises * TotalTime);
    spikeTimes = rand(poissonSpikes,1) * TotalTime;
    deletionFlag = rand(poissonSpikes,1) < NeuronProperties.RateFunction(PositionFunction(spikeTimes));
    SpikeTimes{n} = spikeTimes(deletionFlag);
end

function [NeuronProperties] = generateRandomNeuron(NeuronStats)
    NeuronProperties = NeuronStats;
    NeuronProperties.CenterAngle = rand(1,1) * 2*pi;
    NeuronProperties.RateFunction = @(x) ...  % x is an angle, so we'll make this a Von Mises. 
        NeuronStats.ScaleFactorVonMises * ...
        exp(NeuronStats.PlaceFieldConcentration .* cos(x - NeuronProperties.CenterAngle));

