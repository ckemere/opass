function [x] = SimulatePosition(t, Velocity, TrackLength, EatingTime)
% function [x] = SimulatePosition(t, Velocity, TrackLength, EatingTime)
% Simulates running at a constant velocity on a linear track, pausing at each end to eat

RunningTime = TrackLength/Velocity;
TrackCycle = fix(t / (RunningTime + EatingTime));
Direction = 1 - 2*bitand(1, TrackCycle); % 1 => forward, -1 => backward
TrackPhase = t - TrackCycle * (RunningTime + EatingTime);
TrackPhase(TrackPhase > RunningTime) = RunningTime;
x = Direction.*TrackPhase.*Velocity + (Direction == -1).*TrackLength;

