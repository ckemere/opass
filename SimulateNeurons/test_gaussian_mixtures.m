% clear all;
% run geometricalActionPotentials;

Fs = SamplingRate

[b,a] = butter(2,600/15000,'high');
[bhigh,ahigh] = butter(1,11000/15000);

for i = 1:4
    % x(:,i) = filtfilt(bhigh,ahigh,filtfilt(b,a,SimData(:,i)));
end

% [timepoints, spikes] = detectspikes_thresh(x(1:30*Fs, :), 40, Fs, 40, 9);

D = 2; % Number of PCA components to use
numCh = size(SimData, 2);
for i=1:numCh
    % Each column is an observation
    [U, ~, ~] = svd(spikes(:,:,i), 'econ');
    A{i} = U(:, 1:D);
    pc{i} = spikes(:,:,i)'*U(:,1:D);
end

traindata = [pc{:}];
testdata = traindata(1:600,:);
traindata = traindata(601:end,:);

for k = 1:30
    success = 0;
    while ~success
      try
        success = 1;
        gmfit = gmdistribution.fit(traindata,k);
      catch
        success = 0;
      end
    end
    model{k} = gmfit;
    [~,nlogl(k)] = gmfit.posterior(testdata);

end
