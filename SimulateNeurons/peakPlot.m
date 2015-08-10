function peakPlot(x, cluster)
% function peakPlot(x, cluster)
%
% Plot the points for signals from a tetrode in a pair-wise fashion

if (size(x,2) ~= 4)
    error('Expected Nx4 matrix of points');
end

xmax = 0;
ymax = 0;

k = 1;
for i = 1:4
    for j = i+1:4
        subplot(2,3,k)
        plot(x(:,i), x(:,j), 'o');
        xlabel(sprintf('Dimension %d',i));
        ylabel(sprintf('Dimension %d',j));
        k = k + 1;
        xlim = get(gca,'xlim');
        xmax = max(xmax,xlim(2));
        ylim = get(gca,'ylim');
        ymax = max(ymax,ylim(2));
    end
end

for i = 1:6
    subplot(2,3,i)
    set(gca,'xlim',[0, xmax]);
    set(gca,'ylim',[0, ymax]);
end
