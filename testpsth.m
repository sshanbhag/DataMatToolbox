mintime = -5;
maxtime = 10;
binsize = 1;

% add "dummy" bin - this will be removed 
maxtime2 = maxtime + binsize;
% create vector with maxtime "rounded" up to create even division by binsize
bins = mintime:binsize:(maxtime2 + mod(maxtime2, binsize));
nbins = length(bins);

spiketimes = 10*randn(1, 200);

H = histc(spiketimes, bins);

bins2 = bins(1:(nbins-1))
H2 = H(1:(nbins-1))


bins_old = mintime:binsize:maxtime + mod(maxtime, binsize);
H_old = histc(spiketimes, bins_old);
subplot(211)
bar(bins_old, H_old, 'histc')

subplot(212)
[pH, pb] = psth(spiketimes, binsize, [mintime maxtime])
bar(pb, pH, 'histc')

