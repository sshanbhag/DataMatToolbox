function [spiketimes, spikeisi] = poisson_spiketrain(Nspikes, L)
% returns a length Nspikes vector of spike ISI values in seconds with mean
% ISI of L milliseconds
if Nspikes == 0
	spiketimes = [];
	spikeisi = [];
	return
end

spikeisi = 1 + poissrnd(L, Nspikes+1, 1);
spiketimes = cumsum(0.001 * spikeisi);

	