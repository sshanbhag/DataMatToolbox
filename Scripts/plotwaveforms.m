%------------------------------------------------------------------------
% plotwaveforms
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%	sshanbhag@neomed.edu
%
% takes data from  avgwaveforms and plots it
%------------------------------------------------------------------------
% Created: 10 April 2013 (SJS)
%
% Revisions:
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% clear old bits
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% clear all
% close all

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% File Settings
%------------------------------------------------------------------------
%------------------------------------------------------------------------
[a, b] = system('hostname');
if strncmpi(b, 'YM0360YHEUE-W7', length('YM0360YHEUE-W7'))
	raw_outfile = 'F:\Work\Data\DataWave\batmat\RawWaves.mat';
	avg_outfile = 'F:\Work\Data\DataWave\batmat\AvgWaves.mat';
else
	raw_outfile = '/Users/sshanbhag/Work/Data/DataWave/batmat/RawWaves.mat';
	avg_outfile = '/Users/sshanbhag/Work/Data/DataWave/batmat/AvgWaves.mat';
end


%------------------------------------------------------------------------
%------------------------------------------------------------------------
% load Data
%------------------------------------------------------------------------
% These matfiles contain the following data structures:
% 		(all are of same length)
% 
%	dobjnames	dobj filenames
%	
% 	wprobes		probe for each unit/datafile, vector
% 	wunits		unit ID for each datafile, numeric vector
% 	wunitindx	index into d.Probes for desired unit in wunits
% 
% 	winfo		struct array taken from the d.Probes(probenum), each element having fields:
% 		name			sort waveform name from datawave/neuroshare
% 		cluster		cluster number (e.g., 1, 255)
% 		samprate		snippet sample rate in samples/sec
% 		nclusters	# of clusters found for this probe
% 		nwindows		# of windows (from neuroshare/datawave)
% 
% 	wforms		raw waveform data, cell array
% 	wavg			avg of wforms, cell array
% 	wstddev		stddev of wforms = cell(Nobjects, 1);
% 	invertflag	whether wavform in wavg_invert has been inverted (1 = yes, 0 = no)
% 	wavg_invert inverted (if necessary) average waveform
%------------------------------------------------------------------------

load(raw_outfile, '-MAT');
load(avg_outfile, '-MAT');

% # of files
Nobjects = length(dobjnames);

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% loop through file list and plots things
%------------------------------------------------------------------------
%------------------------------------------------------------------------
for N = 1:Nobjects
	figure(1)
	plot(wforms{N}');
	hold on
		plot(wavg{N}, 'k', 'LineWidth', 2);
		plot(wavg{N} + wstddev{N}, 'LineWidth', 1.1, 'Color', 0.75*[1 1 1]);
		plot(wavg{N} - wstddev{N}, 'LineWidth', 1.1, 'Color', 0.75*[1 1 1]);
	hold off
	tstr = {	dobjnames{N}, sprintf('Probe: %d', P), sprintf('Unit: %d', U) };
	title(tstr, 'Interpreter', 'none')
	
	% find peak value
	[peakval, peakindx] = max(wavg{N});
	
	% plot all traces and mean and mean +/- std. dev
	figure(2)
	subplot(211)
	plot(wavg{N}, 'k');
	hold on
		plot(peakindx, peakval, 'ro');
	hold off
	title(sprintf('wavg(%d)', N));
	
	% check if the max absolute value of the snippet prior to the peak
	% is greater than the peak; if so, invert the trace
	[premaxabs, premaxabsindx] = max(abs(wavg{N}(1:peakindx)));
	if invertflag(N) == 1;
		flipval = -1;
	else
		flipval = 1;
	end
	% plot orig and inverted
	subplot(212)
	plot(wavg_invert{N}, 'k');
	hold on
		plot(peakindx, flipval*peakval, 'ro')
		plot(premaxabsindx, premaxabs, 'go');
	hold off	
	title(sprintf('wavg_invert(%d), invertflag = %d', N, invertflag(N)), 'Interpreter', 'none');
	
	drawnow
	pause
end



