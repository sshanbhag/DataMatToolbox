%------------------------------------------------------------------------
% avgwaveforms
%------------------------------------------------------------------------

%------------------------------------------------------------------------
%  Sharad Shanbhag
%	sshanbhag@neomed.edu
%
% computes avg waveform for each file in unitdata 
%------------------------------------------------------------------------
% Created: 8 April 2013 (SJS)
%
% Revisions:
%	10 Apr 2013 (SJS): completed
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
% Settings for analysis
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% add DataMatToolbox to path
%------------------------------------------------------------------------
addpath(strrep('../DataMatToolbox', '/', filesep));

%------------------------------------------------------------------------
% File Settings
%------------------------------------------------------------------------
[a, b] = system('hostname');
if strncmpi(b, 'YM0360YHEUE-W7', length('YM0360YHEUE-W7'))
	dobjpath = 'F:\Work\Data\DataWave\batmat\dobj';
	unitmatfile = 'F:\Work\Data\DataWave\batmat\UnitInfo.mat';
	raw_outfile = 'F:\Work\Data\DataWave\batmat\RawWaves.mat';
	avg_outfile = 'F:\Work\Data\DataWave\batmat\AvgWaves.mat';
else
	dobjpath = '/Users/sshanbhag/Work/Data/DataWave/batmat/dobj';
	unitmatfile = '/Users/sshanbhag/Work/Data/DataWave/batmat/UnitInfo.mat';
	raw_outfile = '/Users/sshanbhag/Work/Data/DataWave/batmat/RawWaves.mat';
	avg_outfile = '/Users/sshanbhag/Work/Data/DataWave/batmat/AvgWaves.mat';
end
if ~exist('UnitInfo', 'var')
	load(unitmatfile);
end

fnames = UnitInfo(:, finddatacolumn('Filename', unitheader));
probes = UnitInfo(:, finddatacolumn('probeID', unitheader));
units = UnitInfo(:, finddatacolumn('unitnum', unitheader));

% need to do a little cleaning of filenames - remove all " from names
fnames = strrep(fnames, '"', '');

% avoid repeats - find unique fnames without sorting the names
[dobjnames, indices] = findUniqueUnsortedText(fnames);

% # of files
Nobjects = length(dobjnames);

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Analyze Data
%------------------------------------------------------------------------
%------------------------------------------------------------------------

winfo =	repmat(	struct(	'name',			'', ...
									'cluster',		[], ...
									'samprate',		[], ...
									'nclusters',	[], ...
									'nwindows',		[] ...
								), ...
						Nobjects, 1);
					
wprobes = zeros(Nobjects, 1);
wunits = zeros(Nobjects, 1);
wunitindx = zeros(Nobjects, 1);

wfields = fieldnames(winfo(1));
wforms = cell(Nobjects, 1);
wavg = cell(Nobjects, 1);
wstddev = cell(Nobjects, 1);
invertflag = zeros(Nobjects, 1);
wavg_invert = cell(Nobjects, 1);

%------------------------------------------------------------------------
% loop through file list
%------------------------------------------------------------------------
for N = 1:Nobjects
	% read in data using load command (.efd files are mat files)
	load(fullfile(dobjpath, dobjnames{N}), '-MAT');
	
	% convert probe and unit for this file to numbers
	probenum = str2double(probes(indices{N}));
	unitnum = str2double(units(indices{N}));
	if length(unique(probenum)) ~= 1
		error('more than 1 probenum for file %d', N);
	else
		P = probenum(1);
	end
	if length(unique(unitnum)) ~= 1
		error('more than 1 unitnum for file %d', N);
	else
		U = unitnum(1);
	end
	% copy info fields
	for f = 1:length(wfields);
		winfo(N).(wfields{f}) = d.Probes(P).(wfields{f}); 
	end
	% find index for unit waveform
	uindex = find(U == winfo(N).cluster);
	% store probe, unit and unitindex
	wprobes(N) = P;
	wunits(N) = U;
	wunitindx(N) = uindex;

	% convert waveforms cell array to numerical matrix
	wforms{N} = cell2mat(d.Probes(P).wforms{uindex});
	% compute mean and std. dev.
	wavg{N} = mean(wforms{N});
	wstddev{N} = std(wforms{N});
	
	% plot raw and avg and std dev
	figure(1)
	plot(wforms{N}');
	hold on
		plot(wavg{N}, 'k', 'LineWidth', 2);
		plot(wavg{N} + wstddev{N}, 'LineWidth', 1.1, 'Color', 0.75*[1 1 1]);
		plot(wavg{N} - wstddev{N}, 'LineWidth', 1.1, 'Color', 0.75*[1 1 1]);
	hold off
	tstr = {	d.fname, sprintf('Probe: %d', P), sprintf('Unit: %d', U) };
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

	% check if the max absolute value of the snippet prior to the peak
	% is greater than the peak; if so, invert the trace
	[premaxabs, premaxabsindx] = max(abs(wavg{N}(1:peakindx)));
	if peakval < premaxabs
		invertflag(N) = 1;
		flipval = -1;
	else
		flipval = 1;
	end
	wavg_invert{N} = flipval .* wavg{N};
	% plot orig and inverted
	subplot(212)
	plot(wavg_invert{N}, 'k');
	hold on
		plot(peakindx, flipval*peakval, 'ro')
		plot(premaxabsindx, premaxabs, 'go');
	hold off	
	
	drawnow
	clear d
end

save(raw_outfile, '-MAT', 'winfo', 'dobjnames', 'wforms', ...
									'wprobes', 'wunits', 'wunitindx');

save(avg_outfile, '-MAT', 'winfo', 'dobjnames', 'wavg', 'wstddev', ...
									'invertflag', 'wavg_invert', ...
									'wprobes', 'wunits', 'wunitindx');


