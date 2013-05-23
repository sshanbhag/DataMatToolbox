%----------------------------------------------------------------------------
% settings
%----------------------------------------------------------------------------
SAVEFLAG = 1;

%----------------------------------------------------------------------------
% specify files, paths
%----------------------------------------------------------------------------
toolpath = '~/Work/Code/Matlab/dev/Analysis/DataWave/DataMatToolbox';
unitpath =  '~/Work/Data/DataWave/batmat/UnitData.mat';
% path to data objects
dobjpath = '~/Work/Data/DataWave/batmat/dobj';
callinfopath = '~/Work/Data/RepRate/CallInfo.mat';
% output path
outputpath = '~/Work/Data/RepRate';
% adj is machine is PC
if strcmpi(computer, 'PCWIN') | strcmpi(computer, 'PCWIN64')
	homepath = 'F:';
	toolpath = strrep(toolpath, '/', filesep);
	toolpath = strrep(toolpath, '~', homepath);
	unitpath = strrep(unitpath, '/', filesep);
	unitpath = strrep(unitpath, '~', homepath);
	dobjpath = strrep(dobjpath, '/', filesep);
	dobjpath = strrep(dobjpath, '~', homepath);
	callinfopath = strrep(callinfopath, '/', filesep);
	callinfopath = strrep(callinfopath, '~', homepath);
	outputpath = strrep(outputpath, '/', filesep);
	outputpath = strrep(outputpath, '~', homepath);
end

%----------------------------------------------------------------------------
% set up paths if necessary
%----------------------------------------------------------------------------
if ~exist('finddata', 'file')
	disp('adding DataMatToolbox to path')
	addpath(toolpath)
end
%----------------------------------------------------------------------------
% load UnitData information (if UnitData isn't already in workspace)
%----------------------------------------------------------------------------
if ~exist('UnitData', 'var')
	load(unitpath);
end
%----------------------------------------------------------------------------
% load call information
%----------------------------------------------------------------------------
load(callinfopath);

%----------------------------------------------------------------------------
% other settings
%----------------------------------------------------------------------------
spikewin = [-200 800];

%----------------------------------------------------------------------------
% Locate Complete, responsive units, with Type 1 data
%	keep DeleteRow == '0' ('0' means a complete test)  
%	keep Unresponsive Unit == '0' (keep responsive units)
%	keep Type == 5 (reprate, atten = 0)
%	Type == 6 is reprate, atten = 20 dB
%----------------------------------------------------------------------------
fieldname = {'DeleteRow', 'Unresponsive Unit', 'Type'};
searchstr = {'0', '0', '5'};
% completeIndices are indices into UnitData, CompleteData are a subset
% of UnitData that match the search criteria specified by fieldname and
% searchstr
[completeIndices, CompleteData] = finddata(fieldname, searchstr, UnitData, unitheader.fields);

%----------------------------------------------------------------------------
% get filenames
%----------------------------------------------------------------------------
% find Dobj files
filecol = finddatacolumn('Filename', unitheader.fields);
all_filenames = CompleteData(:, filecol);
[filenames, findices, nfiles] = findUniqueUnsortedText(all_filenames);

%----------------------------------------------------------------------------
% some bookkeeping
%----------------------------------------------------------------------------
probecol = finddatacolumn('probeID', unitheader.fields);
unitcol = finddatacolumn('unitnum', unitheader.fields);
stimcol = finddatacolumn('stimulus', unitheader.fields);
attencol = finddatacolumn('atten', unitheader.fields);

%----------------------------------------------------------------------------
%% process files
%----------------------------------------------------------------------------
for n = 1:3

	% load file
	load(fullfile(dobjpath, filenames{n}), 'd');
	d.findVarsAndAtten;

	%------------------------------------------------------------------------
	% parse names to get repetition rate value
	%
	% values with no _<value> portion are repetition rate = 1 Hz
	%------------------------------------------------------------------------
	RateVals = zeros(d.Nvars, 1);
	StimNumber = zeros(d.Nvars, 1);
	for v = 1:d.Nvars
		% get parts of filename
		[tmp, fstr, fext] = fileparts(d.Variables{v});
		% search for '_'
		underscore_locs = ('_' == fstr);
		if any(underscore_locs)
			% portion from underloc to fstr(end - 2) is reprate
			indx = find(underscore_locs);
			RateVals(v) = str2num(fstr( (indx(1)+1):(end - 2) ));
			StimNumber(v) = str2num(fstr( 1:(indx(1)-1) ));
		else
			% if no number is found, default rep rate is 1 Hz
			RateVals(v) = 1;
			StimNumber(v) = str2num(fstr);
		end
	end
	
	%------------------------------------------------------------------------
	% get probe and unit number, get spikes for this atten and stim
	%------------------------------------------------------------------------
	tmp = str2num(cell2mat(CompleteData(findices{n}, probecol)));
	if length(unique(tmp)) ~= 1
		error('probenum error')
	else
		probenum = tmp(1);
	end
	tmp = str2num(cell2mat(CompleteData(findices{n}, unitcol)));
	if length(unique(tmp)) ~= 1
		error('unitnum error')
	else
		unitnum = tmp(1);
	end
	%------------------------------------------------------------------------
	% get atten
	%------------------------------------------------------------------------
	tmp = str2num(cell2mat(CompleteData(findices{n}, attencol)));
	if length(unique(tmp)) ~= 1
		error('atten error')
	else
		atten = tmp(1);
	end

	%------------------------------------------------------------------------
	%% now need to break things up by stimulus and then rep rate
	%------------------------------------------------------------------------
	% find unique stim numbers
	StimList = unique(StimNumber);
	if isempty(StimList)
		error('Empty StimList')
	end
	% loop through list of stims, find out reprates for each stim
	RatesForStim = cell(length(StimList), 1);
	IndicesForStim = cell(length(StimList), 1);
	for s = 1:length(StimList)
		current_stim = StimList(s);
		% locate current stim within master stim number list
		tmpindices = find(current_stim == StimNumber);
		tmpvals = RateVals(tmpindices);
		% store sorted rates and indices
		[RatesForStim{s}, sortindices] = sort(tmpvals);
		IndicesForStim{s} = tmpindices(sortindices);
	end
	
	SpikeData = cell(length(StimList), 1);
	for s = 1:length(StimList)
		tmp = cell(length(IndicesForStim{s}), 1);
		for t = 1:length(tmp)
			tmpspikes = d.getSpikesForStim(	IndicesForStim{s}(t), ...
														probenum,	...
														unitnum,		...
														'offset', spikewin);
			% convert to msec
			for n = 1:length(tmpspikes{1})
				tmpspikes{1}{n} = 0.001 .* tmpspikes{1}{n};
			end
			tmp{t} = tmpspikes{1};
		end
		SpikeData{s} = tmp;
	end
	
	%% plots
	% determine # of plots
	nstim = length(SpikeData);
	nrates = zeros(nstim, 1);
	for s = 1:nstim
		nrates(s) = length(RatesForStim{s});
	end
	% find max n rates
	maxnrates = max(nrates);

	% need to determine unique rate values - to do this, construct list
	% of all rates
	allrates = cell2mat(RatesForStim);
	uniquerates = unique(allrates);

	AllSpikes = cell(maxnrates, nstim);
	for s = 1:nstim
		for r = 1:nrates(s)
			% assign spikes to appropriate rates
			rindex = find(RatesForStim{s}(r) == uniquerates);
			if isempty(rindex)
				error('Unknown rate %d!', RatesForStim{s}(r));
			else
				AllSpikes{rindex, s} = SpikeData{s}{r};
			end
		end
		% if any of AllSpikes are empty, assign "empty" to them
		for r = 1:maxnrates
			if isempty(AllSpikes{r, s})
				AllSpikes{r, s} = 'empty';
			end
		end
	end
	
	plotopts.rowlabels = cell(length(uniquerates), 1);
	for rl = 1:length(uniquerates)
		plotopts.rowlabels{rl} = sprintf('%d Hz', uniquerates(rl));
	end
	plotopts.columnlabels = cell(nstim, 1);
	for s = 1:nstim
		plotopts.columnlabels{s} = CallInfo(StimList(s)).call		
	end

	plotopts.raster_tickmarker = '.';
	plotopts.raster_ticksize = 14;
	plotopts.vertgap = 0.015;
	plotopts.horizgap = 0.04;
	plotopts.plotgap = 0.002;
	
	plotopts.heightscale = .95;

	plotopts.vertoffset = 0.01;
	
	plotopts.timelimits = spikewin;
	
	plotopts.filelabel = d.fname;

	figure
	[H, opts] = rasterpsthmatrix(AllSpikes, plotopts);
	set(gcf, 'Position', [680 81 1100 900]);
	set(gcf, 'Name', d.fname);
	orient(gcf, 'landscape');
	set(gcf, 'PaperPosition', [0.25 0.1 10.5 8.1])
	drawnow
	
	
	if SAVEFLAG
		saveas(gcf, fullfile(outputpath, d.fname), 'fig');
		epsfile = fullfile(outputpath, [d.fname '.eps']);
		pdffile = fullfile(outputpath, [d.fname '.pdf']);
		print(gcf, '-depsc2', epsfile);
		system(sprintf('ps2pdf %s %s', epsfile, pdffile))
	end
	
	
%% end of loop
end	% END n loop



