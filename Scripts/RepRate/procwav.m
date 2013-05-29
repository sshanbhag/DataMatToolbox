%% procwav

%----------------------------------------------------------------------------
%% specify files, paths
%----------------------------------------------------------------------------
% wavfile path
wavpath = '~/Work/Data/RepRate/WavFiles';
% callinfo path
callinfo = '~/Work/Data/RepRate/CallInfo.mat';
% adj is machine is PC
if strcmpi(computer, 'PCWIN') | strcmpi(computer, 'PCWIN64')
	homepath = 'F:';
	wavpath = strrep(wavpath, '/', filesep);
	wavpath = strrep(wavpath, '~', homepath);
	callinfo = strrep(callinfo, '/', filesep);
	callinfo = strrep(callinfo, '~', homepath);
end

% get wavfile names
filelist = dir(fullfile(wavpath, '*.wav'));
nfiles = length(filelist);

% load CallInfo
load(callinfo);

%------------------------------------------------------------------------
%% parse names to get repetition rate value
%
% values with no _<value> portion are repetition rate = 1 Hz
%------------------------------------------------------------------------
RateVals = zeros(nfiles, 1);
StimNumber = zeros(nfiles, 1);
for f = 1:nfiles
	% get parts of filename
	[tmp, fstr, fext] = fileparts(filelist(f).name);
	% search for '_'
	underscore_locs = ('_' == fstr);
	if any(underscore_locs)
		% portion from underloc to fstr(end - 2) is reprate
		indx = find(underscore_locs);
		RateVals(f) = str2num(fstr( (indx(1)+1):(end - 2) ));
		StimNumber(f) = str2num(fstr( 1:(indx(1)-1) ));
	else
		% if no number is found, default rep rate is 1 Hz
		RateVals(f) = 1;
		StimNumber(f) = str2num(fstr);
	end
end

% will be easiest to find durations if CallInfo(:).id is numeric
tmp = struct2cell(CallInfo);
% first row holds id #s
callids = str2double(tmp(1, :));
% 4th row are durations, in seconds
calldurs = str2double(tmp(4, :));

% loop through files
for f = 1:nfiles

	% load wav file
	[wdata, fs] = wavread(fullfile(wavpath, filelist(f).name));

	% get duration for this file's signal
	dur = calldurs(StimNumber(f)==callids)
	
	% compute putative on/offset
	isi = 1/RateVals(f)
	nstim = 4;
	tmp = (0:(nstim-1)) * isi;
	teston = tmp(tmp<1); 
	testoff = teston + dur
	
	% plot signal
	figure(1)
	tvec = ((1:length(wdata)) - 1) ./ fs;
	plot(tvec, wdata);
	
	% find envelope
	env = abs(hilbert(wdata));
	env = abs(wdata);
	% plot overlay 
	hold on
	plot(tvec, env, 'r');
	hold off
	title(filelist(f).name, 'Interpreter', 'none')
	xlabel('time (seconds)')

	% downsample envelope
	env2 = decimate(env, fs*0.0001);
	t2 = ((1:length(env2)) - 1) ./ 10000;
	% plot downsampled envelope
	hold on
	plot(t2, env2, 'g')
	hold off

	% find values above threshold
	nz = env2>0.01;
	% append zeros to start/finish
	nz = [0 nz' 0];
	% diff the matrix
	dnz = diff(nz);
	% onsets are positive 1, offsets are -1, convert to seconds
	onsets = find(dnz==1) ./ 10000;
	offsets = find(dnz==-1) ./ 10000;

	figure(2)
	plot(t2, env2)
	hold on
	plot(((1:length(nz)) - 1) ./ 10000, nz, 'g')
	hold off

	offtol = .1;
	
	% loop through onsets
	onindex = 1;	% index for onsets
	offindex = 1;	% index for offsets
	detindex = 1;	% index for on/off detected vals
	while onindex < length(onsets)
		% check if next offset is < dur
		offFlag = 1;
		while offFlag
			if (offsets(offindex) - onsets(onindex)) ...
					<	(offtol *  dur)
				offindex = offindex + 1;
			else
				offFlag = 0;
				% store onset and offset
				det(detindex, :) = [onsets(onindex) offsets(offindex)];
				detindex = detindex + 1;
			end
		end
		% advance onset index
		onFlag = 1;
		
		while onFlag
			if onsets(onindex) <= offsets(offindex)
				onindex = onindex + 1;
			else
				onFlag = 0;
			end			
		end

	end
	
	L = ylim
	text(teston, L(2)*ones(size(teston)), '.', 'Color', 'g', 'FontSize', 25);
	text(testoff, L(2)*ones(size(testoff)), '.', 'Color', 'r', 'FontSize', 25);
		
	if length(onsets)~= length(offsets)
		error
	end	
	nreps = length(onsets);
	
	fp = 1;
	fprintf(fp, '%s,%d,%d,', filelist(f).name, StimNumber(f), RateVals(f))
	fprintf(fp, '%d,', nreps);
	for r = 1:nreps
		fprintf(fp, '%f,%f,%f,', onsets(r), offsets(r), offsets(r) - onsets(r));
	end
	fprintf(fp, '\n');
	
	pause
end