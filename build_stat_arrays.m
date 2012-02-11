function [Mout, Sout, Mfields] = build_stat_arrays(Sdata, Slist, Clist, Window, BGWindow, BGTime)
%-----------------------------------------------------------------------------
% Sout = build_stat_arrays(Sin)
%-----------------------------------------------------------------------------
% 
% Description
% 
%-----------------------------------------------------------------------------
% Input Arguments:
%-----------------------------------------------------------------------------
% 	Snet		struct array with fields:
% 		file				original data file name
% 		condition		stimulus condition (odor)
% 		unit				unit # (from loadDWfile...)
% 		probe				probe #
% 		cluster			cluster #
% 		BG_timestamps	timestamps for background activity for this unit
% 								(from first 5 seconds of continuous data collected)
% 		atten				stimulus attenuation (dB)
% 		ntrials			# of trials for this stimulus (a.k.a. # sweeps)
% 		spikes			cell array of neural spike times in milliseconds
%-----------------------------------------------------------------------------
% 
% Output Arguments:
% 	Output	output info
%
%-----------------------------------------------------------------------------
% See also: 
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%-----------------------------------------------------------------------------
% Created: 9 February, 2012 (SJS)
%
% Revisions:
%-----------------------------------------------------------------------------
% TO DO:
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
% some defaults
%-----------------------------------------------------------------------------
CONDITIONS = [1 2 3];
WINDOW_MS = 50;
BGWINDOW_MS = 50;
BGTIME_MS = 5000;

%-----------------------------------------------------------------------------
% parse inputs
%-----------------------------------------------------------------------------
if ~exist('Clist', 'var')
	Clist = CONDITIONS;
end
if ~exist('Window', 'var')
	Window{1} = [0 WINDOW_MS];
end
if ~exist('BGWindow', 'var')
	BGWindow = BGWINDOW_MS;
end
if ~exist('BGTime', 'var')
	BGTime = BGTIME_MS;
end

Nwin = length(Window);
for w = 1:Nwin
	Window_size_ms(w) = diff(Window{w});
end

[BGwin, Nbgwin] = compute_windows(BGTime, BGWindow);

%-----------------------------------------------------------------------------
% get # of conditions 
%-----------------------------------------------------------------------------
Nconditions = length(Clist);

%-----------------------------------------------------------------------------
% create the struct for each condition
%-----------------------------------------------------------------------------
for c = 1:Nconditions
	Sc = build_struct_bycondition(Sdata, Slist, Clist(c));
	
	% count spikes
	for u = 1:length(Sc)
		Sc(u).validspikes = cell(Nwin, 1);
		Sc(u).spikecount = cell(Nwin, 1);
		Sc(u).spikerate = cell(Nwin, 1);
		
		for w = 1:Nwin	
			Sc(u).validspikes{w} = cell(Sc(u).ntrials, 1);
			Sc(u).spikecount{w} = zeros(Sc(u).ntrials, 1);
			Sc(u).spikerate{w} = zeros(Sc(u).ntrials, 1);

			for t = 1:Sc(u).ntrials
				% find spikes that are within the current window
				validSpikes = (Window{w}(1) <= Sc(u).spikes{t}) & ...
										(Sc(u).spikes{t} < Window{w}(2));
				% store spiketimes
				if ~isempty(validSpikes)
					Sc(u).validspikes{w}{t} = Sc(u).spikes{t}(validSpikes);
					Sc(u).spikecount{w}(t) = sum(validSpikes);
					Sc(u).spikerate{w}(t) = Sc(u).spikecount{w}(t) / ...
															(0.001 * Window_size_ms(w));
				end
			end	% END t LOOP			
		end	% END w LOOP
		
		
		Sc(u).BG_spiketimes = cell(Nbgwin, 1);
		Sc(u).BG_spikecount = zeros(Nbgwin, 1);
		Sc(u).BG_spikerate = zeros(Nbgwin, 1);
		for b = 1:Nbgwin
			% find spikes that are within the current window
			validSpikes = (BGwin{b}(1) <= Sc(u).BGspikes) & ...
										(Sc(u).BGspikes < BGwin{b}(2));
			% store spiketimes
			if ~isempty(validSpikes)
				Sc(u).BG_spiketimes{b} = Sc(u).BGspikes(validSpikes);
				Sc(u).BG_spikecount(b) = sum(validSpikes);
				Sc(u).BG_spikerate(b) = Sc(u).BG_spikecount(b) / ...
															(0.001 * BGWindow);
			end
		end	% END b LOOP
	end	% END u LOOP
	
	Sout{c} = Sc;
end	% END c LOOP

clear Sc

%-----------------------------------------------------------------------------
% convert to matrix
%
% this will be a 9 X 1 X Nitems cell array
%-----------------------------------------------------------------------------
Mout = cell(Nconditions, 1);
for c = 1:Nconditions
	Smat = struct2cell(Sout{c});
	fnames = fieldnames(Sout{c}(1));
	[nfields, tmp, nunits] = size(Smat);
	
	tmpM = cell(nunits, nfields);
	
	for u = 1:nunits
		for f = 1:nfields
			tmpM{u, f} = Sout{c}(u).(fnames{f});
		end
	end
	
	Mout{c} = tmpM;
end

Mfields = fnames;






















