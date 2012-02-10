% function Mout = build_stat_arrays(Sdata, Slist, Window)
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

% temp
Sdata = bbnData;
Slist = validBBNList;

%-----------------------------------------------------------------------------
% some defaults
%-----------------------------------------------------------------------------
CONDITIONS = [1 2 3];
WINDOW_MS = 50;

%-----------------------------------------------------------------------------
% parse inputs
%-----------------------------------------------------------------------------
if ~exist('Clist', 'var')
	Clist = [1 2 3];
end
if ~exist('Window', 'var')
	Window = [0 WINDOW_MS];
elseif length(Window) == 1
	Window = [0 Window];
end

%-----------------------------------------------------------------------------
% get # of conditions 
%-----------------------------------------------------------------------------
Nconditions = length(Clist);

%-----------------------------------------------------------------------------
% create the struct for each condition
%-----------------------------------------------------------------------------
for c = 1:Nconditions
	Sc{c} = build_struct_bycondition(Sdata, Slist, Clist(c));
end

%-----------------------------------------------------------------------------
% get spike rates for the window
%-----------------------------------------------------------------------------

for l = 1:length(Sc{1});





end
	






























%-----------------------------------------------------------------------------
% convert to matrix
%
% this will be a 9 X 1 X Nitems cell array
%-----------------------------------------------------------------------------
Smat = struct2cell(S);






















