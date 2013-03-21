%% load UnitData information (if UnitData isn't already in workspace)
if ~exist('UnitData', 'var')
	load('/Users/sshanbhag/Work/Data/DataWave/batmat/UnitData.mat');
end
% path to objects
objpath = '/Users/sshanbhag/Work/Data/DataWave/batmat/dobj';

%% Locate Complete, responsive units, with Type 1 data
%	keep DeleteRow == '0' ('0' means a complete test)  
%	keep Unresponsive Unit == '0' (keep responsive units)
%	keep Type == 1
fieldname = {'DeleteRow', 'Unresponsive Unit', 'Type'};
searchstr = {'0', '0', '1'};
[completeIndices, CompleteData] = finddata(fieldname, searchstr, UnitData, unitheader.fields);
 


%% find values for attenuation == '0' within complete data 
% (this could be combined in one step above, but this is done separately 
% as a demonstration
attenval = 0;
[atten0ind, atten0dat] = finddata('atten', num2str(attenval), CompleteData, unitheader.fields);

%% get filenames

% find Dobj files
filecol = finddatacolumn('Filename', unitheader.fields);
filenames = atten0dat(:, filecol);
nfiles = length(filenames);

%% some bookkeeping
probecol = finddatacolumn('probeID', unitheader.fields);
unitcol = finddatacolumn('unitnum', unitheader.fields);
stimcol = finddatacolumn('stimulus', unitheader.fields);

%% load object

% select first file
n = 1;

% load d object from mat file
load(fullfile(objpath, filenames{n}), 'd');

% probe and unit for this data entry
probenum = str2num(atten0dat{n, probecol})
unitnum = str2num(atten0dat{n, unitcol})

% stimulus
stimulus = atten0dat{n, stimcol};

% need to do some manipulating of stimulus in order to get it to
% match up with the way that the object classifies stimuli
%
% main thing is to see if first part of stimulus is 'Noise' - if not, 
% append '.wav' to the stimulus (since it is a wav file)
if isempty(strfind(stimulus, 'Noise'))
	stimulus = [stimulus '.wav'];
end

fprintf('Getting spikes for probe %d, unit %d, stimulus %s, atten %d...\n', ...
					probenum, unitnum, stimulus, attenval);
% get spikes for this probe, unit, stimulus, atten level, over interval [-100
% 900]
S = d.getSpikesForVarAndAtten(probenum, unitnum, stimulus, attenval, [-100 900])



