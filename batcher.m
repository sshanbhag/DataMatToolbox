%------------------------------------------------------------------------
% batcher.m
%------------------------------------------------------------------------
% script for spike rate computation and analysis from DataMat toolbox
%------------------------------------------------------------------------
% Some Notes:
%------------------------------------------------------------------------

clear all

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% define some options
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

inpath = '/Users/sshanbhag/Work/Data/LFHData/TextFiles';

types_to_process = {'BBN', 'LFH'};

outpath = '/Users/sshanbhag/Work/Data/LFHData/MatFiles';

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% File Information
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% get list of files, store info in fstruct struct array
fstruct = dir(inpath);

% make sure there are files
if isempty(fstruct)
	error('%s: no files found in %s', mfilename, inpath);
end

% find files that aren't empty
tmp = struct2cell(fstruct);
bytes = cell2mat(tmp(3, :));
clear tmp

% store valid files
[tmp, valid_findex] =  find(bytes > 0);
nfiles = length(valid_findex);

% loop through valid files
% for findex = 1:nfiles
for findex = 1:2
	% set procFlag to 0
	procFlag = 0;
	
	% build full filename (with path)
	fname = fstruct(valid_findex(findex)).name;
	fullname = fullfile(inpath, fname);

	% check if this is one of the types to process
	for n = 1:length(types_to_process)
		if strfind(fname, types_to_process{n})
			procFlag = 1;
			break;
		end
	end
	
	if procFlag
		% run loadDWfile
	 	[D, S, B] = loadDWfile(fullname, 'PlotData', 'No', 'SaveMatfile', 'Yes', 'OutputPath', outpath);	
	end
end






