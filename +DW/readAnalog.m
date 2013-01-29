function [d, s] = readAnalog(Hfile, Data, varargin)
%------------------------------------------------------------------------
%[d, s] = readAnalog(Hfile, Data)
%------------------------------------------------------------------------
% 	Analog Entities (2)
% 		Continuous, sampled data that represent digitized
% 		analog signals such as position, force, and other experiment signals, as
% 		well as electrode signals such as EKG, EEG and extracellular
% 		microelectrode recordings.
% 		Each index of an analog entity refers to a specific digitized sample.
% 		Each analog entity contains samples from a single channel and the
% 		number of indexes is equal to the number of samples present for that
% 		channel.
%------------------------------------------------------------------------
% Input Arguments:
% 	There are two ways that are used for reading analog data
% 	Analog = readAnalog(Hfile, Data, <entity list (optional)>)
% 		if no entity list (values from Data.AnalogList) is provided,
% 		all analog entities will be read.
% 		Warning: This might use up all memory if the recordings are
% 		of a long duration!
% 
% 	Analog = readAnalog(Hfile, Data, <entity list>, <start point>, < points
%				to read> )
% Output Arguments:
%------------------------------------------------------------------------
% See also: Neuroshare MATLAB API
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 2 January, 2013 (SJS)
% Revisions:
%	10 Jan 2013 (SJS): 
%	 -	renamed to readAnalog
%	 -	added optional input to allow selection of specific analog entities
%	 -	added options to allow selection of start point and # of points
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------
% check entities
if (Data.nAnalog == 0)
    disp('No analog entities available!');
	 d = [];
	 return
end
% check if user asked for a specific list of analog entities
anaList = [];
anaPoints = {};
if length(varargin) == 1
	anaList = varargin{1};
elseif length(varargin) == 3
	anaList = varargin{1};
	anaPoints = {varargin{2}, varargin{3}};
else
	fprintf('error: incorrect input args\n');
	d = [];
	return
end
% if anaList is empty, use everything in Data.AnalogList
if isempty(anaList)
	anaList = Data.AnalogList;
end

% pre-allocate Analog and s arrays
Nlist = length(anaList);
d = repmat(	struct(	'Info',		[], ...
							'startPoint',	[], ...
							'nPoints',	[], ...
							'Data',		[] ), ...
							Nlist, 1	);
s = zeros(Nlist, 1);

if isempty(anaPoints)
	% use all points
	anaPoints{1} = ones(Nlist, 1);
	anaPoints{2} = zeros(Nlist, 1);
	for n = 1:Nlist
		anaPoints{2}(n) = Data.EntityInfo(anaList(n)).ItemCount;
	end
end
% read analog entity information
for n = 1:Nlist
	% get the information about the current analog entity
	[s(n), d(n).Info] = ns_GetAnalogInfo(Hfile, anaList(n));
	% read analog entity data
	[s(n), nRead, d(n).Data] = ...
					ns_GetAnalogData(	Hfile, ...
											anaList(n), ...
											anaPoints{1}(n), ...
											anaPoints{2}(n)	);
	% store # of analog points (items) for this entity
	d(n).nPoints = length(d(n).Data);
	if d(n).nPoints ~= anaPoints{2}(n)
		fprintf('%s warn: nPoints ~= points requested\n', mfilename);
	end
	d(n).startPoint = anaPoints{1}(n);							
end

	