function [d, s] = readNeural(Hfile, Data)
%------------------------------------------------------------------------
%[d, s] = readNeural(Hfile, Data)
%------------------------------------------------------------------------
% 	Neural Event Entities (4)
% 		Timestamps of event and segment entitities that
% 		are known to represent neural action potential firing times.  For
% 		example, if a segment entity contains sorted neural spike waveforms, each
% 		sorted unit is also exported as a neural entity.
% 		Each index of a neural event entity refers to a timestamp for each
% 		neural event.  Each neural event entity contains event times for a
% 		single neural source.  The number of indexes is equal to the number
% 		of entries for that neural event entity.
%------------------------------------------------------------------------
% Input Arguments:
%
% Output Arguments:
%------------------------------------------------------------------------
% See also: Neuroshare MATLAB API
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 2 January, 2013 (SJS)
%
% Revisions:
%	10 Jan 2013 (SJS): renamed to readNeural
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

if Data.nNeural == 0
	disp('No Neural entities available');
	d = [];
else
	NeuralLabels = strvcat(Data.EntityInfo(Data.NeuralList).EntityLabel);

	for n = 1:Data.nNeural
		%{
		% Have to figure out which Neural entities correspond with the selected segment entities
		list = strmatch(EntityInfo(SegmentList(channel(cChannel))).EntityLabel, NeuralLabels, 'exact');
		%}
		% Retrieve the data
		[s, tmp.Info] = ns_GetNeuralInfo(Hfile, Data.NeuralList(n));
		[s, tmp.Data] = ns_GetNeuralData(Hfile, Data.NeuralList(n), 1, Data.EntityInfo(Data.NeuralList(n)).ItemCount);

		indx = 1:(size(tmp.Data, 1) * size(tmp.Data, 2));

		% Get the neural timestamps
		tmp.Data = reshape(tmp.Data, size(tmp.Data, 1) * size(tmp.Data, 2), 1);
		tmp.TimeStamps = {tmp.Data(indx)};
		% Match the neural events with their unit ID
% 		tempID = ones(length(indx) / length(list), 1) * [tmp.Info(:).SourceUnitID];
% 		tmp.Units = {tempID(:)};
		% Remember how many neural events were found
		tmp.ItemCount = length(tmp.TimeStamps);

		d(n) = tmp;
		clear tmp;
	end
end
