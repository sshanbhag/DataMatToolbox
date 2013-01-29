function [d, s] = readEvent(Hfile, Data)
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% 	Event Entities (1)
% 		Discrete events that consist of small time-stamped text
% 		or binary data packets.  These are used to represent data such as trial
% 		markers, experimental events, digital input values, and embedded user
% 		comments.
% 		Each index of an event entity refers to a timestamp and data
% 		combination. The number of indexes is equal to the number of event
% 		entries for that event entity in the data file.
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
%	10 Jan 2013 (SJS): renamed to readEvent
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

if (Data.nEvent == 0)
    disp('No event entities available!');
	 d = [];
else
	s = zeros(Data.nEvent, 1);
	% Read Events in the eventlist
	for n = 1:Data.nEvent
		% get event information
		[s(n, 1), tmp.Info] = ns_GetEventInfo(Hfile, Data.EventList(n));
		% # of Items for this event
		tmp.EventCount = Data.EntityInfo(Data.EventList(n)).ItemCount;
		if tmp.EventCount
			% read data for this event if length ~= 0
			s2 = zeros(tmp.EventCount, 1);
			for m = 1:tmp.EventCount
				[s2(m), tmp.TimeStamp(m), tmp.Data(m), tmp.DataSize(m)] = ...
						ns_GetEventData(Hfile, Data.EventList(n), m);
			end
		else
			% otherwise, assign empty vals to data structure
			tmp.TimeStamp = [];
			tmp.Data = [];
			tmp.DataSize = [];
		end
		d(n) = tmp;
		clear tmp;
	end
end

