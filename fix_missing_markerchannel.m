function tmpEvent = fix_missing_markerchannel(baseEvent, emptyEvent)

DW.DataWaveDefaults;

tmpEvent = baseEvent;

% create bogus marker data
tmpstr = '0';
for c = 2:MARKER_NBASE
	tmpstr = [tmpstr ',' '0']; %#ok<AGROW>
end
% tmpstr = '0,0,,0,0,0,0,0,0,0,0,0,0,0,0,0'
tmpsize = length(tmpstr);

% assign to temp Event data
for n = 1:tmpEvent.EventCount
	tmpEvent.Data{n} = tmpstr;
	tmpEvent.DataSize(n) = tmpsize;
end

% assign data Info bits
if emptyEvent == L
	start_indx = MARKER_NBASE + 1;
	end_indx = 2*MARKER_NBASE;
elseif emptyEvent == R
	start_indx = 1;
	end_indx = MARKER_NBASE;
end

tmpEvent.Info.CSVDesc = '';
loopindx = 1;
for n = start_indx:end_indx
	if loopindx == MARKER_NBASE
		tmpEvent.Info.CSVDesc = [tmpEvent.Info.CSVDesc MARKER_TAGS{n}]; %#ok<*USENS>
	else
		tmpEvent.Info.CSVDesc = [tmpEvent.Info.CSVDesc MARKER_TAGS{n} ', '];
	end
	
	loopindx = loopindx + 1;
end

