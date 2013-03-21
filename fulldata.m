%% load matfile (if FullData isn't already in workspace)
if ~exist('FullData', 'var')
	load('/Users/sshanbhag/Work/Data/DataWave/batmat/FullData6.mat');
end

%% Locate Complete, responsive units, with Type 1 data
%	keep DeleteRow == '0' ('0' means a complete test)  
%	keep Unresponsive Unit == '0' (keep responsive units)
%	keep Type == 1
fieldname = {'DeleteRow', 'Unresponsive Unit', 'Type'};
searchstr = {'0', '0', '1'};
[completeIndices, CompleteData] = finddata(fieldname, searchstr, FullData, header.fields);


%% get unique unit numbers
AllUnitNumbers = zeros(length(CompleteData), 1);
for n = 1:length(CompleteData)
	AllUnitNumbers(n) = str2num(CompleteData{n, 2});
end
% find unique unit numbers
UnitNumbers = unique(AllUnitNumbers);

%% get responsiveness 
attenvals = {'0', '20', '40'};
UnitAttenResp = cell(length(UnitNumbers), 1);
% loop through units
for n = 1:length(UnitNumbers)
	% find Complete data for current unit
	unitrows = find(strcmpi(num2str(UnitNumbers(n)), CompleteData(:, 2)))
	% Allocate responses for atten
	UnitAttenResp{n} = cell(3, 1);
	return
	
	for n = 1:length(attenvals)
		
		
	end
	
	
	
	
end



