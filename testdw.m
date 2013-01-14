clear all
close all

%------------------------------------------------------------
%% file names
%------------------------------------------------------------
DLLName = 'C:\DataWave\DWShared\nsDWFile.dll';
datapath = 'F:\Work\Data\MG';
filename = '12-12-2012--2854_BBNrate.ddf';
% filename = '12-12-2012--2854_FreqScan2.ddf';
% filename = '12-12-2012--2854_strings block.ddf';
% filename = '12-12-2012--2854_RepRate0.ddf';
% filename = '01-03-2013--2961_syllable_block_new_sorted.ddf';

%------------------------------------------------------------
%% init DW struct
%------------------------------------------------------------
D = DW.DWdata(fullfile(datapath,filename));

%------------------------------------------------------------
%% load markers
%------------------------------------------------------------
Events = D.loadMarkers;

%------------------------------------------------------------
%% parse markers
%------------------------------------------------------------
field_names = cell(D.DDF.nEvent, 1);
field_vals = cell(D.DDF.nEvent, 1);
for v = 1:D.DDF.nEvent
	tmp = textscan(Events(v).Info.CSVDesc', '%s', ...
									'Delimiter', ',', ...
									'MultipleDelimsAsOne', 1);
	field_names{v} = tmp{1};
	field_vals{v} = cell(Events(v).EventCount, 1);	
	if Events(v).EventCount
		for n = 1:Events(v).EventCount
			tmp = textscan(Events(v).Data{n}', '%s', ...
										'Delimiter', ',', ...
										'MultipleDelimsAsOne', 1);
			field_vals{v}{n} = tmp{1};
		end	% END n
		

	end	% END if EventCount
end	% END v DDF.nEvent index

for v = 1:length(field_names)
	fprintf('event %d:\n', v)
	disp(field_names{v})
	if ~isempty(field_vals{v})
		disp(field_vals{v}{1})
	end
	fprintf('\n\n');
end
