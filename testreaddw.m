

% plxfilename = '761_11SV_2204_3_sort.plx';
% dwfilename = '761_11SV_2204_3.txt';

datapath = sprintf('/Users/%s/Work/Data/DataWave/TestData/EXPFUNTESTER/', username);
plxfilename = [datapath 'test-single-04.plx'];
dwfilename = [datapath 'test-single-04.txt'];	

if ~exist('plxfiledetails', 'var') | ~exist('stimdata', 'var') | ~exist('rawstim', 'var')
	plxfiledetails = struct(	'noof_files',			1, ...
										'filename',				{{plxfilename}}, ...
										'startinterval',		1, ...
										'noof_intervals',		-1	);

	% load dw info
	[stimdata, vars, plx, rawdata] = loadDWStimData(plxfiledetails.filename{1}, dwfilename);
end



% vars is something like:
% 
% vars = 
% 
%          values: [180x1 double]
%           names: {'BBN level (dB SPL)'}
%     description: 'RIBBN'
% 
% >> vars.values
% 
% ans =
% 
%     10
%     70
%     10
%     40
%      0
%     20
%     30
% 	 .
% 	 .
% 	 .
% 	 etc. ...
% 
% 
% % copy plx data to plx
% plx = rawstim.plxdata;
% % add event timestamps to event_ts field in plx as cell
% plx.event_ts = {plx.event_ts}
% % specify event channel, counts, etc.
% plx.noof_event_chans = 1;
% plx.event_chan_id = [258];
% plx.event_ts_counts = [length(plx.event_ts{1})];
% 
% 
% % create textnum field in rawstim.marker struct to hold stim id numbers
% rawstim.marker.textnum = zeros(size(rawstim.marker.value));
% for n = 1:length(stimdata)
% 	stimindex = stimdata(n).indices;
% 	rawstim.marker.textnum(stimindex) = stim_idnums(n) * ones(size(stimindex));
% end
% 
% vars.values = [ rawstim.marker.textnum rawstim.marker.value ];
% vars.names = {'wavfile_id' 'atten_dB'};
% vars.description = 'Test Stimuli';

spikes = PLX2SpikeCount(plx,'var_matrix','variables',vars,'end', 1);

display_GenData(spikes, 'xdim',2, 'axesdimlist',1, 'minval',10,'title', vars.description);


% display_GenData(spikes, 'xdim',2, 'axesdimlist',1, 'conf95','normalize', 'xvalstep',1, 'minval',0,'title', vars.description);
% thr_conf95_excite = spikeCount2Thr(spikes,'increase','dimname',vars.names{2}, 'conf95',2,'interpolate')
% thr_conf95_inhibit = spikeCount2Thr(spikes,'decrease','dimname',vars.names{2},'conf95',2,'interpolate')
% display_genData(thr_conf95_excite,'xdim',2,'valuelinestyle','o-','title', vars.description);
