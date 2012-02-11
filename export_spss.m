% 1		file
% 2		condition
% 3		unit
% 4		probe
% 5		cluster
% 6		BGspikes
% 7		atten
% 8		ntrials
% 9		spikes
% 10		validspikes
% 11		spikecount
% 12		spikerate
% 13		BG_spiketimes
% 14		BG_spikecount
% 15		BG_spikerate

BGSWEEPTIME_MS = 50;

Windows = {	[0 50] ; ...
				[50 100] ;	...
				[100 150]	;	...
				[150 200]	;	...
			};
WindowSizes = 50;

load statsdata.mat

%---------------------------------------------------------------------------
% data for the BG, BBN
%---------------------------------------------------------------------------
fp = fopen('BBN_SPSS_Background.csv', 'w');
fprintf(fp, 'UnitID,');
fprintf(fp, 'condition,');
fprintf(fp, 'unit,');
fprintf(fp, 'probe,');
fprintf(fp, 'cluster,');
fprintf(fp, 'trial,');
fprintf(fp, 'spikecount,');
fprintf(fp, 'spikerate,');
fprintf(fp, '\n');

Nunits = length(Mbbnstruct{1});

% loop through conditions
for c = 1:3
	
	for u = 1:Nunits
		% make tmp copy of working struct
		M = Mbbnstruct{c}(u);
		
		% build the string of identification/codes
		UnitID = sprintf('%s_Unit_%d', M.file, M.unit);
		tmpstr = sprintf('%s,%d,%d,%d,%d', ...
									UnitID,	...
									M.condition, ...
									M.unit, ...
									M.probe, ...
									M.cluster ...
								);
		
		% # of background values
		nbgvalues = length(M.BG_spikecount);
		% get ntrials vector of random indexes into bgvalues
		tmpi = randi(nbgvalues, [M.ntrials 1]);
		tmpcount = M.BG_spikecount(tmpi);
		tmprate = tmpcount ./ (0.001 * BGSWEEPTIME_MS);
		
		% write to file
		for t = 1:M.ntrials
			fprintf(fp, '%s,%d,%d,%d,\n', tmpstr, t, tmpcount(t), tmprate(t));
		end
		
	end	% END u
end	% END c

if ~any(fp == [0 1])
	fclose(fp);
end

	
%---------------------------------------------------------------------------
% Real Data, BBN
%---------------------------------------------------------------------------
%{
fprintf(fp, 'UnitID,');
fprintf(fp, 'condition,');
fprintf(fp, 'unit,');
fprintf(fp, 'probe,');
fprintf(fp, 'cluster,');
fprintf(fp, 'BGspikes,');
fprintf(fp, 'atten,');
fprintf(fp, 'ntrials,');
fprintf(fp, 'spikes,');
fprintf(fp, 'validspikes,');
fprintf(fp, 'spikecount,');
fprintf(fp, 'spikerate,');
fprintf(fp, 'BG_spiketimes,');
fprintf(fp, 'BG_spikecount,');
fprintf(fp, 'BG_spikerate,');
%}

fp = fopen('BBN_SPSS.csv', 'w');
fprintf(fp, 'UnitID,');
fprintf(fp, 'condition,');
fprintf(fp, 'unit,');
fprintf(fp, 'probe,');
fprintf(fp, 'cluster,');
fprintf(fp, 'atten,');
fprintf(fp, 'window,');
fprintf(fp, 'trial,');
fprintf(fp, 'spikecount,');
fprintf(fp, 'spikerate,');
fprintf(fp, '\n');

Nunits = length(Mbbnstruct{1});

% loop through conditions
for c = 1:3
	
	for u = 1:Nunits
		% make tmp copy of working struct
		M = Mbbnstruct{c}(u);
		
		% build the string of identification/codes
		UnitID = sprintf('%s_Unit_%d', M.file, M.unit);
		tmpstr1 = sprintf('%s,%d,%d,%d,%d', ...
									UnitID,	...
									M.condition, ...
									M.unit, ...
									M.probe, ...
									M.cluster ...
								);
		
		for w = 1:length(M.spikecount)
			% append atten and window # to tmpstr
			tmpstr2 = sprintf('%s,%d,%d', tmpstr1, M.atten, w);
			tmpcount = M.spikecount{w};
			tmprate = M.spikerate{w};
			
			% write to file
			for t = 1:M.ntrials
				fprintf(fp, '%s,%d,%d,%d,\n', tmpstr2, t, tmpcount(t), tmprate(t));
			end
		end	% END w
	end	% END u
end	% END c

if ~any(fp == [0 1])
	fclose(fp);
end


%---------------------------------------------------------------------------
% data for the BG, LFH
%---------------------------------------------------------------------------
fp = fopen('LFH_SPSS_Background.csv', 'w');
fprintf(fp, 'UnitID,');
fprintf(fp, 'condition,');
fprintf(fp, 'unit,');
fprintf(fp, 'probe,');
fprintf(fp, 'cluster,');
fprintf(fp, 'trial,');
fprintf(fp, 'spikecount,');
fprintf(fp, 'spikerate,');
fprintf(fp, '\n');

Nunits = length(Mlfhstruct{1});

% loop through conditions
for c = 1:3
	
	for u = 1:Nunits
		% make tmp copy of working struct
		M = Mbbnstruct{c}(u);
		
		% build the string of identification/codes
		UnitID = sprintf('%s_Unit_%d', M.file, M.unit);
		tmpstr = sprintf('%s,%d,%d,%d,%d', ...
									UnitID,	...
									M.condition, ...
									M.unit, ...
									M.probe, ...
									M.cluster ...
								);
		
		% # of background values
		nbgvalues = length(M.BG_spikecount);
		% get ntrials vector of random indexes into bgvalues
		tmpi = randi(nbgvalues, [M.ntrials 1]);
		tmpcount = M.BG_spikecount(tmpi);
		tmprate = tmpcount ./ (0.001 * BGSWEEPTIME_MS);
		
		% write to file
		for t = 1:M.ntrials
			fprintf(fp, '%s,%d,%d,%d,\n', tmpstr, t, tmpcount(t), tmprate(t));
		end
		
	end	% END u
end	% END c

if ~any(fp == [0 1])
	fclose(fp);
end

	
%---------------------------------------------------------------------------
% Real Data, BBN
%---------------------------------------------------------------------------
%{
fprintf(fp, 'UnitID,');
fprintf(fp, 'condition,');
fprintf(fp, 'unit,');
fprintf(fp, 'probe,');
fprintf(fp, 'cluster,');
fprintf(fp, 'BGspikes,');
fprintf(fp, 'atten,');
fprintf(fp, 'ntrials,');
fprintf(fp, 'spikes,');
fprintf(fp, 'validspikes,');
fprintf(fp, 'spikecount,');
fprintf(fp, 'spikerate,');
fprintf(fp, 'BG_spiketimes,');
fprintf(fp, 'BG_spikecount,');
fprintf(fp, 'BG_spikerate,');
%}

fp = fopen('LFH_SPSS.csv', 'w');
fprintf(fp, 'UnitID,');
fprintf(fp, 'condition,');
fprintf(fp, 'unit,');
fprintf(fp, 'probe,');
fprintf(fp, 'cluster,');
fprintf(fp, 'atten,');
fprintf(fp, 'window,');
fprintf(fp, 'trial,');
fprintf(fp, 'spikecount,');
fprintf(fp, 'spikerate,');
fprintf(fp, '\n');

Nunits = length(Mlfhstruct{1});

% loop through conditions
for c = 1:3
	
	for u = 1:Nunits
		% make tmp copy of working struct
		M = Mlfhstruct{c}(u);
		
		% build the string of identification/codes
		UnitID = sprintf('%s_Unit_%d', M.file, M.unit);
		tmpstr1 = sprintf('%s,%d,%d,%d,%d', ...
									UnitID,	...
									M.condition, ...
									M.unit, ...
									M.probe, ...
									M.cluster ...
								);
		
		for w = 1:length(M.spikecount)
			% append atten and window # to tmpstr
			tmpstr2 = sprintf('%s,%d,%d', tmpstr1, M.atten, w);
			tmpcount = M.spikecount{w};
			tmprate = M.spikerate{w};
			
			% write to file
			for t = 1:M.ntrials
				fprintf(fp, '%s,%d,%d,%d,\n', tmpstr2, t, tmpcount(t), tmprate(t));
			end
		end	% END w
	end	% END u
end	% END c

if ~any(fp == [0 1])
	fclose(fp);
end


