function Rstruct = export_data_by_condition(exFile, exData, exList, spikeWin, excludeCList)
%--------------------------------------------------------------------------
% Rstruct = export_data(exFile, exData, exList, spikeWin)
%--------------------------------------------------------------------------
% 
% Rdata			[length(exList),  9 + 2*length(spikeWin)] cell array
% 
% 	column			Data
% 		1					filename
% 		2					condition id number
% 		3					unit #
% 		4					probe #
% 		5					cluster #
% 		6					background spike count
% 		7					attenuation value
% 		8					Net spike count mean (over entire sweep)
% 		9					Nes spike count std. deviation
% 		
% 		for each value of spikeWin (spike window) there will be 2
% 		columns of data:
% 		
% 			mean count in spike window
% 			std. dev. of spike count within spike window
%--------------------------------------------------------------------------

if ~exist('excludeCList', 'var')
	excludeCList = -1;
end

Nwin = length(spikeWin);
% assume same # of conditions in the exList files
Nconditions = length(exList{1, 1});

% open file
fp = fopen(exFile, 'w');

% write header information
fprintf(fp, 'N_SpikeCountWindows:,%d\n', Nwin);
fprintf(fp, 'Start_ms,End_ms,\n');
for w = 1:Nwin
	fprintf(fp, '%d,%d,\n', spikeWin{w}(1), spikeWin{w}(2) );
end

% write column headers

for c = 1:Nconditions
	prefix = sprintf('C%d_', c);
	fprintf(fp, '%sFile,', prefix);
	fprintf(fp, '%sUnitNum,', prefix);
	fprintf(fp, '%sProbe,', prefix);
	fprintf(fp, '%sCluster,', prefix);
	fprintf(fp, '%sBGcount,', prefix);
	fprintf(fp, '%sdBAtten,', prefix);
	fprintf(fp, '%sNet_mean,', prefix);
	fprintf(fp, '%sNet_sd,', prefix);
	% loop through spike count windows
	for w = 1:Nwin
		fprintf(fp, '%sCount_mean_w%d,', prefix, w);
		fprintf(fp, '%sCount_sd_w%d,', prefix, w);
	end
end
fprintf(fp, '\n');

rindx = 1;
for f = 1:length(exList)
	
	validDataIndices = exList{f, 1};
	validUnitIndices = exList{f, 5};
	validAttenIndices = exList{f, 7};
	
	[Nunits, tmp] = size(validUnitIndices);

	for u = 1:Nunits
		
		lfFlag = 0;
		for c = 1:Nconditions

			tmpD = exData{validDataIndices(c)};

			% columns in validUnitIndices correspond to conditions/files
			vUI = validUnitIndices(u, c);
			tmpU = tmpD.UnitData(vUI);

			if ~any(tmpU.UnitInfo.cluster == excludeCList)
				
				lfFlag = 1;

				attIndex = validAttenIndices(1, c);

				fprintf(fp, '%s,', tmpD.Info.file);
				fprintf(fp, '%d,', tmpU.UnitInfo.unit);
				fprintf(fp, '%d,', tmpU.UnitInfo.probe);
				fprintf(fp, '%d,', tmpU.UnitInfo.cluster);
				fprintf(fp, '%.4f,', tmpU.BG_count);

				fprintf(fp, '%.4f,', tmpD.AttenVals(attIndex));
				fprintf(fp, '%.4f,', tmpU.Net_mean(attIndex));
				fprintf(fp, '%.4f,', tmpU.Net_std(attIndex));

				% loop through spike count windows
				for w = 1:Nwin
					fprintf(fp, '%.4f,', tmpU.Count_mean{attIndex}(w));
					fprintf(fp, '%.4f,', tmpU.Count_std{attIndex}(w));
				end

				% store in Rstruct array
				Rstruct(rindx).file{c} = tmpD.Info.file;
				Rstruct(rindx).condition(c) = tmpD.Info.condition;
				Rstruct(rindx).unit(c) = tmpU.UnitInfo.unit;
				Rstruct(rindx).probe(c) = tmpU.UnitInfo.probe;
				Rstruct(rindx).cluster(c) = tmpU.UnitInfo.cluster;
				Rstruct(rindx).BG_count(c) = tmpU.BG_count;
				Rstruct(rindx).BG_tstamps{c} = tmpU.BG_timestamps;
				Rstruct(rindx).AttenVals(c) = tmpD.AttenVals(attIndex);
				Rstruct(rindx).Net_mean(c) = tmpU.Net_mean(attIndex);
				Rstruct(rindx).Net_std(c) = tmpU.Net_std(attIndex);
				% loop through spike count windows
				for w = 1:Nwin
					Rstruct(rindx).Count_mean{c}(w) = tmpU.Count_mean{attIndex}(w);
					Rstruct(rindx).Count_std{c}(w) = tmpU.Count_std{attIndex}(w);
					Rstruct(rindx).Counts{c, w} = tmpU.counts{attIndex}(w, :);
				end		% END w LOOP
				
				Rstruct(rindx).BG{c} = tmpU.BG;
				Rstruct(rindx).BG_mean{c} = tmpU.BG_mean;
				Rstruct(rindx).BG_std{c} = tmpU.BG_std;
				
				
			end	% end if
		end		% END c LOOP
		if lfFlag
			fprintf(fp, '\n');
			rindx = rindx + 1;
		end
	end		% END u LOOP
end		%END f LOOP
if fp ~= 1
	fclose(fp);
end

