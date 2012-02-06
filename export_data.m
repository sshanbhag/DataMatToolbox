function [Rdata, Rstruct] = export_data(exFile, exData, exList, spikeWin)
%--------------------------------------------------------------------------
% [Rdata, Rstruct] = export_data(exFile, exData, exList, spikeWin)
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


Nwin = length(spikeWin);

fp = fopen(exFile, 'w');

fprintf(fp, 'N_SpikeCountWindows:,%d\n', Nwin);
fprintf(fp, 'Start_ms,End_ms,\n')
for w = 1:Nwin
	fprintf(fp, '%d,%d,\n', spikeWin{w}(1), spikeWin{w}(2) );
end


fprintf(fp, 'File,');
fprintf(fp, 'Condition,');
fprintf(fp, 'UnitNum,');
fprintf(fp, 'Probe,');
fprintf(fp, 'Cluster,');
fprintf(fp, 'BGcount,');
fprintf(fp, 'Atten_dB,');
fprintf(fp, 'Net_mean,');
fprintf(fp, 'Net_sd,');
% loop through spike count windows
for w = 1:Nwin
	fprintf(fp, 'Count_mean_w%d,', w);
	fprintf(fp, 'Count_sd_w%d,', w);
end
fprintf(fp, '\n');


rindx = 0;
for f = 1:length(exList)
	
	validDataIndices = exList{f, 1};
	validUnitIndices = exList{f, 5};
	validAttenIndices = exList{f, 7};
	
	Nconditions = length(validDataIndices);
	[Nunits, tmp] = size(validUnitIndices);

	for c = 1:Nconditions
		tmpD = exData{validDataIndices(c)};
		attIndex = validAttenIndices(1, c);

		for u = 1:Nunits
			% columns in validUnitIndices correspond to conditions/files
			vUI = validUnitIndices(u, c);
			tmpU = tmpD.UnitData(vUI);

			fprintf(fp, '%s,', tmpD.Info.file);
			fprintf(fp, '%d,', tmpD.Info.condition);
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
			fprintf(fp, '\n');
		
			% store in Rdata cell matrix
			rindx = rindx + 1;
			col = 1;
			Rdata{rindx, col} = tmpD.Info.file;	col = col + 1;
			Rdata{rindx, col} = tmpD.Info.condition;	col = col + 1;
			Rdata{rindx, col} = tmpU.UnitInfo.unit;	col = col + 1;
			Rdata{rindx, col} = tmpU.UnitInfo.probe;	col = col + 1;
			Rdata{rindx, col} = tmpU.UnitInfo.cluster;	col = col + 1;
			Rdata{rindx, col} = tmpU.BG_count;	col = col + 1;
			Rdata{rindx, col} = tmpD.AttenVals(attIndex);	col = col + 1;
			Rdata{rindx, col} = tmpU.Net_mean(attIndex);	col = col + 1;
			Rdata{rindx, col} = tmpU.Net_std(attIndex);	col = col + 1;
			% loop through spike count windows
			for w = 1:Nwin
				Rdata{rindx, col} = tmpU.Count_mean{attIndex}(w);	col = col + 1;
				Rdata{rindx, col} = tmpU.Count_std{attIndex}(w);	col = col + 1;
			end		% END w LOOP
		end		% END u LOOP
	end		% END c LOOP
end		%END f LOOP
if fp ~= 1
	fclose(fp);
end


rindx = 0;
for f = 1:length(exList)
	
	validDataIndices = exList{f, 1};
	validUnitIndices = exList{f, 5};
	validAttenIndices = exList{f, 7};
	
	Nconditions = length(validDataIndices);
	[Nunits, tmp] = size(validUnitIndices);

	rindx = rindx + 1;
	for c = 1:Nconditions
		tmpD = exData{validDataIndices(c)};
		attIndex = validAttenIndices(1, c);

		for u = 1:Nunits
			% columns in validUnitIndices correspond to conditions/files
			vUI = validUnitIndices(u, c);
			tmpU = tmpD.UnitData(vUI);
			
			if tmpU.UnitInfo.cluster ~= 0
		
				% store in Rstruct struct array
				Rstruct(rindx).file{c} = tmpD.Info.file;
				Rstruct(rindx).condition(c) = tmpD.Info.condition;
				Rstruct(rindx).unit(c) = tmpU.UnitInfo.unit;
				Rstruct(rindx).probe(c) = tmpU.UnitInfo.probe;
				Rstruct(rindx).cluster(c) = tmpU.UnitInfo.cluster;
				Rstruct(rindx).BG_count(c) = tmpU.BG_count;
				Rstruct(rindx).AttenVals(c) = tmpD.AttenVals(attIndex);
				Rstruct(rindx).Net_mean(c) = tmpU.Net_mean(attIndex);
				Rstruct(rindx).Net_std(c) = tmpU.Net_std(attIndex);
				% loop through spike count windows
				for w = 1:Nwin
					Rstruct(rindx).Count_mean{c}(w) = tmpU.Count_mean{attIndex}(w);
					Rstruct(rindx).Count_std{c}(w)= tmpU.Count_std{attIndex}(w);
				end		% END w LOOP
			end
	
		end		% END u LOOP
	end		% END c LOOP
end		%END f LOOP



