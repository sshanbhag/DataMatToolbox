clear all
close all
clear classes

%------------------------------------------------------------
%% file names
%------------------------------------------------------------
inputpath = '/Users/sshanbhag/Work/Data/DataWave/batmat/convertedDDF';
dobjoutpath = '/Users/sshanbhag/Work/Data/DataWave/batmat/dobjFRA';
figoutpath = '/Users/sshanbhag/Work/Data/DataWave/batmat/FRAplots'
unitmatfile = '/Users/sshanbhag/Work/Data/DataWave/batmat/UnitInfo.mat';

if ~exist('UnitInfo', 'var')
	load(unitmatfile);
end

FRAlist = dir(fullfile(inputpath, '*FRA*.mat'));
nFRA = length(FRAlist)

%% get probe and unit num
% filenames
fnames = UnitInfo(:, finddatacolumn('Filename', unitheader));
probes = UnitInfo(:, finddatacolumn('probeID', unitheader));
units = UnitInfo(:, finddatacolumn('unitnum', unitheader));

un_probes = cell(nFRA, 1);
un_units = cell(nFRA, 1);

fnum = 0
for n = 1:nFRA
	% get filename without .mat
	[~, tmp] = fileparts(FRAlist(n).name)
		
	ind = strfind(tmp, 'FRA')
	
	if isempty(ind)
		error('not found')
	end
	
	searchstr = tmp(1:(ind-2))
	
	% look for searchstr in fnames
	tmpfind = strfind(fnames, searchstr);
	vals = zeros(length(tmpfind), 1);
	for f = 1:length(tmpfind)
		if isempty(tmpfind{f})
			vals(f) = 0;
		else
			vals(f) = tmpfind{f};
		end
	end
	% locate nonzero elements
	indx = find(vals);
	
	% find unique probes
	un_probes{n} = unique(probes(indx));
	% find unique units
	un_units{n} = unique(units(indx));
end

for n = 1:length(FRAlist)
	
	for p = 1:length(un_probes{n})
		
		for u = 1:length(un_units{n})
			
			probenum = str2num(un_probes{n}{p});
			unitnum = str2num(un_units{n}{u});
			
			fprintf('file: %s \t probe: %d \t unit: %d\n', FRAlist(n).name, probenum, unitnum)
		end
	end
end



%------------------------------------------------------------
%% create Data Object
%------------------------------------------------------------
for n = 1:length(FRAlist)
	
	for p = 1:length(un_probes{n})
		
		for u = 1:length(un_units{n})
			
			probenum = str2num(un_probes{n}{p});
			unitnum = str2num(un_units{n}{u});
			
			fprintf('file: %s \t probe: %d \t unit: %d\n', FRAlist(n).name, probenum, unitnum)
	
			matfile = fullfile(inputpath, FRAlist(n).name);
			[tmppath, tmpfile, tmpext] = fileparts(FRAlist(n).name);
			tmpfile = sprintf('%s_Probe%d_Unit%d', tmpfile, probenum, unitnum);

			% load matfile to get D 
			fprintf('%s: loading .mat file %s\n', mfilename, matfile);
			load(matfile);

			%------------------------------------------------------------
			%------------------------------------------------------------
			%% convert to obj
			%------------------------------------------------------------
			%------------------------------------------------------------
			d = DW.FRAdata(D, matfile);
			
			%------------------------------------------------------------
			%------------------------------------------------------------
			%% plot unit waveforms (overlaid)
			%------------------------------------------------------------
			%------------------------------------------------------------
			d.plotUnitWaveforms('probe', probenum, 'unit', unitnum);
			drawnow	
			% save figure
			saveas(gcf, fullfile(figoutpath, [tmpfile '_wforms.fig']), 'fig')
			saveas(gcf, fullfile(figoutpath, [tmpfile '_wforms.jpg']), 'jpg')
			
			%------------------------------------------------------------
			%------------------------------------------------------------
			%% plot FRA
			%------------------------------------------------------------
			%------------------------------------------------------------
			d.plotFRA(probenum, unitnum, [0 800], 1);
			drawnow
			% save figure
			saveas(gcf, fullfile(figoutpath, [tmpfile '_FRAsurf.fig']), 'fig')
			saveas(gcf, fullfile(figoutpath, [tmpfile '_FRAsurf.jpg']), 'jpg')

			%------------------------------------------------------------
			%------------------------------------------------------------
			%% save Data Object
			%------------------------------------------------------------
			%------------------------------------------------------------
			% write the data object.  for debugging use, the D struct 
			% is also saved, but this will ultimately be redundant...
			dobjfile = fullfile(dobjoutpath, ['Dobj_' tmpfile '.mat'])
			save( dobjfile, 'd', 'probenum', 'unitnum');

		end	% END u nit
	end	% END p robe
end	% END n FRAfiles
