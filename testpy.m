inpathname = '/Users/sshanbhag/Work/Data/Bat/BatRestrainedData/878_Sorted_3-28-2014'
infilename = '10-09-2013-878_2948_StringsBlock_Sorted.txt'
infile = fullfile(inpathname, infilename);

% load .txt file as D struct
D = DW.convertPyDDF2MAT(infile)

%% convert to obj

%d = DW.RateData(D, fullfile(inputpath, matfile));
d = DW.RateData(D);



%%
probenum = 1;
unitnum = 1;
psthwin = [-200 799];
binsize = 10;

% plot unit waveforms (overlaid)
d.plotUnitWaveforms('probe', probenum, 'unit', unitnum);

%%


d.plotRasterAndPSTH('probe', probenum, 'unit', unitnum, 'offset', [-200 0])
