function calcbins(spikearray,width,varargin)
%calcbins sums number of spikes that fall into bins of given time width,
% starting offset, and step interval. Results can be returned as a dense or
% sparse array.
%
% For any of the following right hand side (RHS) syntax:
%
% **If a single variable is specified on the left hand side (LHS), a dense
% matrix will be returned - a matrix with every step represented by an
% array element whether or not it contains spikes.
% ex: dense_bins = calcbins(spikearray, width)
% 
% **If two variables are specified on the LHS, a pair of vectors
% representing a sparse matrix will be returned only for bins with spikes.
% The first matrix being sum of spikes, the second being the time the bin
% centers on.
% ex: [sparse_binsum, sparse_bintime] = calcbins(spikearray, width)
%
% calcbins(spikearray,width)
%   calculates bins of width 'width' with offset = width/2 and ending with
%   the last bin needed to count the spike found by max(spikearray).
%
% calcbins(spikearray,width,


%  status = T_PSTHmake opens a dialog asking for the .mat file name of a
%  data file output by loadDWfile.m.  It creates files for each unit
%  following the naming convention [aaa_bbb_ccc_ppp_uuu] where aaa =
%  something, bbb = something else, ccc = another thing, ppp = probe
%  number and uuu = unit. Each file contains the following structure:
%  PSTH_bins
%     sourceFile - .mat file spike times came from
%     probe
%     unit
%     binSums - an array of sums of spike times falling into 1ms bins
%        where the binSum(1) is the sum from time 0 to 1ms and binSum(end)
%        is the last bin required to represent the last spike time recorded

  fprintf(1,'Calculating bins. (should only occur once unless binsize changes)\n');
  obj.psthbinsexist = 1;
  numunits = size(obj.D.UnitData,2);
  for i = 1:numunits
    if obj.D.UnitData(i).sorted == 0
      obj.D.UnitData(i).psthbins = [];
    else
      
%The following (old version) loads two matrices containing sparse bins and
%correlated counts but does it by unit.
%      ts = obj.D.UnitData(i).timestamp * 1e-6;
%      intobin = floor(ts/obj.psthbinsize);
%      [obj.D.UnitData(i).psthbinsum, obj.D.UnitData(i).psthbins] = hist(intobin,unique(intobin));
%      obj.D.UnitData(i).psthbins = obj.D.UnitData(i).psthbins * obj.psthbinsize;

%The following loads sparse and nonsparse matrices in Stimulus with
%structure of Stimulus.Spiketimes - {Nunits x Nreps}
      ts = obj.D.UnitData(i).timestamp * 1e-6;
      intobin = floor(ts/obj.psthbinsize);
      [obj.D.UnitData(i).psthbinsum, obj.D.UnitData(i).psthbins] = hist(intobin,unique(intobin));
%Unsparsify
      allbins(1:numel(obj.D.UnitData(i).psthbins))=0;
      
      obj.D.UnitData(i).psthbins = obj.D.UnitData(i).psthbins * obj.psthbinsize;
    end%if
  end%for
end