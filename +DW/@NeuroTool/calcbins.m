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
%   calculates bins of width 'width' with offset = width/2 and step size
%   equal to width ending with the bin needed to count the spike found by
%   max(spikearray).
%
% calcbins(spikearray,width,offset)
%   performs exactly the same but begins at time offset specified
%
% calcbins(spikearray,width,offset,step)
%   step specifies the distance between center points of each bin such that
%   if a bin's width is greater than the step, an overlap will occur. Using
%   this, it is possible to implement a sliding window. Note: if step ig
%   greater than width, spikes can fall between bins and not be counted.
%
% See also: N
  
  fprintf(1,'Calculating bins.\n');
  %defaults
  offset = width / 2;
  step = width;
  switch nargin
    case 2
      %already did offset and step
    case 3
      offset = varargin(1);
    case 4
      offset = varargin(1);
      step = varargin(2);
    otherwise
      error('N:calcbins:badinput','Bad input to N.calcbins! Go to the corner. See help');
  end%switch

  
  intobin = floor((spikearray + ((width / 2) - offset)) / width);
  
  
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