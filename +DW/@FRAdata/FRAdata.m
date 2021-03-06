%-----------------------------------------------------------------------------
% FRAdata.m < Data
%-----------------------------------------------------------------------------
% DataMat Toolbox
% DW Package
% Class Definition
%-----------------------------------------------------------------------------
% FRAdata implements analysis specific to Frequency Response Area data
% and is a subclass of the Data class in the DW (datawave) package
%-----------------------------------------------------------------------------
% See also: Data (class)
%-----------------------------------------------------------------------------

%-----------------------------------------------------------------------------
% initial coding and design:
%	Tony Slagle
%	tonyslagle@gmail.com
%
% Continuing development: 
%	Sharad J. Shanbhag
%	sshanbhag@neomed.edu
%-----------------------------------------------------------------------------
% Created: 19 February, 2013 as subclass of Data (SJS)
%
% Revisions:
%-----------------------------------------------------------------------------

%*****************************************************************************
%*****************************************************************************
%*****************************************************************************
% class definition
%*****************************************************************************
%*****************************************************************************
%*****************************************************************************
classdef FRAdata < DW.Data
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define protected properties
	%------------------------------------------------------------------------
	properties (SetAccess = protected)
		Frequencies
		AttenLevels		
		Nfreqs
		sortedFreqs
		sortFreqsX
		sortedAtten
		sortAttX
		attcount
		Natten
		SpikeTimes
		SpikeCount
		MeanCount
		StdDevCount
	end
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define other properties
	%------------------------------------------------------------------------
	%
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
  
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	% Define methods
	%------------------------------------------------------------------------
	methods
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		% Constructor
		%------------------------------------------------------------------------
		% Initializes the object
		%------------------------------------------------------------------------
		function obj = FRAdata(varargin)
		%---------------------------------------------------------------------	
		%	FRAdata(<fileName>) 
		%	Constructor method
		%	opens file called fileName (char) or opens
		%  a dialog box to get a filename if the fileName provided does not exist.
		%---------------------------------------------------------------------	

			%--------------------------------------------------------
			% call superclass constructor
			%--------------------------------------------------------
			obj = obj@DW.Data(varargin);
			%--------------------------------------------------------		
			% parse input and verify
			%--------------------------------------------------------
			if nargin == 0
				return
			end
		end	% END FRAdata CONSTRUCTOR
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		% Data input/reading/plotting methods
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		
		%------------------------------------------------------------------------
		function plotFRA(obj, probenum, unitnum, frawin, varargin)
		%------------------------------------------------------------------------
		%	FRAdata.plotFRA(probenum, unitnum, frawin)
		%------------------------------------------------------------------------
			
			%------------------------------------------------
			% compute FRA for desired unit and window
			%------------------------------------------------
			obj.computeFRA(probenum, unitnum, frawin);
			%------------------------------------------------
			% plot as color patch
			%------------------------------------------------
			% create figure and subplot
			figure
			subplot(211);
			% plot as color patch.  need to flip sorted atten to get
			% higher atten (lower amplitude tones) at bottom of plot and
			% lower atten (higher amp) at top, per FRA plot convention
			% the y axis labels will be in reverse order, but we'll take care
			% of that later
			xdata = log10(obj.sortedFreqs);
			ydata = fliplr(obj.sortedAtten{1});
			pcolor(xdata, ydata, obj.MeanCount);
			% show color legend
			colorbar
			% deal with labels and title
			xlabel('Log Frequency (kHz)');
			ylabel('Attenuation (dB)');
 			title(	{	obj.fname, ...
							sprintf('Avg Spike Count, [%d-%d] ms window', ...
														frawin(1), frawin(2)), ...
							sprintf('Probe: %d  Unit %d', ...
															probenum, unitnum) }, ...
						'Interpreter', 'none');
			% re-do X tick labels to that they're more readable
			xt = get(gca, 'XTick');
			xtl = cell(length(xt), 1);
			get(gca, 'XTickLabel');
			for n = 1:length(xt);
				xtl{n} = sprintf('%.0f', 0.001 * 10^xt(n));
			end
			set(gca, 'XTickLabel', xtl);
			% correct the Y tick labels, as promised
 			set(gca, 'YTickLabel', flipud(get(gca, 'YTickLabel')))
			%------------------------------------------------
			% create subplot and plot waterfall 
			% (in fashion similar to color patch)
			%------------------------------------------------
			subplot(212);
			waterfall(xdata, ydata, obj.MeanCount);
			% labels, again deal with log labels
			xlabel('Log Frequency (kHz)');
			ylabel('Attenuation (dB)');
			zlabel('Spike Count');
			xt = get(gca, 'XTick');
			xtl = cell(length(xt), 1);
			get(gca, 'XTickLabel');
			for n = 1:length(xt)
				xtl{n} = sprintf('%.0f', 0.001 * 10^xt(n));
			end
			set(gca, 'XTickLabel', xtl);
			set(gca, 'YTickLabel', flipud(get(gca, 'YTickLabel')))
		end	% END plotFRA
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function S = getFRAspikes(obj, varargin)
		%------------------------------------------------------------------------
		%
		%------------------------------------------------------------------------
		%	'probe', probenum
		%		'probe' option will select probenum for plotting
		%	'unit', unitnumber
		%		'unit' selects unit to display
		%	'window', [tstart tend]
		%		species time window for spikes (re: start of sweep) in millisec
		%
		%	out is a struct array (with # of elements in array == # freqs)
		%		out(f).spikes = {# atten vals, 1} cell array of spike times (usec)
		%		out(f).name = char string with spike name
		%		out(f).freq = stimulus frequency (Hz)
		%		out(f).atten = [# atten vals, 1] attenuation values
		%
		%------------------------------------------------------------------------
			
			%------------------------------------------------
			% process inputs
			%------------------------------------------------
			probenum = 1;
			unitnum = 0;
			spwin = [];
			if ~isempty(varargin)
				a = 1;
				while a <= length(varargin)
					switch upper(varargin{a})	
						case 'PROBE'
							probenum = varargin{a+1};
							a = a + 2;
						case 'UNIT'
							unitnum = varargin{a+1};
							a = a + 2;
						case 'WINDOW'
							spwin = varargin{a+1};
							a = a + 2;
						otherwise
							error('%s: unknown option %s', mfilename, varargin{a});
					end
				end				
			end
			if ~between(probenum, 1, length(obj.Probes))
				error('%s: probe must be in range [1:%d]', mfilename, ...
																	length(obj.Probes));
			end			
			%----------------------------------------------------------
			% check if Frequencies and AttenLevels have been found
			%----------------------------------------------------------
			if isempty(obj.Frequencies) || isempty(obj.AttenLevels)
				% if not, find 'em
				obj.findFreqAndAtten;
			end			
			
			%----------------------------------------------------------
			% get the spikes struct for probe and unit, and store it
			%----------------------------------------------------------
			S = obj.getSpikesForProbe(probenum, 'unit', unitnum);
			
			%----------------------------------------------------------
			% get spikes within analysis window
			%----------------------------------------------------------
			% if no values present for window, we're done
			if isempty(spwin)
				return
			end
			% otherwise, loop through S...
			
			for s = 1:length(S)
				% assign attenuations and frequency
				
				% and loop through spikes cell array...
				for c = 1:length(S(s).spikes)
					% and loop through reps
					for r = 1:length(S(s).spikes{c})
						% and select spikes
						tmp = find_valid_timestamps( ...
														S(s).spikes{c}{r}, ...
														1000 * spwin(1), ...
														1000 * spwin(2) );
						S(s).spikes{c}{r} = tmp{1};
					end
				end
			end
			
		end	% END getFRA
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
	
		function S = getSortedFRASpikes(obj, probenum, unitnum, frawin)
		%------------------------------------------------------------------------
		%	FRAdata.getSortedFRASpikes(probenum, unitnum, frawin)
		%------------------------------------------------------------------------
			fprintf('%s: getting spikes for probe %d, unit %d, [%d-%d] window\n', ...
									mfilename, probenum, unitnum, frawin(1), frawin(2));
			%------------------------------------------------
			% sort frequencies and atten
			%------------------------------------------------
			if isempty(obj.sortedFreqs)
				obj.sortFrequenciesAndAtten;
			end
			
			%------------------------------------------------
			% get the spikes
			%------------------------------------------------
			Spikes = obj.getFRAspikes('probe', probenum, 'unit', unitnum, ...
										'window', frawin);

			%------------------------------------------------
			% get the spike times for each sorted freq and 
			% level and count # of spikes
			%------------------------------------------------
			% allocate obj.SpikeTimes and FRAspikecount
			S.SpikeTimes = cell(obj.Natten, obj.Nfreqs);
			S.SpikeCount = cell(obj.Natten, obj.Nfreqs);
			S.MeanCount = zeros(obj.Natten, obj.Nfreqs);
			S.StdDevCount = zeros(obj.Natten, obj.Nfreqs);
			% loop through frequencies (sorted)
			for f = 1:obj.Nfreqs
				fIndx = obj.sortFreqsX(f);
				% loop through attenuation from low to high
				for a = 1:obj.Natten
					aIndx = obj.sortAttX{f}(a);
					% store the spike times in the obj.SpikeTimes cell array
					S.SpikeTimes{a, f} = Spikes(fIndx).spikes{aIndx};
					% add up number of spikes for this freq and atten combination
					% and store each rep in obj.SpikeCount cell matrix
					S.SpikeCount{a, f} = zeros(length(S.SpikeTimes{a, f}), 1);
					for rep = 1:length(S.SpikeTimes{a, f})
						S.SpikeCount{a, f}(rep) = length(S.SpikeTimes{a, f}{rep});
					end	% END rep LOOP
					% compute mean and std dev spike count
					S.MeanCount(a, f) = mean(S.SpikeCount{a, f});
					S.StdDevCount(a, f) = std(S.SpikeCount{a, f});
				end	% END a LOOP
			end	% END f LOOP			
			%------------------------------------------------
			% compute mean and std dev spike count
			%------------------------------------------------
			% loop through frequencies (sorted)
			for f = 1:obj.Nfreqs
				% loop through attenuation from low to high
				for a = 1:obj.Natten
					S.MeanCount(a, f) = mean(S.SpikeCount{a, f});
					S.StdDevCount(a, f) = std(S.SpikeCount{a, f});
				end	% END a LOOP
			end	% END f LOOP
			S.freqs = obj.sortedFreqs;
			S.atten = obj.sortedAtten;
		end	% end METHOD
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		% Analysis
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		
		
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function sortFrequenciesAndAtten(obj)
		%------------------------------------------------------------------------
		%	FRAdata.sortFrequenciesAndAtten()
		%------------------------------------------------------------------------
		% sort frequencies and atten
		%------------------------------------------------------------------------
			fprintf('%s: sorting FRA freqs and atten...', mfilename);
			
			%------------------------------------------------
			% find frequencies and atten if necessary
			%------------------------------------------------
			if isempty(obj.Frequencies) || (obj.Nfreqs == 0)
				obj.findFreqAndAtten;
			end
				%------------------------------------------------
			% sort frequencies and atten
			%------------------------------------------------
			% count number of different tone frequencies
			obj.Nfreqs = length(obj.Frequencies);
			% sort freqs from low to high, keeping indices
			[obj.sortedFreqs, obj.sortFreqsX] = sort(obj.Frequencies);
			% sort the atten levels - a little more complicated due to 
			% there being a list of atten vals for each individual frequency
			obj.sortedAtten = cell(obj.Nfreqs, 1);
			obj.sortAttX = cell(obj.Nfreqs, 1);
			obj.attcount = zeros(obj.Nfreqs, 1);
			% loop through the frequencies
			for f = 1:obj.Nfreqs
				[obj.sortedAtten{f}, obj.sortAttX{f}] = sort(obj.AttenLevels{f});
				obj.attcount(f) = length(obj.sortedAtten{f});
			end
			% find max # of atten levels
			obj.Natten = max(obj.attcount);	
		
		end	% END sorfFrequenciesAndAtten(obj)
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------

		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		function computeFRA(obj, probenum, unitnum, frawin)
		%------------------------------------------------------------------------
		%	FRAdata.computeFRA(probenum, unitnum, frawin)
		%------------------------------------------------------------------------
			fprintf('%s: computing FRA for probe %d, unit %d, [%d-%d] window\n', ...
									mfilename, probenum, unitnum, frawin(1), frawin(2));
			%------------------------------------------------
			% sort frequencies and atten
			%------------------------------------------------
			if isempty(obj.sortedFreqs) || (obj.Nfreqs == 0)
				obj.sortFrequenciesAndAtten;
			end
			
			%------------------------------------------------
			% get the spikes
			%------------------------------------------------
			Spikes = obj.getFRAspikes('probe', probenum, 'unit', unitnum, ...
										'window', frawin);

			%------------------------------------------------
			% get the spike times for each sorted freq and 
			% level and count # of spikes
			%------------------------------------------------
			% allocate obj.SpikeTimes and FRAspikecount
			obj.SpikeTimes = cell(obj.Natten, obj.Nfreqs);
			obj.SpikeCount = cell(obj.Natten, obj.Nfreqs);
			obj.MeanCount = zeros(obj.Natten, obj.Nfreqs);
			obj.StdDevCount = zeros(obj.Natten, obj.Nfreqs);
			% loop through frequencies (sorted)
			for f = 1:obj.Nfreqs
				fIndx = obj.sortFreqsX(f);
				% loop through attenuation from low to high
				for a = 1:obj.Natten
					aIndx = obj.sortAttX{f}(a);
					% store the spike times in the obj.SpikeTimes cell array
					obj.SpikeTimes{a, f} = Spikes(fIndx).spikes{aIndx};
					% add up number of spikes for this freq and atten combination
					% and store each rep in obj.SpikeCount cell matrix
					obj.SpikeCount{a, f} = zeros(length(obj.SpikeTimes{a, f}), 1);
					for rep = 1:length(obj.SpikeTimes{a, f})
						obj.SpikeCount{a, f}(rep) = length(obj.SpikeTimes{a, f}{rep});
					end	% END rep LOOP
					% compute mean and std dev spike count
					obj.MeanCount(a, f) = mean(obj.SpikeCount{a, f});
					obj.StdDevCount(a, f) = std(obj.SpikeCount{a, f});
				end	% END a LOOP
			end	% END f LOOP			
			%------------------------------------------------
			% compute mean and std dev spike count
			%------------------------------------------------
			% loop through frequencies (sorted)
			for f = 1:obj.Nfreqs
				% loop through attenuation from low to high
				for a = 1:obj.Natten
					obj.MeanCount(a, f) = mean(obj.SpikeCount{a, f});
					obj.StdDevCount(a, f) = std(obj.SpikeCount{a, f});
				end	% END a LOOP
			end	% END f LOOP
	
		end	% END computeFRA
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
	
		%------------------------------------------------------------------------
		function varargout = findFreqAndAtten(obj)
		%------------------------------------------------------------------------
		% ��� assumption: stimuli delivered from only 1 speaker !!!
		%------------------------------------------------------------------------
		
			if isempty(obj.Stimuli)
				error('%s: StimulusList is not initialized', mfilename);
			end
			
			% get channels for data
			C = obj.Stimuli.getChannelAsNum;
			obj.Frequencies = zeros(length(obj.Stimuli.GroupList), 1);
			obj.AttenLevels = cell(length(obj.Stimuli.GroupList), 1);
			% loop through groups
			for g = 1:length(obj.Stimuli.GroupList)
				% get indices into C (and S) for this group
				sind = obj.Stimuli.GroupList{g};
				freqs = zeros(size(sind));
				atts = zeros(size(sind));
				fprintf('Group %d:\n', g);
				% loop through indices
				for s = 1:length(sind)
					% get the stim object
					sobj = obj.Stimuli.S{sind(s), C(sind(s))};
					fprintf('\t%.0f\t\t%.2f\n', sobj.Freq, sobj.Attenuation);
					freqs(s) = sobj.Freq;
					atts(s) = sobj.Attenuation;
				end
				if any(freqs(1) ~= freqs)
					fprintf('Grouping Error!\n');
				end
				obj.Frequencies(g) = freqs(1);
				obj.AttenLevels{g} = atts;
			end
		
			% assign output
			if nargout
				varargout{1} = C;
			end

		end	% END findFreqAndAtten
		%------------------------------------------------------------------------
		%------------------------------------------------------------------------
		
	end	% End of methods
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
	%------------------------------------------------------------------------
end	% End of classdef
%*****************************************************************************
%*****************************************************************************
%*****************************************************************************

