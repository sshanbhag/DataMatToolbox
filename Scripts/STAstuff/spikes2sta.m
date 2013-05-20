function varargout = spikes2sta(spikes, varargin)
%------------------------------------------------------------------------
% [input, category, site] = spikes2sta(spikes)
%------------------------------------------------------------------------
% 
% Takes cell array of spikes and converts to format that can be used
% with the Spike Train Analysis Toolbox (STAtoolkit)
% 
%------------------------------------------------------------------------
% Input Arguments:
%
%	spikes		{# stimuli or conditions, # reps per stim, # sites} cell array
% 
%		spikes array (cell array of spiketimes) is assumed to be organized as:
% 
%			[stim1, rep1]	[stim1, rep2]	...	[stim1, rep P]
%		 	[stim2, rep1]
%			.
%			.
%			.
%			[stim M, rep1]				.	.	.		[stim M, rep P]
% 	
%		where 
%		 	M = # of stimuli / stimulus categories / variable categories
%		 	P = # of reps per stimulus
%			N = # of simultaneously recorded sites (usually 1)
% 
% 	Options:
% 	
% 		'stimlabels'		{M, 1} cell array of strings, stimulus names
% 		'sitenames'			{N, 1} cell array of strings, site names
% 		'sitetags'			{'episodic'} (default) or {'continuous'}
%
%------------------------------------------------------------------------
% 
% Output Arguments:
% 	input		struct
%------------------------------------------------------------------------
%
%------------------------------------------------------------------------
% input struct
%------------------------------------------------------------------------
% Name			Description									Type		Size	
% M				Number of stimulus classes				int32		1	
%					in experiment
% N				Number of sites in experiment			int32		1	
% categories	Array of structures with				category	M
%					categorized response data
% sites			Array of structures with recording	site		N
%					site information
%------------------------------------------------------------------------
%
%------------------------------------------------------------------------
% site struct elements:
%------------------------------------------------------------------------
% Holds information about the physical location(s) from which the data was
% obtained. In our framework, a data set should only contain multiple sites only
% if the sites were recorded simultaneously. Data from sites that were recorded
% sequentially should appear in separate input data structures.
%------------------------------------------------------------------------
% Name					Description									Type				Size	
% label					Recording site label						cell array		1	
% recording_tag		Can be either continuous				cell array		1	
%							or episodic
% time_scale			The scale factor required to			double			1	
%							convert the time mesurements 
%							in seconds
% time_resolution		The temporal resolution of the		double			1	
%							recording (prior to scaling 
%							by time_scale)					
% si_unit				The pluralized international			char				variable	
%							system of units (SI) base or 
%							derived unit of the sampled data
% si_prefix				The international system of units	double			1
%							(SI) prefix (a power of ten)		
%------------------------------------------------------------------------
% label
%     A text label that describes the site.
% recording_tag
%     Can be either continuous or episodic.
% time_scale
%     The scale factor required to convert the time measurements in seconds. If
%     time measurements were reported in milliseconds, then time_scale=0.001.
% time_resolution
%     The temporal resolution of the recording (prior to scaling by time_scale).
% si_unit
%     For continuous data, the pluralized international system of units (SI)
%     base or derived unit of the sampled data (e.g., amperes, volts). This
%     field is not required for episodic data.
% si_prefix
%     For continuous data, the international system of units (SI) prefix, which
%     is a power of ten (e.g., 10-3). Named prefixes are not allowed. This field
%     is not required for episodic data.
%------------------------------------------------------------------------
%
%------------------------------------------------------------------------
% category
%------------------------------------------------------------------------
% The basis for the analysis is to determine how much information an ensemble of
% spike trains convey about some discrete stimulus attribute. The spike trains
% are partitioned into categories according to their stimulus attributes.
%------------------------------------------------------------------------
% Name		Description									Type				Size	
% label		Stimulus class label						cell array		1	
% P			Number of response trials				int32				1	
% trials		Array of structures containing		trial				P×input.N
%				information about trials		
%------------------------------------------------------------------------
%
%------------------------------------------------------------------------
% trial
%------------------------------------------------------------------------
% Holds information about single response trial.
%------------------------------------------------------------------------
% Name				Description									Type			Size	
% start_time		Start time for the trial				double		1	
% end_time			End time for the trial					double		1	
% Q					Number of data points in the list	int32			1	
% list				An array with the data					double		Q
%------------------------------------------------------------------------
% start_time
%     The start time of the recording. This time is scaled by the site's
%     time_scale to give the start time in seconds.
% end_time
%     The end time of the recording. This time is scaled by the site's
%     time_scale to give the end time in seconds.
%------------------------------------------------------------------------
%
%------------------------------------------------------------------------
% See also: 
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% Sharad J. Shanbhag
% sshanbhag@neomed.edu
%------------------------------------------------------------------------
% Created: 3 April, 2013 (SJS)
%
% Revisions:
%	11 Apr 2013 (SJS): changed handling of start_time and end_time
%------------------------------------------------------------------------
% TO DO:
%------------------------------------------------------------------------

%------------------------------------------------------------------------
% get some values from examining the spikes matrix
%------------------------------------------------------------------------
[M, P, N] = size(spikes);

%------------------------------------------------------------------------
%------------------------------------------------------------------------
sitenames = {'default'};
sitetags = {'episodic'};
time_scale = 1e-6;
time_resolution = 1e6;
si_unit = 'none';
si_prefix = 1;
start_time = zeros(M, P, N);
end_time = (1/time_scale).*ones(M, P, N);
	
%------------------------------------------------------------------------
% parse input options
%------------------------------------------------------------------------
nargs = length(varargin);
if nargs
	argn = 1;
	while argn <= nargs
		switch(upper(varargin{argn}))
			case 'STIMLABELS'
				stimlabels = varargin{argn + 1};
				argn = argn + 2;
			case 'SITENAMES'
				sitenames = varargin{argn + 1};
				argn = argn + 2;
			case	'SITETAGS'
				sitetags = varargin{argn + 1};
				argn = argn + 2;
			case 'TIMESCALE'
				time_scale = varargin{argn + 1};
				argn = argn + 2;
			case 'TIMERESOLUTION'
				time_resolution = varargin{argn + 1};
				argn = argn + 2;
			case 'SI_UNIT'
				si_unit = varargin{argn + 1};
				argn = argn + 2;
			case 'SI_PREFIX'
				si_prefix = varargin{argn + 1};
				argn = argn + 2;
			case 'START_TIME'
				start_time = varargin{argn + 1};
				argn = argn + 2;
			case 'END_TIME'
				end_time  = varargin{argn + 1};
				argn = argn + 2;
			otherwise
				error('%s: unknown option %s\n', mfilename, varargin{argn});
		end	% END switch
	end	% END while argn
end	% END if nargs

if ~exist('stimlabels', 'var')
	fprintf('%s: Using indices as stimulus labels...\n', mfilename);
	stimlabels = sprintf('%d', 1:M);
end

if max(size(start_time)) == 1
	start_time = start_time .* ones(M, P, N);
end
if max(size(end_time)) == 1
	end_time = end_time .* ones(M, P, N);
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% site struct elements:
%------------------------------------------------------------------------
% Holds information about the physical location(s) from which the data was
% obtained. In our framework, a data set should only contain multiple sites only
% if the sites were recorded simultaneously. Data from sites that were recorded
% sequentially should appear in separate input data structures.
%------------------------------------------------------------------------
% Name					Description									Type				Size	
% label					Recording site label						cell array		1	
% recording_tag		Can be either continuous				cell array		1	
%							or episodic
% time_scale			The scale factor required to			double			1	
%							convert the time mesurements 
%							in seconds
% time_resolution		The temporal resolution of the		double			1	
%							recording (prior to scaling 
%							by time_scale)					
% si_unit				The pluralized international			char				variable	
%							system of units (SI) base or 
%							derived unit of the sampled data
% si_prefix				The international system of units	double			1
%							(SI) prefix (a power of ten)		
%------------------------------------------------------------------------
% label
%     A text label that describes the site.
% recording_tag
%     Can be either continuous or episodic.
% time_scale
%     The scale factor required to convert the time mesurements in seconds. If
%     time measurements were reported in milliseconds, then time_scale=0.001.
% time_resolution
%     The temporal resolution of the recording (prior to scaling by time_scale).
% si_unit
%     For continuous data, the pluralized international system of units (SI)
%     base or derived unit of the sampled data (e.g., amperes, volts). This
%     field is not required for episodic data.
% si_prefix
%     For continuous data, the international system of units (SI) prefix, which
%     is a power of ten (e.g., 10-3). Named prefixes are not allowed. This field
%     is not required for episodic data.
%------------------------------------------------------------------------
site = repmat(	struct(	'label', {}, ...
								'recording_tag', {}, ...
								'time_scale', [], ...
								'time_resolution', [], ...
								'si_unit', [], ...
								'si_prefix', ''	), ...
					N, 1);
for n = 1:N
	site(n).label = {sitenames{n}};
	site(n).recording_tag = {sitetags{n}};
	site(n).time_scale = time_scale(n);
	site(n).time_resolution = time_resolution(n);
	site(n).si_unit = si_unit;
	site(n).si_prefix = si_prefix;
end
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% category
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% The basis for the analysis is to determine how much information an ensemble of
% spike trains convey about some discrete stimulus attribute. The spike trains
% are partitioned into categories according to their stimulus attributes.
%------------------------------------------------------------------------
% Name		Description									Type				Size	
% label		Stimulus class label						cell array		1	
% P			Number of response trials				int32				1	
% trials		Array of structures containing		trial				P×input.N
%				information about trials		
%------------------------------------------------------------------------
category = repmat(struct('label', '', 'P', int32(P), 'trials', []), M, 1);
for m = 1:M
	category(m).label = {stimlabels(m)};
	category(m).trials = repmat(	struct(	'start_time', [], ...
														'end_time', [], ...
														'Q', [], ...
														'list', [] ), ...
											P, 1);
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% trial
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% Holds information about single response trial.
%------------------------------------------------------------------------
% Name				Description									Type			Size	
% start_time		Start time for the trial				double		1	
% end_time			End time for the trial					double		1	
% Q					Number of data points in the list	int32			1	
% list				An array with the data					double		Q
%------------------------------------------------------------------------
% start_time
%     The start time of the recording. This time is scaled by the site's
%     time_scale to give the start time in seconds.
% end_time
%     The end time of the recording. This time is scaled by the site's
%     time_scale to give the end time in seconds.
%------------------------------------------------------------------------
%------------------------------------------------------------------------
n = 1;
for m = 1:M
	for p = 1:P
		category(m).trials(p).start_time = start_time(m, p, n);
		category(m).trials(p).end_time = end_time(m, p, n);
		category(m).trials(p).Q = int32(length(spikes{m, p}));
		category(m).trials(p).list = spikes{m, p};
	end
end

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% input struct
%------------------------------------------------------------------------
% Name			Description									Type		Size	
% M				Number of stimulus classes				int32		1	
%					in experiment
% N				Number of sites in experiment			int32		1	
% categories	Array of structures with				category	M
%					categorized response data
% sites			Array of structures with recording	site		N
%					site information
%------------------------------------------------------------------------
input = struct('M', int32(M), 'N', int32(N), 'categories', category, 'sites', site);

if any(nargout == (0:3))
	varargout{1} = input;
end

if any(nargout == (2:3))
	varargout{2} = site;
end

if any(nargout == 3)
	varargout{3} = category;
end
