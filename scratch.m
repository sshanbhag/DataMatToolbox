function varargout = scratch(varargin)
	nargout
	
	varargout{1} = 'arg1';

end





%{ 
load F:\Work\Data\DataWave\DDFtest\d.mat

for stim = 1:d.Stimuli.N
	fprintf('%d\t\t%s\n', stim, d.Stimuli.S{stim, 2}.Filename)
	
	
end
return
% check channel for this stimulus
c = d.Stimuli.Channel(stim)

switch c
	
	case 'L'
		
		
	case 'R'
		
		[m, mlist] = d.Stimuli.S{stim, 2}.match(d.Stimuli.S(:, 2))
		
	case 'B'
		
end

%}


%{
for stim = 1:d.Stimuli.N
% for stim = 1:1
	fprintf('%d sweeps for stimulus %d\n', d.Stimuli.Nsweeps(stim), stim);
	Sweepstart{stim} = zeros(d.Stimuli.Nsweeps(stim), 1);
	Sweepend{stim} = zeros(d.Stimuli.Nsweeps(stim), 1);
	PreSweep{stim} = zeros(d.Stimuli.Nsweeps(stim), 1);
	PostSweep{stim} = zeros(d.Stimuli.Nsweeps(stim), 1);

	for sweep = 1:d.Stimuli.Nsweeps(stim)
		mindx = d.Stimuli.MarkerList{stim}(sweep);
		fprintf('\t%.0f\t\t%.0f\t%s\t%.0f\n', d.Markers(mindx).OutputTimestampL, ...
			d.Markers(mindx).OutputTimestampR, d.Markers(mindx).WavFilenameR, ...
			d.Markers(mindx).AttenuationR);
		
		if mindx == 1
			preindx = 0;
		else
			preindx = mindx - 1;
		end
		if mindx == Nmarkers
			postindx = Nmarkers;
		else
			postindx = mindx + 1;
		end
		
		if d.Stimuli.Channel{stim} == 'L'
			Sweepstart{stim}(sweep) = d.Markers(mindx).OutputTimestampL;
			Sweepend{stim}(sweep) = d.Markers(postindx).OutputTimestampL;
		elseif d.Stimuli.Channel{stim} == 'R'
			Sweepstart{stim}(sweep) = d.Markers(mindx).OutputTimestampR;
			Sweepend{stim}(sweep) = d.Markers(postindx).OutputTimestampR;			
		elseif d.Stimuli.Channel{stim} == 'B'
			% use min value of l and r
			Sweepstart{stim}(sweep) = min([d.Markers(mindx).OutputTimestampL ...
														d.Markers(mindx).OutputTimestampR]);
			Sweepend{stim}(sweep) = min([d.Markers(postindx).OutputTimestampL ...
														d.Markers(postindx).OutputTimestampR]);
		else
			error('%s: bad channel value %s', mfilename, d.Stimuli.Channel{stim});
		end
		% correct final sweep time
		if postindx == Nmarkers
			Sweepend{stim}(sweep) = Sweepstart{stim}(sweep) + SWEEPDUR; 
		end
		PreSweep{stim}(sweep) = preindx;
		PostSweep{stim}(sweep) = postindx;

	end	% END sweep
	
	
	
end	% END stim
%}