
ngroups = length(d.Stimuli.GroupList);
for g = 1:ngroups
	figure
	% get all spikes for stimulus group
	Sindx = d.Stimuli.GroupList{g}

	spikes = d.getSpikesForStimulusGroup(g, 1);
	for n = 1:length(Sindx)
		subplot(length(Sindx), 1, n)
		rasterplot(spikes{n}, [0 1200*1000])
		s = Sindx(n);
		tstr1 = fullfile(d.Stimuli.S{s, 2}.Filepath, d.Stimuli.S{s, 2}.Filename);
		tstr2 = sprintf('Atten = %d', d.Stimuli.S{s, 2}.Attenuation);
		title({tstr1, tstr2}, 'Interpreter', 'none');
	end

end