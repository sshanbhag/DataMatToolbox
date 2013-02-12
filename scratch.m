
S = d.getSpikesForProbe(1);
ngroups = length(S);
allspikes = cell(3, ngroups);
% loop through groups
for g = 1:ngroups
	figure(g)
	% get Stimulus List indices for this group
	Sindx = d.Stimuli.GroupList{g};
	% loop through stim indices
	for n = 1:length(Sindx)
		subplot(length(Sindx), 1, n)
		% convert spiketimes to milliseconds
		spikes = cell(length(S(g).spikes{n}), 1);
 		for t = 1:length(S(g).spikes{n})
 			spikes{t} = 0.001*S(g).spikes{n}{t};
 		end
		rasterplot(spikes, [0 1200])
		s = Sindx(n);
		tstr1 = fullfile(d.Stimuli.S{s, 2}.Filepath, d.Stimuli.S{s, 2}.Filename);
		tstr2 = sprintf('Atten = %d', d.Stimuli.S{s, 2}.Attenuation);
		title({tstr1, tstr2}, 'Interpreter', 'none');
		allspikes{n, g} = spikes;
	end

end

%% test stimulus group sort algorithm
%{
S = { 'a', 'b', 'c', 'd', 'e', 'a', 'b', 'c', 'd', 'e'};
N = length(S);

GroupList = {};
runFlag = 1;
nUnique = 0;spik
% search all indices initially
masterList = 1:N;
searchList = 1:N;
%----------------------------------
% loop
%----------------------------------
while runFlag
	% pick the element to compare against the full list
	stim = searchList(1);
	fprintf('%d:\t', stim);
	% compare
	comp = logical(strcmpi(S{stim}, S(searchList)));

	% store unique indices
	nUnique = nUnique + 1;
	searchList(comp)
	GroupList{nUnique} = masterList(comp); %#ok<*SAGROW>
	% eliminate them from the list to search
	searchList = searchList(~(comp));
	% and shrink the unsearched masterList as well
	masterList = masterList(~(comp));
	% Check if we're done
	if isempty(searchList)
		runFlag = 0;
	end
end
%}
