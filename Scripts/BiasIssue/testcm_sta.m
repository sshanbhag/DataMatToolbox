
clear all
close all


%------------------------------------------------------------------------
%------------------------------------------------------------------------
% set path to toolboxes
%------------------------------------------------------------------------
%------------------------------------------------------------------------
% STAkitMAC
tmp = which('metric');
if isempty(tmp)
	fprintf('%s: setting up STAtoolkit paths...', mfilename)
	addpath(genpath('/Users/sshanbhag/Work/Code/Matlab/stable/Toolbox/STAkitMAC'))
	fprintf('... done \n\n')
end
clear tmp;

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% reset random # gen
%------------------------------------------------------------------------
%------------------------------------------------------------------------
RandObj = RandStream.getDefaultStream;
RandObj.reset;

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% set values
%------------------------------------------------------------------------
%------------------------------------------------------------------------
Stimuli = 2:15;
Reps = 5:5:100;
Iterations = 100;

% options for STA toolkit
opts.unoccupied_bins_strategy = -1;
opts.entropy_estimation_method = {'plugin', 'tpmc'};
opts.possible_words = 'recommended';

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% pre allocate arrays
%------------------------------------------------------------------------
%------------------------------------------------------------------------
H.mean = zeros(length(Stimuli), length(Reps));
H.sd = zeros(length(Stimuli), length(Reps));
Hplugin.mean = zeros(length(Stimuli), length(Reps));
Hplugin.sd = zeros(length(Stimuli), length(Reps));
Htpmc.mean = zeros(length(Stimuli), length(Reps));
Htpmc.sd = zeros(length(Stimuli), length(Reps));

%------------------------------------------------------------------------
%------------------------------------------------------------------------
% run!
%------------------------------------------------------------------------
%------------------------------------------------------------------------

%------------------------------------------
% loop through nstimuli list
%------------------------------------------
for s = 1:length(Stimuli)
	Nstimuli = Stimuli(s);
	%------------------------------------------
	% loop through nreps list
	%------------------------------------------
	for r = 1:length(Reps)
		% get current # of reps
		Nreps = Reps(r);
		%------------------------------------------
		% allocate arrays to store information values
		%------------------------------------------
		h = zeros(1, Iterations);
		hplugin = zeros(1, Iterations);
		htpmc = zeros(1, Iterations);
		%------------------------------------------
		% loop for # of iterations
		%------------------------------------------
		for n = 1:Iterations
			%------------------------------------------------------------------
			% create matrix of counts for confusion matrix
			% confusion matrix has # stimuli along rows, test along column
			%------------------------------------------------------------------
			% need to generate using monte carlo method (to ensure that 
			% votes along each row add up to consistent value)
			%------------------------------------------------------------------
			CM = zeros(Nstimuli);
			% loop through rows (responses)
			for cm_s = 1:Nstimuli
				% loop for number of trials
				for cm_r = 1:Nreps
					rbin = randi(Nstimuli, 1);
					CM(cm_s, rbin) = CM(cm_s, rbin) + 1;
				end
			end
			%------------------------------------------------------------------
			% compute information from confusion matrix, store in H.
			%------------------------------------------------------------------
			h(n) = cminfo(CM);
			%------------------------------------------------------------------
			% compute information using STA toolkit
			%------------------------------------------------------------------
			% first, convert CM to 2d hist
			[table, opts] = matrix2hist2d(CM, opts);
			% then use info2d to get information measures
			[table, opts] = info2d(table, opts);
			% store values
			hplugin(n) = table.information(1).value;
			htpmc(n) = table.information(2).value;
			
		end	% END Iterations Loop
		% compute means and standard deviations
		H.mean(s, r) = mean(h);
		H.sd(s, r) = std(h);
		Hplugin.mean(s, r) = mean(hplugin);
		Hplugin.sd(s, r) = std(hplugin);
		Htpmc.mean(s, r) = mean(htpmc);
		Htpmc.sd(s, r) = std(htpmc);
		
	end	% END Nresp Loop
end	% END Nstimuli Loop

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% save results matrices
%------------------------------------------------------------------------
%------------------------------------------------------------------------

save('testcm_sta_results.mat', '-MAT')

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% display results matrices
%------------------------------------------------------------------------
%------------------------------------------------------------------------
display_matrix(H.mean, '%.2f');
display_matrix(H.sd, '%.2f');
display_matrix(Hplugin.mean, '%.2f');
display_matrix(Hplugin.sd, '%.2f');
display_matrix(Htpmc.mean, '%.2f');
display_matrix(Htpmc.sd, '%.2f');

%------------------------------------------------------------------------
%------------------------------------------------------------------------
%% plot results
%------------------------------------------------------------------------
%------------------------------------------------------------------------

corder = [	0.00	0.00	1.00	;	...
				0.00	1.00	0.00	;	...
				1.00	0.00	0.00	;	...
				0.00	0.25	0.75	;	...
				0.00	0.75	0.25	;	...
				0.25	0.75	0.00	;	...
				0.75	0.25	0.00	;	...
				0.00	0.50	0.50	;	...
				0.50	0.00	0.50	;	...
				1.00	0.00	1.00	;	...
				1.00	1.00	0.00	;	...
				0.00	1.00	1.00	;	...
				0.00	0.00	0.00	;	...
			];

figure(1)
set(gcf, 'DefaultAxesColorOrder', corder);
subplot(211)
plot(Reps, H.mean', '.-')
legend(num2str(Stimuli'))
legend(gca, 'Location', 'Best')
legend boxoff
ylabel('H avg (bits)')
title({'Random Confusion Matrix', 'Victor, Purpura Method'})

subplot(212)
plot(Reps, H.sd', '.-');
xlabel('number of trials per stimulus')
ylabel('H std. dev. (bits)')
title('colors indicate # of stimuli')

figure(2)
set(gcf, 'DefaultAxesColorOrder', corder);
subplot(211)
plot(Reps, Hplugin.mean', '.-')
legend(num2str(Stimuli'))
legend(gca, 'Location', 'Best')
legend boxoff
ylabel('H avg (bits)')
title({'Random Confusion Matrix', 'STA plugin Method'})

subplot(212)
plot(Reps, Hplugin.sd', '.-');
xlabel('number of trials per stimulus')
ylabel('H std. dev. (bits)')
title('colors indicate # of stimuli')

figure(3)
set(gcf, 'DefaultAxesColorOrder', corder);
subplot(211)
plot(Reps, Htpmc.mean', '.-')
legend(num2str(Stimuli'))
legend(gca, 'Location', 'Best')
legend boxoff
ylabel('H avg (bits)')
title({'Random Confusion Matrix', 'STA tpmc Method'})

subplot(212)
plot(Reps, Htpmc.sd', '.-');
xlabel('number of trials per stimulus')
ylabel('H std. dev. (bits)')
title('colors indicate # of stimuli')

%%
figure(1)
subplot(211); ylim([-0.2 1.8]); legend(gca, 'Location', 'NorthEast'); legend boxoff
subplot(212); ylim([0 0.2]);
figure(2)
subplot(211); ylim([-0.2 1.8]); legend(gca, 'Location', 'NorthEast'); legend boxoff
subplot(212); ylim([0 0.2]);
figure(3)
subplot(211); ylim([-0.2 1.8]); legend(gca, 'Location', 'NorthEast'); legend boxoff
subplot(212); ylim([0 0.2]);
