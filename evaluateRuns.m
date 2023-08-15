mask = logical(TuringUnstable);
numNoPatterns = sum(Patterns <= 1e-5);
numWhereTuringNotEnough = sum(Patterns(mask) <= 1e-5);

disp([modelName, ' in ', num2str(dims), 'D.'])
disp([num2str(NumRuns), ' runs in total. Of those, ', num2str(numNoPatterns), ' ended up with no patterns.']);
disp([num2str(sum(mask)), ' of the total were Turing unstable. ', 'Of those, ', num2str(numWhereTuringNotEnough), ' ended up with no patterns.'])
disp(['So, Turing was not enough in ', num2str(100*numWhereTuringNotEnough/sum(mask),'%.2f'), '% of cases.'])
