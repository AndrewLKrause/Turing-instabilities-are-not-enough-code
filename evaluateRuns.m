mask = logical(TuringUnstable);
numWhereTuringNotEnough = sum(Patterns(mask) <= 1e-5);

disp(modelName)
disp([num2str(NumRuns), ' runs. ', num2str(sum(mask)), ' were Turing unstable. ', 'Of those, ', num2str(numWhereTuringNotEnough), ' resulted in no patterns.'])
disp(['So, Turing was not enough in ', num2str(100*numWhereTuringNotEnough/sum(mask),'%.2f'), '% of cases.'])