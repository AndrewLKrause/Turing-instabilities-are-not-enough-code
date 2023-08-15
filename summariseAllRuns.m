files = dir("*.mat");
for i = 1 : length(files)
	load(files(i).name);
	summariseRuns;
	fprintf("\n\n");
end
