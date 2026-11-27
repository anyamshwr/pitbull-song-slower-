% SETUP  Run once per MATLAB session from the repository root.
% Adds the code folders to the path. All scripts read audio/ and write out/
% relative to this folder, so keep it as MATLAB's current folder.
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'src'), fullfile(root, 'tests'), fullfile(root, 'experiments'));
cd(root);
fprintf('HP-TSM paths set. Try:  test_integration   or   run_pipeline\n');
