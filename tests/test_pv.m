% TEST_PV  Acceptance tests for tsm_pv (Person B).
% Run from the project root. Needs audio/synth_tone.wav  (run make_test_signals once).
% Set lockModes = false while only the basic version exists; add true after locking is written.

fs = 22050;
if ~exist(fullfile('audio','synth_tone.wav'), 'file'), make_test_signals(); end
x  = load_audio(fullfile('audio','synth_tone.wav'));
lockModes = [false true];                      % later: [false true]

for lock = lockModes
    fprintf('\n=== lock = %d ===\n', lock);

    % 1. length and shape
    for alpha = [0.5 1 1.2 1.8 3]
        y  = tsm_pv(x, alpha, 4096, lock);
        ok = iscolumn(y) && length(y) == round(alpha*length(x));
        fprintf('alpha %.1f: length %d (expected %d)  %s\n', ...
                alpha, length(y), round(alpha*length(x)), passfail(ok));
    end

    % 2. identity at alpha = 1  (ignore first/last 4096 samples)
    y   = tsm_pv(x, 1, 4096, lock);
    err = max(abs(y(4097:end-4096) - x(4097:end-4096)));
    fprintf('identity error = %.3g  %s  (expect ~1e-12)\n', err, passfail(err < 1e-6));

    % 3. pitch preservation: first note is 440 Hz, lasts 2 s in the input
    for alpha = [1.2 1.8 3]
        y   = tsm_pv(x, alpha, 4096, lock);
        seg = y(round(0.5*fs) + (1:fs));                       % 1 s inside the first note
        win = 0.5 - 0.5*cos(2*pi*(0:length(seg)-1)'/length(seg));
        S   = abs(fft(seg .* win));
        [~, i] = max(S(1:floor(length(seg)/2)));
        f   = (i-1) * fs / length(seg);
        fprintf('alpha %.1f: dominant frequency %.1f Hz  %s  (expect 440)\n', ...
                alpha, f, passfail(abs(f-440) <= 2));
    end
end

function s = passfail(ok)
    if ok, s = 'PASS'; else, s = 'FAIL'; end
end
