% TEST_WSOLA  Acceptance tests for tsm_wsola (Person C, baseline).
% Run from the project root:   >> test_wsola

clear; clc;
fs = 22050;
if ~exist(fullfile('audio','synth_tone.wav'), 'file'), make_test_signals(); end
tone = load_audio(fullfile('audio','synth_tone.wav'));
clk  = load_audio(fullfile('audio','synth_clicks.wav'));
mix  = load_audio(fullfile('audio','synth_mix.wav'));
if ~exist('out', 'dir'), mkdir('out'); end

%% Test 1: length and shape ---------------------------------------------
for alpha = [0.5 1 1.2 1.8 3]
    y = tsm_wsola(mix, alpha);
    assert(iscolumn(y) && isreal(y), 'Test 1 FAILED: output must be a real column');
    assert(length(y) == round(alpha*length(mix)), ...
           sprintf('Test 1 FAILED: alpha %.1f gives %d samples', alpha, length(y)));
end
fprintf('Test 1  length / shape for alpha in {0.5,1,1.2,1.8,3} : ok\n');

%% Test 2: identity at alpha = 1 ----------------------------------------
for sig = {tone, mix}
    s   = sig{1};
    y   = tsm_wsola(s, 1);
    err = max(abs(y(1025:end-1024) - s(1025:end-1024)));
    fprintf('Test 2  identity error            : %.2e  (must be < 1e-9)\n', err);
    assert(err < 1e-9, 'Test 2 FAILED: alpha = 1 does not reproduce the input');
end

%% Test 3: tone at 1.8 is smooth (WSOLA vs plain OLA) -------------------
alpha = 1.8;
yw  = tsm_wsola(tone, alpha);
yo  = tsm_ola(tone, alpha, 1024);                  % same frame size, no search
seg = @(y) y(round(0.5*fs) + (1:fs));              % 1 s inside the first note
lw  = tonal_leak(seg(yw), fs);
lo  = tonal_leak(seg(yo), fs);
fprintf('Test 3  energy outside harmonics: WSOLA %.4f   OLA(N=1024) %.4f\n', lw, lo);
assert(lw < lo, 'Test 3 FAILED: WSOLA should be cleaner than OLA on a tone');
audiowrite(fullfile('out','wsola_tone_a1.8.wav'), 0.99*yw/max(abs(yw)), fs);

%% Test 4: clicks at 1.8 stutter (paper Fig. 3f) ------------------------
% Count separate bursts per expected click: a doubled click shows up as two
% peaks in the envelope within the ~1000-sample frame region.
yc   = tsm_wsola(clk, alpha);
tExp = alpha * (0.25 : 0.5 : 4 - 0.3);
nb   = bursts_per_click(yc, tExp, fs);
fprintf('Test 4  bursts per click (1 = clean): %s\n', mat2str(nb));
fprintf('        clicks repeated or doubled: %d of %d\n', sum(nb > 1), numel(nb));
audiowrite(fullfile('out','wsola_clicks_a1.8.wav'), 0.99*yc/max(abs(yc)), fs);
ym = tsm_wsola(mix, alpha);
audiowrite(fullfile('out','wsola_mix_a1.8.wav'), 0.99*ym/max(abs(ym)), fs);
fprintf('        saved out/wsola_*_a1.8.wav for listening\n');

fprintf('\nALL WSOLA TESTS PASSED\n');


% ---------------------------------------------------------------- helpers
function r = tonal_leak(seg, fs)
    n   = length(seg);
    win = 0.5 - 0.5*cos(2*pi*(0:n-1)'/n);
    S   = abs(fft(seg .* win)).^2;
    S   = S(1:floor(n/2));
    f   = (0:floor(n/2)-1)' * fs / n;
    near = false(size(f));
    for h = [440 880 1320], near = near | abs(f - h) <= 5; end
    r = sum(S(~near)) / sum(S);
end

function nb = bursts_per_click(y, tc, fs)
% Envelope = 2 ms moving max of |y|; count rising crossings of 30% of the
% local peak in a 100 ms window starting 10 ms before the expected onset.
    nb = zeros(size(tc));
    for k = 1:numel(tc)
        i0  = max(1, round((tc(k) - 0.01)*fs));
        i1  = min(length(y), i0 + round(0.1*fs));
        env = movmax(abs(y(i0:i1)), round(0.002*fs));
        above = env > 0.3*max(env);
        nb(k) = sum(diff([0; above]) == 1);
    end
end
