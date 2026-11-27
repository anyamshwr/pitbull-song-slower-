% TEST_OLA  Acceptance tests for tsm_ola (Person C, Block 3).
% Run from the project root:   >> test_ola

clear; clc;
fs = 22050;
playSound = false;                % true to listen to tests 4 and 5
if ~exist(fullfile('audio','synth_clicks.wav'), 'file'), make_test_signals(); end
clk  = load_audio(fullfile('audio','synth_clicks.wav'));
tone = load_audio(fullfile('audio','synth_tone.wav'));
if ~exist('out', 'dir'), mkdir('out'); end

%% Test 1: length and shape ---------------------------------------------
for alpha = [0.5 1 1.2 1.8 3]
    y = tsm_ola(clk, alpha);                       % two-argument call must work
    assert(iscolumn(y) && isreal(y), 'Test 1 FAILED: output must be a real column');
    assert(length(y) == round(alpha*length(clk)), ...
           sprintf('Test 1 FAILED: alpha %.1f gives %d samples', alpha, length(y)));
end
fprintf('Test 1  length / shape for alpha in {0.5,1,1.2,1.8,3} : ok\n');

%% Test 2: identity at alpha = 1 ----------------------------------------
y   = tsm_ola(clk, 1);
err = max(abs(y(301:end-300) - clk(301:end-300)));
fprintf('Test 2  identity error            : %.2e  (must be < 1e-9)\n', err);
assert(err < 1e-9, 'Test 2 FAILED: alpha = 1 does not reproduce the input');

%% Test 3: transient timing at alpha = 1.8 ------------------------------
alpha = 1.8;
y     = tsm_ola(clk, alpha);
t0    = 0.25 : 0.5 : 4 - 0.3;                      % onsets in synth_clicks
tExp  = alpha * t0;
tFound = onset_times(y, tExp, fs);
fprintf('Test 3  expected onsets (s): %s\n', mat2str(tExp, 3));
fprintf('        found onsets    (s): %s\n', mat2str(tFound, 3));
assert(all(abs(tFound - tExp) < 0.015), 'Test 3 FAILED: click onsets moved by > 15 ms');

%% Test 4: crispness, N = 256 vs N = 2048 -------------------------------
% A sharp click has its energy packed into a few ms. Measure the fraction
% of each click's energy (in a 100 ms window) that falls inside the
% loudest 10 ms: higher = crisper. Stuttering spreads it out.
y256  = tsm_ola(clk, alpha, 256);
y2048 = tsm_ola(clk, alpha, 2048);
c256  = compactness(y256,  tExp, fs);
c2048 = compactness(y2048, tExp, fs);
cOrig = compactness(clk,   t0,   fs);
fprintf('Test 4  energy in loudest 10 ms:  original %.3f   N=256 %.3f   N=2048 %.3f\n', ...
        cOrig, c256, c2048);
assert(c256 > c2048, 'Test 4 FAILED: N = 256 should be crisper than N = 2048');
audiowrite(fullfile('out','ola_clicks_a1.8_N256.wav'),  0.99*y256 /max(abs(y256)),  fs);
audiowrite(fullfile('out','ola_clicks_a1.8_N2048.wav'), 0.99*y2048/max(abs(y2048)), fs);

%% Test 5: OLA on a tone (evidence for the report) ----------------------
% OLA breaks the phase at each frame join. Measure it: the stretched tone
% should keep the 440 Hz peak but leak energy into other frequencies.
yt  = tsm_ola(tone, alpha);
leakOla  = tonal_leak(yt(round(0.5*fs) + (1:fs)), fs);
leakOrig = tonal_leak(tone(round(0.5*fs) + (1:fs)), fs);
fprintf('Test 5  energy outside the 3 harmonics: original %.4f   OLA a=1.8 %.4f\n', ...
        leakOrig, leakOla);
audiowrite(fullfile('out','ola_tone_a1.8.wav'), 0.99*yt/max(abs(yt)), fs);
fprintf('        saved out/ola_tone_a1.8.wav (should sound rough / warbling)\n');

if playSound
    soundsc(y256, fs);  pause(length(y256)/fs + 0.5);
    soundsc(y2048, fs); pause(length(y2048)/fs + 0.5);
    soundsc(yt, fs);
end
fprintf('\nALL OLA TESTS PASSED\n');


% ---------------------------------------------------------------- helpers
function tf = onset_times(y, tExp, fs)
% First sample above 20% of the local peak, searched +-60 ms around each
% expected onset.
    tf = zeros(size(tExp));
    for k = 1:numel(tExp)
        i0  = max(1, round((tExp(k) - 0.06)*fs));
        i1  = min(length(y), round((tExp(k) + 0.06)*fs));
        seg = abs(y(i0:i1));
        j   = find(seg > 0.2*max(seg), 1);
        tf(k) = (i0 + j - 2) / fs;
    end
end

function c = compactness(y, tc, fs)
% Mean over clicks of (energy in the loudest 10 ms) / (energy in 100 ms).
    c = zeros(size(tc));
    n10 = round(0.010*fs);
    for k = 1:numel(tc)
        i0  = max(1, round((tc(k) - 0.02)*fs));
        i1  = min(length(y), i0 + round(0.1*fs));
        e   = y(i0:i1).^2;
        c(k) = max(movsum(e, n10)) / sum(e);
    end
    c = mean(c);
end

function r = tonal_leak(seg, fs)
% Fraction of spectrum energy NOT within +-5 Hz of 440, 880, 1320 Hz.
    n   = length(seg);
    win = 0.5 - 0.5*cos(2*pi*(0:n-1)'/n);
    S   = abs(fft(seg .* win)).^2;
    S   = S(1:floor(n/2));
    f   = (0:floor(n/2)-1)' * fs / n;
    near = false(size(f));
    for h = [440 880 1320], near = near | abs(f - h) <= 5; end
    r = sum(S(~near)) / sum(S);
end
