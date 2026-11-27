% TEST_INTEGRATION  Day-4 acceptance tests (shared context, section 10).
% Run from the project root after all blocks are in the folder.

clear; clc;
fs = 22050;
if ~exist(fullfile('audio','synth_mix.wav'), 'file'), make_test_signals(); end
x    = load_audio(fullfile('audio','synth_mix.wav'));
tone = load_audio(fullfile('audio','synth_tone.wav'));
if ~exist('out', 'dir'), mkdir('out'); end

%% 1. Contract -----------------------------------------------------------
[xh, xp] = hps(x);
for alpha = [0.5 1 1.2 1.8 3]
    Lout = round(alpha*length(x));
    a = tsm_pv(xh, alpha);  b = tsm_ola(xp, alpha);  c = tsm_wsola(x, alpha);
    assert(iscolumn(a) && length(a) == Lout, 'Test 1 FAILED: tsm_pv contract');
    assert(iscolumn(b) && length(b) == Lout, 'Test 1 FAILED: tsm_ola contract');
    assert(iscolumn(c) && length(c) == Lout, 'Test 1 FAILED: tsm_wsola contract');
end
fprintf('Test 1  contract (pv, ola, wsola; 5 alphas)  : ok\n');

%% 2. HPS reconstruction -------------------------------------------------
e = max(abs(xh + xp - x));
fprintf('Test 2  hps reconstruction error  : %.2e\n', e);
assert(e < 1e-9, 'Test 2 FAILED');

%% 3. Identity of the full pipeline -------------------------------------
y = tsm_hp(x, 1.0);
e = max(abs(y(4097:end-4096) - x(4097:end-4096)));
fprintf('Test 3  tsm_hp identity error     : %.2e  (must be < 1e-6)\n', e);
assert(e < 1e-6, 'Test 3 FAILED: bug is in the split or a normalisation');

%% 4. Pitch --------------------------------------------------------------
yt  = tsm_pv(tone, 1.8);
seg = yt(round(0.5*fs) + (1:fs));
win = 0.5 - 0.5*cos(2*pi*(0:fs-1)'/fs);
S   = abs(fft(seg .* win));
[~, i] = max(S(1:floor(fs/2)));
f = (i-1) * fs / fs;
fprintf('Test 4  pitch after 1.8x stretch   : %.1f Hz, length ratio %.3f\n', ...
        f, length(yt)/length(tone));
assert(abs(f - 440) <= 2, 'Test 4 FAILED');

%% 5. Transient timing through the full pipeline ------------------------
alpha = 1.8;
yhp   = tsm_hp(x, alpha);
tExp  = alpha * (0.25 : 0.5 : 4 - 0.3);
tf    = zeros(size(tExp));
for k = 1:numel(tExp)                     % onset = first rise in the HF part
    i0 = round((tExp(k) - 0.06)*fs);  i1 = round((tExp(k) + 0.06)*fs);
    d  = abs(diff(yhp(i0:i1), 2));        % 2nd difference: tone ~0, noise clicks large
    j  = find(d > 0.3*max(d), 1);
    tf(k) = (i0 + j - 1) / fs;
end
fprintf('Test 5  expected onsets (s): %s\n', mat2str(tExp, 3));
fprintf('        found onsets    (s): %s\n', mat2str(tf, 3));
assert(all(abs(tf - tExp) < 0.015), 'Test 5 FAILED');

%% 6. Listening files + objective transient sharpness --------------------
methods = {'HP',  @(s) tsm_hp(s, alpha); ...
           'PV',  @(s) tsm_pv(s, alpha); ...
           'WS',  @(s) tsm_wsola(s, alpha); ...
           'OLA', @(s) tsm_ola(s, alpha, 256)};
fprintf('Test 6  click sharpness: share of each click''s high-passed energy\n');
fprintf('        (100 ms window) inside its loudest 10 ms. Original = 0.84.\n');
for m = 1:size(methods, 1)
    ym = methods{m, 2}(x);
    d  = diff(ym, 2);                     % high-pass: removes the steady tone
    cf = zeros(size(tExp));
    for k = 1:numel(tExp)                 % window starts 5 ms before the onset
        i0 = round(tExp(k)*fs) - round(0.005*fs);
        e  = d(i0 : i0 + round(0.1*fs)).^2;
        cf(k) = max(movsum(e, round(0.01*fs))) / sum(e);
    end
    fprintf('        %-4s %.3f\n', methods{m, 1}, mean(cf));
    audiowrite(fullfile('out', sprintf('cmp_synth_mix_%s_a1.8.wav', methods{m,1})), ...
               0.99*ym/max(abs(ym)), fs);
end
fprintf('        saved out/cmp_synth_mix_<method>_a1.8.wav  -> listen: HP crisp + steady,\n');
fprintf('        PV smeared clicks, OLA rough tone, WS doubled clicks\n');

fprintf('\nALL INTEGRATION TESTS PASSED\n');
