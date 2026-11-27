% TEST_HPS  Acceptance tests for hps.m (Person A, Block 1).
% Run from the project root:   >> test_hps
%
% Tests 1-4 use only synthetic signals and must all pass.
% Test 5 (real music) is listening/inspection; list your clips in realClips.

clear; clc;

playSound = false;                 % set true to listen to the synth_mix split
realClips = {};                    % e.g. {'drums_loop.wav', 'pop_clip.wav', 'voice.wav'}

fs = 22050;
if ~exist(fullfile('audio', 'synth_mix.wav'), 'file'), make_test_signals(); end
x     = load_audio(fullfile('audio', 'synth_mix.wav'));
tone  = load_audio(fullfile('audio', 'synth_tone.wav'));
click = load_audio(fullfile('audio', 'synth_clicks.wav'));

%% Test 1: reconstruction xh + xp = x -----------------------------------
[xh, xp] = hps(x);                 % defaults N = 1024, l = 5
err = max(abs(xh + xp - x));
fprintf('Test 1  reconstruction error      : %.2e  (must be < 1e-9)\n', err);
assert(err < 1e-9, 'Test 1 FAILED: xh + xp ~= x');

%% Test 2: lengths and shapes -------------------------------------------
assert(iscolumn(xh) && iscolumn(xp), 'Test 2 FAILED: outputs must be column vectors');
assert(length(xh) == length(x) && length(xp) == length(x), 'Test 2 FAILED: wrong length');
assert(isreal(xh) && isreal(xp), 'Test 2 FAILED: outputs must be real');
xr = x(:).';                       % a row vector input must also give columns
[a1, a2] = hps(xr);
assert(iscolumn(a1) && iscolumn(a2) && length(a1) == length(x), 'Test 2 FAILED: row-vector input');
[s1, s2] = hps(x(1:1000));         % very short signal
assert(length(s1) == 1000 && max(abs(s1 + s2 - x(1:1000))) < 1e-9, 'Test 2 FAILED: short input');
[z1, z2] = hps(zeros(5000, 1));    % silence
assert(all(z1 == 0) && all(z2 == 0), 'Test 2 FAILED: silence must stay silent');
fprintf('Test 2  lengths / shapes / edge cases : ok\n');

%% Test 3: separation quality on synth_mix ------------------------------
% synth_mix is a weighted sum of synth_tone and synth_clicks. Find the two
% weights by least squares, then measure how much of each component ended
% up in the right output.
c    = [tone, click] \ x;
tcmp = c(1) * tone;                % tone as it appears inside the mix
kcmp = c(2) * click;               % clicks as they appear inside the mix
assert(norm(x - tcmp - kcmp) / norm(x) < 1e-3, 'Test 3 setup: mix is not tone + clicks');

toneInXh   = (tcmp' * xh) / (tcmp' * tcmp);   % ideal 1
toneInXp   = (tcmp' * xp) / (tcmp' * tcmp);   % ideal 0
clicksInXp = (kcmp' * xp) / (kcmp' * kcmp);   % ideal 1
clicksInXh = (kcmp' * xh) / (kcmp' * kcmp);   % ideal 0 (small negative is cross-talk)
fprintf('Test 3  tone   -> xh %.4f   tone   -> xp %.4f\n', toneInXh, toneInXp);
fprintf('        clicks -> xp %.4f   clicks -> xh %.4f\n', clicksInXp, clicksInXh);
assert(toneInXh   > 0.98, 'Test 3 FAILED: tone is not in xh (filter directions swapped?)');
assert(clicksInXp > 0.90, 'Test 3 FAILED: clicks are not in xp (filter directions swapped?)');

% single-component sanity checks: tone alone, clicks alone
[th, tp] = hps(tone);   [ch, cp] = hps(click);
fracTone  = (th' * th) / (tone' * tone);
fracClick = (cp' * cp) / (click' * click);
fprintf('        pure tone: %.4f of energy in xh;  pure clicks: %.4f in xp\n', fracTone, fracClick);
assert(fracTone > 0.98 && fracClick > 0.95, 'Test 3 FAILED: single-component check');

% Where does the rest of xp come from? (expected: the start of the signal,
% the 440 -> 659 Hz note change at t = 2 s, and the end of the clip.)
blk = round(0.05 * fs);   nb = floor(length(xp) / blk);
Eb  = sum(reshape(xp(1:nb*blk).^2, blk, nb), 1);
[~, order] = sort(Eb, 'descend');
fprintf('        loudest 50 ms blocks of xp (s): %s\n', mat2str(sort((order(1:8) - 1) * 0.05), 3));

%% Test 4: other parameters still reconstruct exactly -------------------
for N = [512 1024 2048]
    for l = [3 5 9]
        [ha, pa] = hps(x, N, l);
        e = max(abs(ha + pa - x));
        assert(e < 1e-9, sprintf('Test 4 FAILED for N=%d, l=%d (err %.2e)', N, l, e));
    end
end
fprintf('Test 4  N in {512,1024,2048}, l in {3,5,9} reconstruct to < 1e-9 : ok\n');

%% Test 5: real clips (inspection) --------------------------------------
for k = 1:numel(realClips)
    xr = load_audio(fullfile('audio', realClips{k}));
    [rh, rp] = hps(xr);
    e = max(abs(rh + rp - xr));
    assert(e < 1e-9, ['Test 5 FAILED reconstruction on ' realClips{k}]);
    fprintf('Test 5  %-18s recon %.1e   energy in xp: %5.1f %%\n', ...
            realClips{k}, e, 100 * (rp' * rp) / (xr' * xr));
    if playSound
        fprintf('        playing xh then xp ...\n');
        soundsc(rh, fs);  pause(length(rh)/fs + 0.5);
        soundsc(rp, fs);  pause(length(rp)/fs + 0.5);
    end
end
if isempty(realClips)
    fprintf('Test 5  skipped (fill in realClips to run it)\n');
end

if playSound
    soundsc(xh, fs);  pause(length(xh)/fs + 0.5);
    soundsc(xp, fs);
end

fprintf('\nALL HPS TESTS PASSED\n');
