% EXP_RUN_ALL  Experiments + report figures (Person C).
% Put clips in audio/ and list them below (file names without folder).
% Run from the project root:   >> exp_run_all
%
% For every clip it:
%   1. stretches with HP, PV, WS (WSOLA), OLA at alpha = 1.2 and 1.8
%      -> out/results/<clip>_<method>_a<alpha>.wav
%   2. times every run           -> out/runtime.csv  + table in the command window
%   3. saves report figures      -> out/fig/<clip>_*.png
%        _hps.png         spectrogram of x, xh, xp     (paper Fig. 2 style)
%        _waveforms.png   waveform of x, xh, xp        (paper Fig. 1 style)
%        _transient.png   zoom on one drum hit, every method at 1.8 (paper Fig. 3)
%        _spec_HP_PV.png  spectrogram of HP vs PV output at 1.8

clear; clc; close all;

clips  = {'synth_mix.wav', 'pitbull_clip.wav'};   % <-- add your clips here
alphas = [1.2 1.8];
fs     = 22050;

for d = {'out', fullfile('out','results'), fullfile('out','fig')}
    if ~exist(d{1}, 'dir'), mkdir(d{1}); end
end

methods = {'HP',  @(x,a) tsm_hp(x, a); ...
           'PV',  @(x,a) tsm_pv(x, a); ...
           'WS',  @(x,a) tsm_wsola(x, a); ...
           'OLA', @(x,a) tsm_ola(x, a, 256)};

rt = {};                                         % runtime rows
for c = 1:numel(clips)
    if ~exist(fullfile('audio', clips{c}), 'file')
        fprintf('skipping %s (not in audio/)\n', clips{c});  continue;
    end
    [~, name] = fileparts(clips{c});
    x   = load_audio(fullfile('audio', clips{c}));
    dur = length(x) / fs;
    fprintf('\n=== %s  (%.1f s) ===\n', name, dur);

    % ---- 1 + 2: stretch every way, time it, save audio -----------------
    Y = struct();
    for a = alphas
        for m = 1:size(methods, 1)
            tic;  y = methods{m, 2}(x, a);  sec = toc;
            y = 0.99 * y / max(abs(y));          % same loudness for listening
            audiowrite(fullfile('out','results', ...
                sprintf('%s_%s_a%.1f.wav', name, methods{m,1}, a)), y, fs);
            rt(end+1, :) = {name, methods{m,1}, a, sec, dur};   %#ok<SAGROW>
            fprintf('  %-4s a=%.1f  %6.2f s   (%5.1f s per minute of audio)\n', ...
                    methods{m,1}, a, sec, 60*sec/dur);
            if a == 1.8, Y.(methods{m,1}) = y; end
        end
    end

    % ---- 3a: HPS spectrograms ------------------------------------------
    [xh, xp] = hps(x);
    f = figure('Visible', 'off', 'Position', [100 100 900 750]);
    sigs = {x, xh, xp};
    ttl  = {'Original x', 'Harmonic part x_h (sustained notes)', ...
            'Percussive part x_p (drum hits)'};
    ref = [];                                    % same colour scale for all 3
    for k = 1:3
        subplot(3, 1, k);  ref = plot_spec(sigs{k}, fs, ref);  title(ttl{k});
    end
    save_fig(f, name, 'hps');

    % ---- 3b: waveforms ---------------------------------------------------
    f = figure('Visible', 'off', 'Position', [100 100 900 600]);
    tt = (0:length(x)-1)' / fs;
    for k = 1:3
        subplot(3, 1, k);  plot(tt, sigs{k});  ylim([-1 1]);  xlim([0 dur]);
        title(ttl{k});  ylabel('amplitude');
    end
    xlabel('time (s)');
    save_fig(f, name, 'waveforms');

    % ---- 3c: transient zoom, every method at 1.8 -------------------------
    % Plotted high-passed (2nd difference) so steady notes drop out and
    % only the hits remain. All panels share one amplitude scale.
    t0 = strongest_hit(x, fs);                   % time of the clearest drum hit
    yl = max(abs(diff(x, 2)));
    f  = figure('Visible', 'off', 'Position', [100 100 900 900]);
    subplot(5, 1, 1);  plot_zoom(x, t0, fs, yl);
    title(sprintf('Original, hit at %.2f s', t0));
    for m = 1:size(methods, 1)
        subplot(5, 1, m+1);  plot_zoom(Y.(methods{m,1}), 1.8*t0, fs, yl);
        title(sprintf('%s, \\alpha = 1.8 (hit expected at %.2f s)', methods{m,1}, 1.8*t0));
    end
    xlabel('time relative to the hit (ms)');
    save_fig(f, name, 'transient');

    % ---- 3d: spectrogram HP vs PV at 1.8 ---------------------------------
    f = figure('Visible', 'off', 'Position', [100 100 900 550]);
    subplot(2, 1, 1);  ref = plot_spec(Y.HP, fs, []);  title('HP, \alpha = 1.8 (hits stay as thin vertical lines)');
    subplot(2, 1, 2);  plot_spec(Y.PV, fs, ref);  title('PV, \alpha = 1.8 (hits smeared in time)');
    save_fig(f, name, 'spec_HP_PV');
end

% ---- runtime table --------------------------------------------------------
fid = fopen(fullfile('out','runtime.csv'), 'w');
fprintf(fid, 'clip,method,alpha,seconds,clip_seconds,seconds_per_minute\n');
for r = 1:size(rt, 1)
    fprintf(fid, '%s,%s,%.1f,%.3f,%.2f,%.2f\n', rt{r,:}, 60*rt{r,4}/rt{r,5});
end
fclose(fid);
fprintf('\nSaved: out/results/*.wav, out/fig/*.png, out/runtime.csv\n');
fprintf('Paper reference (HP, MATLAB): ~4 s per minute at 1.2, ~5 s at 1.8\n');


% ===================================================================== helpers
function ref = plot_spec(y, fs, ref)
% Magnitude spectrogram in dB with our own STFT (N = 1024, hop 256).
% ref = top of the colour scale in dB; pass [] to use this signal's max.
    N = 1024;  H = 256;
    w = 0.5 - 0.5*cos(2*pi*(0:N-1)'/N);
    y = [y(:); zeros(N,1)];
    T = floor((length(y) - N)/H) + 1;
    S = zeros(N/2+1, T);
    for t = 0:T-1
        F = fft(y(t*H + (1:N)) .* w);
        S(:, t+1) = abs(F(1:N/2+1));
    end
    S = 20*log10(S + 1e-6);
    imagesc((0:T-1)*H/fs, (0:N/2)*fs/N, S);
    axis xy;  ylim([0 5000]);
    if isempty(ref), ref = max(S(:)); end
    caxis([ref-70 ref]);  colorbar;
    ylabel('frequency (Hz)');  xlabel('time (s)');
end

function plot_zoom(y, tc, fs, yl)
% 200 ms of the high-passed waveform centred on time tc, x axis in ms.
    d  = [0; diff(y, 2); 0];                     % same length as y
    i0 = max(1, round((tc - 0.1)*fs));
    i1 = min(length(y), round((tc + 0.1)*fs));
    plot(((i0:i1) - 1)/fs*1000 - tc*1000, d(i0:i1));
    xlim([-100 100]);  ylim([-yl yl]);  ylabel('high-passed');
end

function t0 = strongest_hit(x, fs)
% Time of the sharpest onset: biggest jump in high-passed energy,
% ignoring the first and last 0.3 s.
    d   = diff(x, 2);
    env = movsum(d.^2, round(0.005*fs));
    jump = [0; max(diff(env), 0)];
    edge = round(0.3*fs);
    jump([1:edge, end-edge:end]) = 0;
    [~, i] = max(jump);
    t0 = i / fs;
end

function save_fig(f, name, tag)
    print(f, fullfile('out','fig', sprintf('%s_%s.png', name, tag)), '-dpng', '-r150');
    close(f);
end
