function make_test_signals()
    fs = 22050; dur = 4;
    t = (0:round(dur*fs)-1)' / fs;
    if ~exist('audio', 'dir'), mkdir('audio'); end

    % sustained harmonic notes: A4 then E5, with 2 extra harmonics
    f0 = 440 * ones(size(t));  f0(t >= 2) = 659.26;
    ph = 2*pi*cumsum(f0) / fs;
    tone = 0.5*sin(ph) + 0.25*sin(2*ph) + 0.12*sin(3*ph);

    % drum-like noise bursts every 0.5 s starting at 0.25 s
    rng(0);
    clicks = zeros(size(t));
    n = 300;
    for s = 0.25 : 0.5 : dur - 0.3
        i = round(s*fs) + 1;
        clicks(i:i+n-1) = clicks(i:i+n-1) + randn(n, 1) .* exp(-(0:n-1)'/60);
    end

    mix = 0.5*tone/max(abs(tone)) + 0.5*clicks/max(abs(clicks));
    wr = @(name, s) audiowrite(fullfile('audio', name), 0.9*s/max(abs(s)), fs);
    wr('synth_tone.wav', tone);
    wr('synth_clicks.wav', clicks);
    wr('synth_mix.wav', mix);
end
