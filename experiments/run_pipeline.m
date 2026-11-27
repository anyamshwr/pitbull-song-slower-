if ~exist('out', 'dir'), mkdir('out'); end
name = 'synth_mix';
x = load_audio(fullfile('audio', [name '.wav']));
for alpha = [1.0 1.2 1.8]
    y = tsm_hp(x, alpha);
    assert(length(y) == round(alpha*length(x)), 'length mismatch');
    audiowrite(fullfile('out', sprintf('%s_hp_a%.1f.wav', name, alpha)), y, 22050);
    fprintf('alpha=%.1f  ok  (%d samples)\n', alpha, length(y));
end
