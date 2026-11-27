function x = load_audio(path)
    [x, fs] = audioread(path);
    x = mean(x, 2);                      % mono
    if fs ~= 22050
        x = resample(x, 22050, fs);      % Signal Processing Toolbox
    end
    x = x(:);
    x = 0.9 * x / max(abs(x));           % peak normalise
end
