function [xh, xp] = hps(x, N, l)
% HPS  Harmonic-percussive separation by median filtering (Fitzgerald 2010,
%      as used in Driedger, Mueller, Ewert, IEEE SPL 2014, Section II).
%
%   [xh, xp] = hps(x, N, l)
%
%   x   mono signal (any orientation; treated as a column)
%   N   frame size            (default 1024)
%   l   median half-length    (default 5, so the filter length is 2l+1 = 11)
%
%   xh  harmonic part,   column vector, same length as x
%   xp  percussive part, column vector, same length as x
%   Guarantee: max(abs(xh + xp - x)) < 1e-9
%
%   Hop is H = N/4, window is a periodic Hann (our choices, see shared context).

    if nargin < 2 || isempty(N), N = 1024; end
    if nargin < 3 || isempty(l), l = 5;    end

    x = x(:);                         % force column vector
    L = length(x);
    H = N / 4;                        % hop size (75% overlap)
    K = N/2 + 1;                      % one-sided spectrum: bins 0 ... N/2
    w = 0.5 - 0.5*cos(2*pi*(0:N-1)'/N);   % periodic Hann window

    % ---- Step 1: padding -------------------------------------------------
    % N zeros in front and N+H behind, so every original sample is covered
    % by four frames and the edges reconstruct exactly.
    xpad = [zeros(N,1); x; zeros(N+H,1)];
    T    = floor((length(xpad) - N)/H) + 1;      % number of frames

    % ---- Step 2: STFT (own implementation) -------------------------------
    % Rows = frequency, columns = time.
    X = zeros(K, T);
    for t = 0:T-1
        fr = xpad(t*H+1 : t*H+N) .* w;           % windowed frame
        F  = fft(fr);
        X(:, t+1) = F(1:K);                      % keep bins 0 ... N/2
    end

    % ---- Step 3: magnitude and median filters ----------------------------
    Y  = abs(X);
    Yh = movmedian(Y, 2*l+1, 2);     % along columns = along TIME (horizontal)
    Yp = movmedian(Y, 2*l+1, 1);     % along rows    = along FREQ (vertical)

    % ---- Step 4: binary masks (paper's comparison operators) -------------
    Mh = (Yh >= Yp);
    Mp = (Yp >  Yh);                 % Mh + Mp = 1 everywhere

    % ---- Step 5: mask the COMPLEX spectrogram ----------------------------
    Xh = X .* Mh;
    Xp = X .* Mp;

    % ---- Step 6: inverse STFT (least-squares overlap-add) ----------------
    yh = istft_lsee(Xh, N, H);
    yp = istft_lsee(Xp, N, H);

    % ---- Step 7: remove padding, return columns of length L --------------
    xh = yh(N+1 : N+L);
    xp = yp(N+1 : N+L);
    xh = xh(:);
    xp = xp(:);
end


function ypad = istft_lsee(X, N, H)
% ISTFT_LSEE  Least-squares inverse STFT (Griffin & Lim 1984).
%   X is a K x T one-sided spectrogram (K = N/2+1). The result is the
%   padded time signal of length (T-1)*H + N; the caller removes the padding.
    T = size(X, 2);
    w = 0.5 - 0.5*cos(2*pi*(0:N-1)'/N);          % same window as analysis

    ypad = zeros((T-1)*H + N, 1);
    wsum = zeros((T-1)*H + N, 1);

    for t = 0:T-1
        Xcol = X(:, t+1);
        % rebuild the full conjugate-symmetric spectrum from bins 0 ... N/2
        full = [Xcol; conj(Xcol(end-1:-1:2))];
        fr   = real(ifft(full));                  % frame in time domain
        idx  = t*H+1 : t*H+N;
        ypad(idx) = ypad(idx) + fr .* w;          % synthesis window + add
        wsum(idx) = wsum(idx) + w.^2;             % analysis w and synthesis w
    end

    ypad = ypad ./ max(wsum, 1e-8);               % least-squares normalisation
end
