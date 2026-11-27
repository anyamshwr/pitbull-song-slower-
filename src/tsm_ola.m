function y = tsm_ola(x, alpha, N)
% TSM_OLA  Plain overlap-add time-scale modification (Person C, Block 3).
%   y = tsm_ola(x, alpha)        N = 256
%   y = tsm_ola(x, alpha, N)
%
%   Paper Eq. (5) with y_t = x_t: frames are copied unchanged from the input,
%   re-spaced to the synthesis hop, windowed once, overlap-added and
%   normalised by the summed window. Used on the percussive part x_p only.
%
%   Output: column vector with exactly round(alpha*length(x)) samples.

    if nargin < 3 || isempty(N), N = 256; end

    x    = x(:);
    L    = length(x);
    Hs   = N/2;                                   % synthesis hop (fixed)
    Lout = round(alpha * L);
    w    = 0.5 - 0.5*cos(2*pi*(0:N-1)'/N);        % periodic Hann

    % ---- shared synthesis recipe (section 6 of shared context) ----
    xpad = [zeros(N/2,1); x; zeros(2*N,1)];
    T    = ceil((Lout + N)/Hs) + 1;               % number of frames
    ypad = zeros((T-1)*Hs + N, 1);
    wsum = zeros(size(ypad));

    for t = 0:T-1
        % 1. analysis frame at a_t = round(t*Hs/alpha), zero-padded to N
        a   = round(t * Hs / alpha);              % 0-based start in xpad
        fr  = zeros(N,1);
        if a < length(xpad)                       % else: frame is all zeros
            seg = xpad(a+1 : min(a+N, length(xpad)));
            fr(1:length(seg)) = seg;
        end

        % 2. window once, place at t*Hs, accumulate the window (not w.^2)
        idx = t*Hs + (1:N);
        ypad(idx) = ypad(idx) + fr .* w;
        wsum(idx) = wsum(idx) + w;
    end

    % 3. normalise, trim, fix length
    ypad = ypad ./ max(wsum, 1e-3);
    y    = ypad(N/2+1 : min(N/2+Lout, length(ypad)));
    if length(y) < Lout
        y = [y; zeros(Lout - length(y), 1)];
    end
end
