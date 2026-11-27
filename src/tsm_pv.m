function y = tsm_pv(x, alpha, N, lock)
% TSM_PV  Phase vocoder time-scale modification (Person B's block).
%   y = tsm_pv(x, alpha)            N = 4096, lock = true
%   y = tsm_pv(x, alpha, N, lock)
%   Output: column vector with exactly round(alpha*length(x)) samples.
%
%   STATUS: basic vocoder (lock = false) is complete.
%           Identity phase locking is a TODO (see marked block below).

    if nargin < 3, N = 4096; end
    if nargin < 4, lock = true; end

    x     = x(:);
    L     = length(x);
    Hs    = N/4;                                  % synthesis hop (fixed)
    K     = N/2 + 1;                              % one-sided bins k = 0..K-1
    Lout  = round(alpha * L);
    w     = 0.5 - 0.5*cos(2*pi*(0:N-1)'/N);       % periodic Hann
    omega = 2*pi*(0:K-1)' / N;                    % bin centre freq, rad/sample

    % ---- shared synthesis recipe (section 6 of shared context) ----
    xpad = [zeros(N/2,1); x; zeros(2*N,1)];
    T    = ceil((Lout + N)/Hs) + 1;               % number of frames
    a    = round((0:T-1) * Hs / alpha);           % 0-based analysis starts
    ypad = zeros((T-1)*Hs + N, 1);
    wsum = zeros(size(ypad));

    phiPrev    = [];                              % analysis phase of previous frame
    phiOutPrev = [];                              % output phase of previous frame

    for t = 1:T                                   % MATLAB t  <->  frame number t-1
        % 1. cut and window the frame (zero-pad if it runs past the end)
        s   = a(t);
        fr  = zeros(N,1);
        seg = xpad(s+1 : min(s+N, length(xpad)));
        fr(1:length(seg)) = seg;

        X   = fft(fr .* w);
        X   = X(1:K);
        mag = abs(X);
        phi = angle(X);

        % 2. phase update
        if t == 1
            phiOut = phi;
        else
            Ha   = a(t) - a(t-1);                          % ACTUAL analysis hop
            dphi = phi - phiPrev - omega*Ha;               % deviation from expected advance
            dphi = mod(dphi + pi, 2*pi) - pi;              % wrap to [-pi, pi)
            inst = omega + dphi / Ha;                      % instantaneous frequency

            if ~lock
                phiOut = phiOutPrev + Hs*inst;             % basic vocoder
            else
                 
                     isPeak = (mag >= movmax(mag, 5)) & (mag > 0);   % bigger than 2 neighbours each side
pk = find(isPeak);

if isempty(pk)
    phiOut = phiOutPrev + Hs*inst;               % no peaks (silence): basic update
else
    if numel(pk) >= 2
        nearest = interp1(pk, pk, (1:K)', 'nearest', 'extrap');  % nearest peak for every bin
    else
        nearest = pk(1) * ones(K,1);
    end
    phiPeak = zeros(K,1);
    phiPeak(pk) = phiOutPrev(pk) + Hs*inst(pk);  % advance phase at the peaks only
    phiOut  = phiPeak(nearest) + (phi - phi(nearest));  % lock the other bins to their peak
 end
                     
            end
        end
        phiPrev    = phi;
        phiOutPrev = phiOut;

        % 3. resynthesis: new spectrum -> time frame -> window -> overlap-add
        Y    = mag .* exp(1i*phiOut);
        full = [Y; conj(Y(end-1:-1:2))];              % conjugate-symmetric spectrum
        yf   = real(ifft(full)) .* w;
        idx  = (t-1)*Hs + (1:N);
        ypad(idx) = ypad(idx) + yf;
        wsum(idx) = wsum(idx) + w.^2;                 % windowed twice -> w.^2
    end

    % normalise, trim, fix length
    ypad = ypad ./ max(wsum, 1e-3);
    y    = ypad(N/2+1 : min(N/2+Lout, length(ypad)));
    if length(y) < Lout
        y = [y; zeros(Lout - length(y), 1)];
    end
end
