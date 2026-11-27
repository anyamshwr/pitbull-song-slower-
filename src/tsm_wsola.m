function y = tsm_wsola(x, alpha)
% TSM_WSOLA  Waveform-similarity overlap-add (Verhelst & Roelands 1993).
%   Baseline for the comparison (Person C).
%   y = tsm_wsola(x, alpha)
%
%   Like OLA, but each analysis frame may move up to +-Delta samples from its
%   nominal position. We pick the position whose waveform is most similar to
%   the natural continuation of the previously chosen frame, so periodic
%   signals line up across frame joins.
%
%   Parameters (our choice, paper gives none): N = 1024, Hs = N/2, Delta = N/4.
%   Output: column vector with exactly round(alpha*length(x)) samples.

    N     = 1024;
    Hs    = N/2;                                  % synthesis hop
    Delta = N/4;                                  % search tolerance

    x    = x(:);
    L    = length(x);
    Lout = round(alpha * L);
    w    = 0.5 - 0.5*cos(2*pi*(0:N-1)'/N);        % periodic Hann

    % ---- shared synthesis recipe, with extra back padding for the search ----
    xpad = [zeros(N/2,1); x; zeros(3*N,1)];
    T    = ceil((Lout + N)/Hs) + 1;
    ypad = zeros((T-1)*Hs + N, 1);
    wsum = zeros(size(ypad));
    Lx   = length(xpad);

    pPrev = 0;
    for t = 0:T-1
        a = round(t * Hs / alpha);                % nominal 0-based start

        if t == 0
            p = a;                                % first frame: no search
        else
            % natural continuation of the previous chosen frame in the INPUT
            ref = grab_frame(xpad, pPrev + Hs, N);

            % candidate starts p in [a-Delta, a+Delta], clamped to p >= 0
            pmin = max(a - Delta, 0);
            pmax = a + Delta;
            seg  = grab_frame(xpad, pmin, pmax - pmin + N);   % covers all candidates

            % dot product of every candidate with ref, in one go
            c = conv(seg, flipud(ref), 'valid');            % length pmax-pmin+1
            % energy of every candidate (normalised cross-correlation)
            e = movsum(seg.^2, N, 'Endpoints', 'discard');
            score = c ./ sqrt(max(e, 1e-12));               % norm(ref) is constant

            if all(e < 1e-12)
                p = a;                            % silence: keep nominal position
            else
                % among equal scores prefer the one closest to the nominal
                % position (makes alpha = 1 an exact identity)
                cand = (pmin:pmax)';
                best = max(score);
                tie  = find(score >= best - 1e-12 * max(abs(best), 1));
                [~, j] = min(abs(cand(tie) - a));
                p = cand(tie(j));
            end
        end

        fr  = grab_frame(xpad, p, N);
        idx = t*Hs + (1:N);
        ypad(idx) = ypad(idx) + fr .* w;          % windowed once
        wsum(idx) = wsum(idx) + w;
        pPrev = p;
    end

    ypad = ypad ./ max(wsum, 1e-3);
    y    = ypad(N/2+1 : min(N/2+Lout, length(ypad)));
    if length(y) < Lout
        y = [y; zeros(Lout - length(y), 1)];
    end
end


function fr = grab_frame(xpad, s, n)
% GRAB_FRAME  n samples of xpad starting at 0-based index s, zero-padded at the end.
    fr = zeros(n, 1);
    if s < length(xpad)
        seg = xpad(s+1 : min(s+n, length(xpad)));
        fr(1:length(seg)) = seg;
    end
end
