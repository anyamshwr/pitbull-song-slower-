function y = tsm_hp(x, alpha)
% Harmonic-percussive TSM (Driedger, Mueller, Ewert 2014)
    x = x(:);
    Lout = round(alpha * length(x));

    [xh, xp] = hps(x, 1024, 5);          % Block 1 (Person A)
    yh = tsm_pv(xh, alpha, 4096, true);  % Block 2 (Person B)
    yp = tsm_ola(xp, alpha, 256);        % Block 3 (Person C)

    yh = fixlen(yh, Lout);
    yp = fixlen(yp, Lout);
    y  = yh + yp;

    m = max(abs(y));
    if m > 1, y = 0.99 * y / m; end      % avoid clipping when saving
end

function z = fixlen(z, L)
    z = z(:);
    if length(z) >= L, z = z(1:L); else, z = [z; zeros(L - length(z), 1)]; end
end
