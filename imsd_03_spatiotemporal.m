CH1 = A;

F = fftn(CH1);
G0 = F .* conj(F);
G0 = ifftn(G0);
G0 = fftshift(fftshift(G0,1),2);
STCorr = real(G0)/(mean(CH1(:))^2 * numel(CH1)) - 1;

STCorr_a=STCorr(P/2+1+int64(a1/px_size),P/2+1+int64(a1/px_size),:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Preparazione ---
% [ny,nx,nt] = size(STCorr_a);
ny=size(STCorr_a,1);
nx=size(STCorr_a,2);
nt=TauLimit;

xv = a1;
yv = a1; % se i due assi sono uguali; altrimenti usa b1

[X,Y] = meshgrid(xv,yv);
xy = cat(3,X,Y);

centerIdx = floor(ny/2)+1;
centerJdx = floor(nx/2)+1;

frame0 = STCorr_a(:,:,1);
amp0 = max(frame0(:)) - min(frame0(:));
offset0 = median(frame0(:));
sigma0 = max(range(xv),range(yv))/10;

start0 = [amp0, 0, 0, sigma0, offset0]; % [A, x0, y0, sigma, offset]

rmax = max([abs(xv(:)); abs(yv(:))]);
lb = [0, -rmax, -rmax, 0, min(frame0(:))-abs(offset0)*10];
ub = [inf,  rmax,  rmax, max(range(xv),range(yv)), max(frame0(:))+abs(offset0)*10];

options = optimoptions('lsqcurvefit', ...
    'Algorithm','trust-region-reflective', ...
    'Display','off', ...
    'MaxIterations',2000, ...
    'MaxFunctionEvaluations',10000, ...
    'ScaleProblem','jacobian');

modelFunG = @(pg,xy) pg(1) .* exp(-((xy(:,:,1)-pg(2)).^2 + (xy(:,:,2)-pg(3)).^2) ./ (pg(4).^2)) + pg(5);

Gaussb = zeros(ny,nx,nt);
x_G = zeros(nt,1);
y_G = zeros(nt,1);
s_G = zeros(nt,1);
o_G = zeros(nt,1);

startG = start0;

for i = 1:nt
    dataFrame = STCorr_a(:,:,i);

    % aggiorna start dai parametri precedenti (se disponibili)
    if i > 1 && all(isfinite([x_G(i-1), y_G(i-1), s_G(i-1)]))
        startG(2) = x_G(i-1);
        startG(3) = y_G(i-1);
        startG(4) = sqrt(s_G(i-1));
        startG(5) = o_G(i-1);
    end

    % TypicalX non deve avere zeri
    tx = startG;
    tx(abs(tx) < eps) = 1;   
    % options = optimoptions(options,'TypicalX',tx);

    try
        [G,resnorm, residual, exitflag, output, lambda, J] = lsqcurvefit(modelFunG, startG, xy, dataFrame, lb, ub, options);
        if exitflag <= 0
            warning('Frame %d: fit non convergente. Fallback a start0.', i);
            [G,~,~,~] = lsqcurvefit(modelFunG, start0, xy, dataFrame, lb, ub, options);
        end
    catch ME
        warning('Frame %d: errore nel fit (%s). Uso parametri precedenti.', i, ME.message);
        if i > 1
            G = [Gaussb(centerIdx,centerJdx,i-1), x_G(i-1), y_G(i-1), sqrt(s_G(i-1)), o_G(i-1)];
        else
            G = start0;
        end
    end

    ci = nlparci(G, residual(:), 'jacobian', J);

% Errore stimato come metà della larghezza dell'intervallo (circa ±1σ)
param_errors = (ci(:,2) - ci(:,1)) / 2;

    Gaussb(:,:,i) = modelFunG(G, xy);
    amp_G(i,1)  = G(1);
    x_G(i,1) = G(2);
    y_G(i,1) = G(3);
    s_G(i,1) = (G(4)).^2;
    o_G(i,1) = G(5);

    amp_G_err(i,1) = param_errors(1);
    x_G_err(i,1)   = param_errors(2);
    y_G_err(i,1)   = param_errors(3);
    s_G_err(i,1)   = param_errors(4).*G(4).*2;
    o_G_err(i,1)   = param_errors(5);

    % estrazione sezioni
    idx_x = centerJdx + round(x_G(i) / px_size);
    idx_y = centerIdx + round(y_G(i) / px_size);

    idx_x = max(1, min(nx, idx_x));
    idx_y = max(1, min(ny, idx_y));

    Gauss_y(:,i) = Gaussb(:, idx_x, i);
    Gauss_x(:,i) = Gaussb(idx_y, :, i)';
end

R2_G = zeros(nt,1);

for i = 1:TauLimit
    Gfit = Gaussb(:,:,i);
    Gdata = STCorr_a(:,:,i);
    SSres = sum((Gdata(:)-Gfit(:)).^2);
    SStot = sum((Gdata(:)-mean(Gdata(:))).^2);
    R2_G(i) = 1 - SSres/SStot;
end

spatiotemporalCorrelationViewer(time_vec, frame_time, px_size, a1, b1, TauLimit, STCorr_a(:,:,1:TauLimit), Gaussb,...
    amp_G, s_G, o_G, R2_G, amp_G_err, s_G_err, o_G_err)