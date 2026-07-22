CH1 = A;
F = fftn(CH1);
G0 = F .* conj(F);
G0 = ifftn(G0);
G0 = fftshift(fftshift(G0,1),2);
STCorr = real(G0)/(mean(CH1(:))^2 * numel(CH1)) - 1;
STCorr_a = STCorr(P/2+1+int64(a1/px_size),P/2+1+int64(a1/px_size),:);

ny=size(STCorr_a,1); nx=size(STCorr_a,2); nt=TauLimit;
xv = a1; yv = a1;
[X,Y] = meshgrid(xv,yv);
xy = cat(3,X,Y);
centerIdx = floor(ny/2)+1;
centerJdx = floor(nx/2)+1;

frame0 = STCorr_a(:,:,1);
amp0 = max(frame0(:)) - min(frame0(:));
offset0 = median(frame0(:));
sigma0 = max(range(xv),range(yv))/10;
start0 = [amp0, 0, 0, sigma0, offset0];
rmax = max([abs(xv(:)); abs(yv(:))]);
lb = [0, -rmax, -rmax, 0, min(frame0(:))-abs(offset0)*10];
ub = [inf,  rmax,  rmax, max(range(xv),range(yv)), max(frame0(:))+abs(offset0)*10];
options = optimoptions('lsqcurvefit', 'Algorithm','trust-region-reflective', ...
    'Display','off', 'MaxIterations',2000, 'MaxFunctionEvaluations',10000, 'ScaleProblem','jacobian');
modelFunG = @(pg,xy) pg(1) .* exp(-((xy(:,:,1)-pg(2)).^2 + (xy(:,:,2)-pg(3)).^2) ./ (pg(4).^2)) + pg(5);

if ~exist('fitt2comp','var') || ~fitt2comp
%%% ==================== STANDARD single-Gaussian (original) ====================
Gaussb = zeros(ny,nx,nt);
x_G = zeros(nt,1); y_G = zeros(nt,1); s_G = zeros(nt,1); o_G = zeros(nt,1);
startG = start0;
for i = 1:nt
    dataFrame = STCorr_a(:,:,i);
    if i > 1 && all(isfinite([x_G(i-1), y_G(i-1), s_G(i-1)]))
        startG(2)=x_G(i-1); startG(3)=y_G(i-1); startG(4)=sqrt(s_G(i-1)); startG(5)=o_G(i-1);
    end
    try
        [G,~, residual, exitflag, ~, ~, J] = lsqcurvefit(modelFunG, startG, xy, dataFrame, lb, ub, options);
        if exitflag <= 0
            [G,~, residual, ~, ~, ~, J] = lsqcurvefit(modelFunG, start0, xy, dataFrame, lb, ub, options);
        end
    catch
        if i > 1, G = [Gaussb(centerIdx,centerJdx,i-1), x_G(i-1), y_G(i-1), sqrt(s_G(i-1)), o_G(i-1)];
        else,      G = start0; end
        residual = zeros(numel(dataFrame),1); J = zeros(numel(dataFrame),5);
    end
    try, ci = nlparci(G, residual(:), 'jacobian', J); param_errors = (ci(:,2)-ci(:,1))/2;
    catch, param_errors = nan(5,1); end
    Gaussb(:,:,i) = modelFunG(G, xy);
    amp_G(i,1)=G(1); x_G(i,1)=G(2); y_G(i,1)=G(3); s_G(i,1)=G(4)^2; o_G(i,1)=G(5);
    amp_G_err(i,1)=param_errors(1); x_G_err(i,1)=param_errors(2); y_G_err(i,1)=param_errors(3);
    s_G_err(i,1)=param_errors(4)*G(4)*2; o_G_err(i,1)=param_errors(5);
end
R2_G = zeros(nt,1);
for i = 1:nt
    Gfit=Gaussb(:,:,i); Gdata=STCorr_a(:,:,i);
    R2_G(i)=1-sum((Gdata(:)-Gfit(:)).^2)/sum((Gdata(:)-mean(Gdata(:))).^2);
end
spatiotemporalCorrelationViewer(time_vec, frame_time, px_size, a1, b1, TauLimit, ...
    STCorr_a(:,:,1:TauLimit), Gaussb, amp_G, s_G, o_G, R2_G, amp_G_err, s_G_err, o_G_err)

else
%%% ==================== TWO-COMPONENT fit (fitt2comp=true) ====================
%%% remove the zero-lag central shot-noise spike (as in imsd_02)
for i = 1:nt
    fr = STCorr_a(:,:,i);
    nb = fr(centerIdx-1:centerIdx+1, centerJdx-1:centerJdx+1); nb(2,2) = -inf;
    fr(centerIdx,centerJdx) = max(nb(:));
    STCorr_a(:,:,i) = fr;
end
r2v = X(:).^2 + Y(:).^2;

%%% instrumental waist sigma_T from the tau=0 correlation
g1 = @(p,rr) p(1)*exp(-rr/p(2)^2) + p(3);
q0 = lsqcurvefit(g1, [amp0, sigma0, offset0], r2v, reshape(STCorr_a(:,:,1),[],1), ...
                 [0 0.02 -10], [inf max(range(xv),range(yv)) 10], options);
sigma_T = abs(q0(2));
fprintf('   two-component fit: sigma_T (waist) = %.3f um\n', sigma_T);

%%% per-lag: G = A_trap*exp(-r^2/sigma_T^2) + A_mob*exp(-r^2/sigma_D^2) + offset
g2 = @(p,rr) p(1)*exp(-rr/sigma_T^2) + p(2)*exp(-rr/p(3)^2) + p(4);
Gaussb = zeros(ny,nx,nt);
A_trap = zeros(nt,1); amp_G = zeros(nt,1); s_G = zeros(nt,1); o_G = zeros(nt,1);
A_trap_err = zeros(nt,1); amp_G_err = zeros(nt,1); s_G_err = zeros(nt,1); o_G_err = zeros(nt,1);
prev = [amp0*0.5, amp0*0.5, sigma_T*2, offset0];
for i = 1:nt
    yv2 = reshape(STCorr_a(:,:,i),[],1);
    p0 = [max(prev(1),1e-7), max(prev(2),1e-7), min(max(prev(3),sigma_T*1.1),rmax*2), median(yv2)];
    try
        [p,~,resid,~,~,~,Jc] = lsqcurvefit(g2, p0, r2v, yv2, [0 0 sigma_T -10], [inf inf rmax*2 10], options);
        prev = p;
        try, ci = nlparci(p, resid, 'jacobian', Jc); pe = (ci(:,2)-ci(:,1))/2; catch, pe = nan(4,1); end
    catch
        p = prev; pe = nan(4,1);
    end
    A_trap(i)  = p(1); amp_G(i) = p(2); s_G(i) = p(3)^2; o_G(i) = p(4);
    A_trap_err(i) = pe(1); amp_G_err(i) = pe(2); s_G_err(i) = pe(3)*p(3)*2; o_G_err(i) = pe(4);
    Gaussb(:,:,i) = reshape(g2(p, r2v), ny, nx);
end

%%% map to the standard variable names -> the MOBILE curve is the "primary"
%%% (the two-component Gaussians are centred: x_G, y_G = 0; o_G is the real
%%%  fitted offset so the viewer's ylim on G_infty is valid)
x_G = zeros(nt,1); y_G = zeros(nt,1);
x_G_err = zeros(nt,1); y_G_err = zeros(nt,1);
R2_G = zeros(nt,1);
for i = 1:nt
    Gfit=Gaussb(:,:,i); Gdata=STCorr_a(:,:,i);
    R2_G(i)=1-sum((Gdata(:)-Gfit(:)).^2)/sum((Gdata(:)-mean(Gdata(:))).^2);
end

%%% trapped component (flat iMSD) - used by imsd_04_output
s_G_trap = sigma_T^2 * ones(nt,1);

spatiotemporalCorrelationViewer(time_vec, frame_time, px_size, a1, b1, TauLimit, ...
    STCorr_a(:,:,1:TauLimit), Gaussb, amp_G, s_G, o_G, R2_G, amp_G_err, s_G_err, o_G_err)
end
