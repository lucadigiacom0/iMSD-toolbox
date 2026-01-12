%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% spatial autocorrelation function %%%%%%%%%

%%% fitting parameters
startGxy=[0.4 0.1 0];
lGxy=[0 0 -1];
uGxy=[inf inf inf];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = (-R/2:px_size:(R-px_size)/2);
b = (-R/2:px_size:(R-px_size)/2);

disp ('     calculating g(\xi, eta, 0)_n ...')

%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:N
    fprintf('     computing g_fluct(xi, eta, 0)_%d\n', k);
    I = double(A(:,:,k));
    meanI = mean(I(:));
    dI = I - meanI;
    F = fft2(dI);
    ACF = real(ifft2(F .* conj(F)));
    ACF = fftshift(ACF);
    Z(:,:,k) = ACF / (numel(I) * meanI^2);
end

Z_mean=mean(Z,3);

midX = ceil(size(Z,1)/2);
midY = ceil(size(Z,2)/2);
zeroLag(k) = Z(midX,midY,k);

[MX, MY, numFrames] = size(Z);
cx = round((MX + 1)/2);
cy = round((MY + 1)/2);
for k = 1:numFrames
    Zk = Z(:,:,k);
    mask = true(MX, MY);
    mask(cx, cy) = false;
    maxVal = max(Zk(mask), [], 'all');
    Z(cx, cy, k) = maxVal;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%5

a1=linspace(-px_size*SpatialLimit/2, px_size*SpatialLimit/2, SpatialLimit);
b1 = a1;
Z2=Z(P/2+1+int64(a1/px_size),P/2+1+int64(a1/px_size),:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

%%% fitting spatial autocorrelation function
for i=1:N

[X Y]=meshgrid(a1,a1);
xy(:,:,1)=X;
xy(:,:,2)=Y;
options=optimset('display','off');
modelFunGxy=@(pgxy,xy)(pgxy(1)*exp(-((xy(:,:,1)).^2+(xy(:,:,2)).^2)./pgxy(2).^2)+pgxy(3));
%options = optimset('Display','off');
 Gxy=lsqcurvefit(modelFunGxy,startGxy,xy,Z2(:,:,i),lGxy,uGxy,options);
 Gaussxy(:,:,i)=modelFunGxy(Gxy,xy);
 a_Gxy(i,1)=Gxy(1);
 s_Gxy(i,1)=Gxy(2);
 o_Gxy(i,1)=Gxy(3);
 
 ydata = Z2(:,:,i);
 yfit = Gaussxy(:,:,i);
 SSres = sum((ydata(:) - yfit(:)).^2);            % residual sum of squares
 SStot = sum((ydata(:) - mean(ydata(:))).^2);     % total sum of squares
 R2 = 1 - SSres / SStot;
 R2_Gxy(i,1) = R2;

 Gxy_x(:,i)= Gaussxy(:,length(a1)/2+1,i);
 Gxy_y(:,i)= Gaussxy(length(a1)/2+1,:,i);
 
end
Z1=mean(Z2,3);
for i=1:N
    Gx(:,i)= Z2(:,length(a1)/2+1,i);
    Gy(:,i)= Z2(length(a1)/2+1,:,i);
end

omega2=mean(s_Gxy,1);
err_omega2=std(s_Gxy,1);
omega=omega2;
err_omega=err_omega2;

spatialCorrelationViewer(a1,b1,Z2,Gaussxy,px_size,a_Gxy,s_Gxy,o_Gxy,R2_Gxy);