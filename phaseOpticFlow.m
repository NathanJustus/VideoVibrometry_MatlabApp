function [Vx,Vy,t] = phaseOpticFlow(vr,scale,crop,t0,nf,fs,nscales)

% notes:
%   add in setting to do single image frames
%   add in setting to control number of scales and orientations

t = (0:nf-1)/fs; % time vector

vr.CurrentTime = t0;

I = im2double(im2gray(vr.readFrame));
I = imcrop(I,crop);
I = imresize(I,scale);

[r,c] = size(I);
G = gaborFilterBank(nscales,4,r,c); % 3 scales and 4 orientations
r1 = imgaborFilter(I,G);

w = [1 4 6 4 1]/16; % weighting kernel for pixel neighborhood

M11 = zeros(size(r1));
M12 = zeros(size(r1));
M22 = zeros(size(r1));
B1 = zeros(size(r1));
B2 = zeros(size(r1));
drdx = zeros(size(r1));
drdy = zeros(size(r1));
Vx = zeros([size(I) nf]);
Vy = zeros([size(I) nf]);

H = waitbar(0,'Please wait...','name','Video Processing Progress');
for i = 1:nf
    I = im2double(im2gray(vr.readFrame));
    I = imcrop(I,crop);
    I = imresize(I,scale);
    
    r2 = imgaborFilter(I,G); % apply filter bank
    dphidt = angle(r2.*conj(r1))*fs; % time derivative using complex form, automatically unwraps
    dphidt = ampWeightedBlur(dphidt,abs(r2),2); % amplitude weighted spatial blur
    
    % spatial gradients of phase
    % see the van hulle paper that uses spatial gradients rather than the
    % spatio-temporal gradients seen in fleet & jepson
    for j = 1:length(G)
        % derivative of complex images, see phase-based disparity paper for
        % the demodulated kernels using spatial frequency and a cos/sin
        % filter for bandpass sampling
        k0 = G(j).centerFrequency;
        h = [-exp(-1i*2*k0) 8*exp(-1i*k0) 0 -8*exp(1i*k0) exp(1i*2*k0)]/12;
        drdx(:,:,j) = conv2(r2(:,:,j),h,'same') + 1i*k0*r2(:,:,j);
        drdy(:,:,j) = conv2(r2(:,:,j),h.','same') + 1i*k0*r2(:,:,j);
    end
    % spatial derivatives of phase using phase difference
    dphidx = imag(drdx.*conj(r2))./(abs(r2).^2);
    dphidy = imag(drdy.*conj(r2))./(abs(r2).^2);
    
    % set up least squares problem for each pixel for flow estimation
    % NOTE:
    % attempt to do the actual LS-solve for a 5x5 neighborhood rather than
    % setting derivatives to zero
    % solve the system grad(phi)*[u v] = -dphi/dt
    dpdx2 = dphidx.^2;
    dpdxdpdy = dphidx.*dphidy;
    dpdy2 = dphidy.^2;
    dpdxdpdt = dphidx.*dphidt;
    dpdydpdt = dphidy.*dphidt;
    % generate pixel neighborhood summations
    for j = 1:length(G)
        % convolutions act like summing over the pixel neighborhoods
        % weighting kernel is used to set neighborhood
        M11(:,:,j) = conv2(w,w,dpdx2(:,:,j),'same');
        M12(:,:,j) = conv2(w,w,dpdxdpdy(:,:,j),'same');
        M22(:,:,j) = conv2(w,w,dpdy2(:,:,j),'same');
        B1(:,:,j) = conv2(w,w,dpdxdpdt(:,:,j),'same');
        B2(:,:,j) = conv2(w,w,dpdydpdt(:,:,j),'same');
    end
    % sum up the least squares entries over all orientations and scales
    m11 = sum(M11,3);
    m12 = sum(M12,3);
    m22 = sum(M22,3);
    b1 = -sum(B1,3);
    b2 = -sum(B2,3);
    % apply cramers rule in batch for the 2x2 system to obtain pixel
    % velocities
    Vx(:,:,i) = (b1.*m22-m12.*b2)./(m11.*m22-m12.^2);
    Vy(:,:,i) = (m11.*b2-b1.*m12)./(m11.*m22-m12.^2);
    
    % shift registers
    r1 = r2;
    
    % progress bar update
    progress = round(i/nf*100);
    waitbar(progress/100,H,[num2str(progress) '% complete'])
end


close(H)

