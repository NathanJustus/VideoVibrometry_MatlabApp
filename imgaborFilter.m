function r = imgaborFilter(I,G)
% takes input filter bank generated using gaborFilterBank() and image I
% (double representation) and outputs (complex) array r where the third
% dimension is each orientation and scale corresponding to the elements in
% filterbank G

I = fftshift(fft2(I));
r = complex(zeros([size(I) length(G)]));
for i = 1:length(G)
    R = G(i).FFTkernel.*I;
    r(:,:,i) = (ifft2(ifftshift(R)));
end