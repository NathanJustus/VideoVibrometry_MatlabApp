function G = gaborFilterBank(n_scales,n_orients,r,c)
% generates gabor filters in the frequency domain for a specific image size
%
% Inputs
%   nscales     number of spatial scales in the filter bank
%   n_orients   number of irnetations in the filter bank
%   r           size of first dimension of desired image to be filtered
%   c           size of second dimension of desired image to be filtered
%
% Outputs
%   G           structure containing the gabor filter bank

% TODO: allow different aspect ratios
%       control aspect ratio based on number of orientations
%       input checking

% r and c and the rows and column sizes of the image to be filtered
% generate frequency vectors in rad/px from -pi to pi
if mod(c,2) == 0 % x-axis
    u = (-c/2:c/2-1)/c*2*pi;
else
    u = (-(c-1)/2:(c-1)/2)/c*2*pi;
end

if mod(r,2) == 0 % y-axis
    v = ((-r/2:r/2-1)/r)'*2*pi;
else
    v = ((-(r-1)/2:(r-1)/2)/r)'*2*pi;
end

% set bandwidth using sigma, aspect ratio is 1 in this case
sigma = sqrt(-(pi/3-pi/2)^2/log(1/2)); % sets overlap at half maximums
scale = 1:n_scales;
orientation = linspace(0,pi,n_orients+1);
orientation = orientation(1:end-1);
frequency = pi./(2.^(scale));

G = struct();
ind = 1;
for i = 1:length(scale)
    for j = 1:length(orientation)
        [u0,v0] = pol2cart(orientation(j),frequency(i));
        G(ind).FFTkernel = exp(-((u-u0).^2./sigma^2+(v-v0).^2./sigma^2));
        G(ind).orientation = rad2deg(orientation(j)); % deg
        G(ind).centerFrequency = frequency(i); % rad/px
        G(ind).wavelength = 2*pi/frequency(i); % px
        G(ind).sigma = sigma;
        G(ind).fx = u; % x frequency vector (rad/px)
        G(ind).fy = v'; % y frequency vector(rad/px)
        ind = ind+1;
    end
    sigma = sigma/2; % reduce the bandwidth by half at each scale
end