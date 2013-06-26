% GABORCONVOLVE - function for convolving image with log-Gabor filters
%
% Usage: [EO, BP] = gaborconvolve(im,  nscale, norient, minWaveLength, mult, ...
%			    sigmaOnf, Lnorm, feedback)
%
% Arguments:
% The convolutions are done via the FFT.  Many of the parameters relate 
% to the specification of the filters in the frequency plane.  
%
%   Variable       Suggested   Description
%   name           value
%  ----------------------------------------------------------
%    im                        Image to be convolved.
%    nscale          = 4;      Number of wavelet scales.
%    norient         = 6;      Number of filter orientations.
%    minWaveLength   = 3;      Wavelength of smallest scale filter.
%    mult            = 2;      Scaling factor between successive filters.
%    sigmaOnf        = 0.65;   Ratio of the standard deviation of the
%                              Gaussian describing the log Gabor filter's
%                              transfer function in the frequency domain
%                              to the filter center frequency. 
%    Lnorm            0        Optional integer indicating what norm the
%                              filters should be normalized to.  A value of 1
%                              will produce filters with the same L1 norm, 2
%                              will produce filters with matching L2
%                              norm. the default value of 0 results in no
%                              normalization (the filters have unit height
%                              Gaussian transfer functions on a log frequency
%                              scale) 
%    feedback         0/1      Optional parameter.  If set to 1 a message
%                              indicating which orientation is being
%                              processed is printed on the screen.
%
% Returns:
%
%   EO - 2D cell array of complex valued convolution results
%
%        EO{s,o} = convolution result for scale s and orientation o.
%        The real part is the result of convolving with the even
%        symmetric filter, the imaginary part is the result from
%        convolution with the odd symmetric filter.
%
%        Hence:
%        abs(EO{s,o}) returns the magnitude of the convolution over the
%                     image at scale s and orientation o.
%        angle(EO{s,o}) returns the phase angles.
%
%   BP - Cell array of bandpass images corresponding to each scale s.
%   
%
% Notes on filter settings to obtain even coverage of the spectrum
% dthetaOnSigma 1.5
% sigmaOnf  .85   mult 1.3
% sigmaOnf  .75   mult 1.6     (bandwidth ~1 octave)
% sigmaOnf  .65   mult 2.1
% sigmaOnf  .55   mult 3       (bandwidth ~2 octaves)
%                                                       
% For maximum speed the input image should be square and have a 
% size that is a power of 2, but the code will operate on images
% of arbitrary size.  
%
%
% The determination of mult given sigmaOnf is entirely empirical
% What I do is plot out the sum of the filters in the frequency domain
% and see how even the coverage of the spectrum is.
% If there are concentric 'gaps' in the spectrum one needs to
% reduce mult and/or reduce sigmaOnf (which increases filter bandwidth)
%
% If there are 'gaps' radiating outwards then one needs to reduce
% dthetaOnSigma (increasing angular bandwidth of the filters)
%

% For details of log-Gabor filters see: 
% D. J. Field, "Relations Between the Statistics of Natural Images and the
% Response Properties of Cortical Cells", Journal of The Optical Society of
% America A, Vol 4, No. 12, December 1987. pp 2379-2394

% Copyright (c) 2001-2010 Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% http://www.csse.uwa.edu.au/
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% May   2001 - Original version
% April 2010 - Reworked to tidy things up. Return of bandpass images added.

function [EO, BP, S] = gaborconvolve_3d(im, nscale, num_theta, num_phi, minWaveLength, mult, sigmaOnf, dThetaSigma, dPhiSigma)

% Pre compute data
thetaSigma = pi / num_theta / dThetaSigma;
phiSigma = (2 * pi) / num_phi / dPhiSigma;

EO = cell(nscale, num_theta, num_phi);          
imagefft = fftn(im);

S = 0;
if ~ndims(im) == 3
    error('3D image required');
end

[height, width, depth] = size(im);					

% Pre-compute some stuff to speed up filter construction
% Set up X and Y matrices with ranges normalised to +/- 0.5
% The following code adjusts things appropriately for odd and even values
% of rows and columns.
xrange = linspace(-0.5, 0.5, width);
yrange = linspace(-0.5, 0.5, height);
zrange = linspace(-0.5, 0.5, depth);

[x,y,z] = meshgrid(xrange, yrange, zrange);

radius = sqrt(x.^2 + y.^2 + z.^2);       % Matrix values contain *normalised* radius from centre.
theta = atan2(y, x);               % Matrix values contain polar angle.
m_ab = abs(mean(radius(:)));
phi = acos(z ./ (radius + m_ab));

radius = ifftshift(radius);       % Quadrant shift radius and theta so that filters
theta  = ifftshift(theta);
phi  = ifftshift(phi);        % are constructed with 0 frequency at the corners.
radius(1, 1, :) = 1;                 % Get rid of the 0 radius value at the 0
                                  % frequency point (now at top-left corner)
                                  % so that taking the log of the radius will 
                                  % not cause trouble.

sintheta = sin(theta);
costheta = cos(theta);
sinphi = sin(phi);
cosphi = cos(phi);
clear phi; clear theta;    % save a little memory

% Filters are constructed in terms of two components.
% 1) The radial component, which controls the frequency band that the filter
%    responds to
% 2) The angular component, which controls the orientation that the filter
%    responds to.
% The two components are multiplied together to construct the overall filter.

% Construct the radial filter components...
% First construct a low-pass filter that is as large as possible, yet falls
% away to zero at the boundaries.  All log Gabor filters are multiplied by
% this to ensure no extra frequencies at the 'corners' of the FFT are
% incorporated. This keeps the overall norm of each filter not too dissimilar.
lp = lowpassfilter_3d([height, width, depth], .45, 15);   % Radius .45, 'sharpness' 15

logGabor = cell(1, nscale);
BP = cell(nscale,1);

for s = 1:nscale
    wavelength = minWaveLength * mult^(s - 1);
    fo = 1.0 / wavelength;                  % Centre frequency of filter.
    logGabor{s} = exp((-(log(radius / fo)) .^ 2) / (2 * log(sigmaOnf)^2));  
    logGabor{s} = logGabor{s} .* lp;        % Apply low-pass filter
    logGabor{s}(1, 1, :) = 0;                  % Set the value at the 0 frequency point of the filter back to zero (undo the radius fudge).      
    
    BP{s} = ifft2(imagefft .* logGabor{s});
end

% The main loop - For each orientation.
for e = 1:num_theta
    eleAngle = (e - 1) * pi / num_theta;     % Calculate filter angle.

    % Pre-compute filter data specific to this orientation
    % For each point in the filter matrix calculate the angular distance from the
    % specified filter orientation.  To overcome the angular wrap-around problem
    % sine difference and cosine difference values are first computed and then
    % the atan2 function is used to determine angular distance.
    dts = sintheta * cos(eleAngle) - costheta * sin(eleAngle);     % Difference in sine.
    dtc = costheta * cos(eleAngle) + sintheta * sin(eleAngle);     % Difference in cosine.
    dtheta = abs(atan2(dts, dtc));                                 % Absolute angular distance.
    
    for a = 1:num_phi  
        aziAngle = a * 2 * pi / num_phi;
        dps = sinphi * cos(aziAngle) - cosphi * sin(aziAngle);     % Difference in sine.
        dpc = cosphi * cos(aziAngle) + sinphi * sin(aziAngle);     % Difference in cosine.
        dphi = abs(atan2(dps, dpc));                               % Absolute angular distance.

        spread = exp(((-dphi .^ 2) / (2 * phiSigma^2)) + ((-dtheta .^ 2) / (2 * thetaSigma^2)));
        
        % Reset filter wavelength.
        wavelength = minWaveLength;  
        
        for s = 1:nscale,
            filter = fftshift(logGabor{s} .* spread);
            S = S + filter .* conj(filter);
            EO{s, e, a} = ifft2(imagefft .* filter);

            % Finally calculate Wavelength of next filter and process the next scale
            wavelength = wavelength * mult;       
        end
    end
end
% TODO: Is this correct?? Why?
% S = flipdim(S, 3);