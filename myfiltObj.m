function filtObj = myfiltObj(cutoffFrequencies, nPoles, samplingFrequency)
%% Returns a discrete-time filter object for a Butterworth bandpass filter
%
% File History:
%   Butterworth Bandpass filter designed using FDESIGN.BANDPASS.
%   M-File generated by MATLAB(R) 7.6 and the Signal Processing Toolbox 6.9.
%   Generated on: 31-May-2011 13:51:44

% Construct an FDESIGN object and call its BUTTER method.
h  = fdesign.bandpass('N,Fc1,Fc2', nPoles, cutoffFrequencies(1), ...
                                   cutoffFrequencies(2), samplingFrequency);
filtObj = design(h, 'butter');