function [cleanSignal,filt] = filterFP(rawSignal,sampleRate,cutoff,order,filtType)
%filterFP - Filter raw photometry signals
%
%   [cleanSignal,filt] = filterFP(rawSignal,sampleRate,cutoff,order,filtType)
%   [cleanSignal,filt,adjustedSamplingRate] = filterFP(rawSignal,sampleRate,cutoff,order,filtType,downsamplefactor)
%
%   Description: This function allows the user to implement a lowpass or
%   highpass butterworth filter. A butterworth filter was choosen because the passband and stopband are
%   maximally flat. However, the roll-off is broad; you must increase frequency order if you want
%   a steeper roll-off. This function will also allow the user to visualize the bode
%   plot of the designed filter. Optionally, the function will also allow the
%   user to downsample the signal afterwards.
%
%
%   Input:
%   rawSignal = Unadjusted photometry signal
%   sampleRate = sampling rate used to acquire photometry signal
%   cutoff = Cut off frequency desired by the user
%   order = Desired filter order -- Note: Higher order provides a steeper
%   roll-off
%   filtType = Choose between 'lowpass' or 'highpass'
%   downsamplefactor = Optional final input that allows user to downsample
%   the signal. This optional will also provide the user with the downsampled
%   sampling rate as an additional output
%
%   Output:
%   cleanSignal = filtered signal
%   filt = filter object to present filter/system properties
% 
%   Author: Pratik Mistry

if size(rawSignal,1) == 1
    rawSignal = rawSignal';
end

switch filtType
    case 'lowpass'
        filt=designfilt('lowpassiir','FilterOrder',order,'HalfPowerFrequency',cutoff,'SampleRate',sampleRate,'DesignMethod','butter');
        %fvtool(filt);
    case 'highpass'
        filt=designfilt('highpassiir','FilterOrder',order,'HalfPowerFrequency',cutoff,'SampleRate',sampleRate,'DesignMethod','butter');
        %fvtool(filt);
end
cleanSignal=filtfilt(filt,rawSignal);
end
