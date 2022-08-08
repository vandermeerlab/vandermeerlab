function [data] = processFP(params, data)
%Process Fiber Photometry
%
%   [data] = processFP(data,params)
%
%   Description: This function is designed to process fiber photometry data
%   for the lab. The function performs demodulation (if selected),
%   filtering, baselining, and downsampling for all photometry traces in
%   the recording. The parameters for the analysis are found in the params
%   structure, which is created from a user-created scripted based on the
%   processParam.m file.
%
%   Input:
%   - data - A data structure specific to the Tritsch Lab. Created using
%   the convertH5_FP script
%   - params - A structure created from a variant of the processParams
%   script
%
%   Output:
%   - data - Updated data structure containing processed data
%
%   Author: Pratik Mistry 2019

nAcq = length(data.acq);
lpCut = params.FP.lpCut; filtOrder = params.FP.filtOrder;

dsRate = params.dsRate;

interpType = params.FP.interpType;
fitType = params.FP.fitType; winSize = params.FP.winSize;
winOv = params.FP.winOv;
basePrc = params.FP.basePrc;

rawFs = data.acq.Fs;
Fs = rawFs/dsRate;
data.final.Fs = Fs;

for n = 1:nAcq
    L = size(data.acq(n).time,1);
    nFP = length(data.acq(n).FP);
    data.final(n).FP = cell(nFP,1);
    data.final(n).nbFP = cell(nFP,1);
    data.final(n).FPbaseline = cell(nFP,1);
    L = length(1:dsRate:L);
    for x = 1:nFP
        rawFP = data.acq(n).FP{x};
        nbFP = filterFP(rawFP,rawFs,lpCut,filtOrder,'lowpass');
        nbFP = downsample(nbFP, dsRate);
        [FP,baseline] = baselineFP(nbFP,interpType,fitType,basePrc,winSize,winOv,Fs);
        data.final(n).FP{x} = FP;
        data.final(n).nbFP{x} = nbFP;
        data.final(n).FPbaseline{x} = baseline;
    end
    timeVec = [1:L]/Fs;
    data.final(n).time = timeVec';
end