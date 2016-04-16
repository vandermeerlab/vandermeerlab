function [err,confMat] = DecodeErrorZ(cfg_in,P,z)
% function [err,confMat] = DecodeErrorZ(cfg_in,P,z)
%
% computes decoding error and confusion matrices based on decoded posterior 
% tsd P and true value tsd z
%
% NOTE: assumes that P.data is [nBins x nSamples]
%
% MvdM 2016-02-23

cfg_def = [];
cfg_def.method = 'max'; % {'max','com','p(true'}

cfg = ProcessConfig(cfg_def,cfg_in);

nBins = size(P.data,1);

% for each time sample in P, find true Z
true_z_idx = nearest_idx3(P.tvec,z.tvec);
true_z = z.data(true_z_idx);

% but only keep values of P with no NaNs
keep_idx = ~isnan(nansum(P.data));
true_z = true_z(keep_idx);
P.data = P.data(:,keep_idx);

switch cfg.method
    case 'max'
        [~,decoded_z] = max(P.data);
        err = abs(true_z-decoded_z);
        
        % also make confusion matrix based on decoded_z
        confMat.thr = confusionmat(true_z,decoded_z,'order',1:nBins);
        confMat.thr = confMat.thr./repmat(nansum(confMat.thr,2),[1 nBins]); % normalize
                
end

% confusion matrix with full posterior
confMat.full = nan(nBins);
for iB = 1:nBins
    this_z_idx = find(true_z == iB);
    if length(this_z_idx) > 1
        confMat.full(:,iB) = nanmean(P.data(:,this_z_idx),2);
    elseif length(this_z_idx) == 1 % if only one sample, don't nanmean
        confMat.full(:,iB) = P.data(:,this_z_idx);
    end
end

temp_err = nan(size(P.tvec));
temp_err(keep_idx) = err;
err = tsd(P.tvec,temp_err);