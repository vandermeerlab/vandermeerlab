function [slopes,intercepts] = circ_regress(cfg_in,data_in)
% function [slopes,intercepts] = circ_regress(cfg_in,data_in)
%
% moving window circular regression
%
% method: incrementally adds values to data_in (from cfg.offset_delta to
% 2*pi-cfg_offset_delta) and finds the best linear regression fit
%
% INPUTS
%
% - data_in: [N x 1] circular data, domain [0,2*pi>
%
% OUTPUTS
%
% - slopes: linear regression fit slopes
% - intercepts: linear regression fit intercepts (at the first point in the
%   window)
% 
% CONFIGS
%
% cfg_def.window_size = 11; % size of moving window (in samples)
%    NOTE: this value must be odd
% cfg_def.offset_delta = pi/64; % step size
%
% MvdM 2018

cfg_def = [];
cfg_def.window_size = 11;
cfg_def.offset_delta = pi/64;

cfg = ProcessConfig(cfg_def,cfg_in);

if mod(cfg.window_size,2) ~= 1
    error('cfg.window_size must be odd');
end

data_in = data_in(:);

if (max(data_in) >= 2*pi) || (min(data_in) < 0)
    error('data_in must contain values within the interval [0,2*pi>');
end

nPoints = length(data_in);
hw = floor(cfg.window_size / 2); % half window length

fv = 0:cfg.offset_delta:2*pi-cfg.offset_delta; % offsets to fit
nFits = length(fv);

for iP = nPoints:-1:1
    
    this_window(1) = max(1,iP-hw);
    this_window(2) = min(nPoints,iP+hw);
    
    this_data = data_in(this_window(1):this_window(2)); % data to fit on this iteration
    
    xval = this_window(1):this_window(2); xval = xval(:);

    for iF = nFits:-1:1
        
        yval = modpi(this_data + fv(iF));
        
        [coeff,s] = fast_polyfit(xval,yval,1); % obtain linear fit
        
        this_fit_slope(iF) = coeff(1); this_fit_int(iF) = coeff(2);
        this_fit_r(iF) = s.normr;
        
    end % of data offsets
   
    % find best fit and store
    [~,keep_idx] = min(this_fit_r);
    slopes(iP) = this_fit_slope(keep_idx);
    intercepts(iP) = modpi(this_fit_int(keep_idx) - fv(keep_idx));
    
end % of data points


function out = modpi(in)

out = mod(in,2*pi);