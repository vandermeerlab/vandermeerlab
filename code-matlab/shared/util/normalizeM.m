function out = normalizeM(in, varargin)
% function out = normalizeM(in, method)
%
% row-wise normalization of input matrix

switch length(varargin)

    case 0
        method = 'zscore'; % default
    
    case 1
        method = varargin{1};
        
    otherwise
            error('Not yet supported');
            
end

nCols = size(in, 2);
switch method
    
    case 'zscore'
        
        m = nanmean(in, 2); m = repmat(m, [1 nCols]);
        s = nanstd(in, [], 2); s = repmat(s, [1 nCols]);
        out = (in - m) ./ s;
        
    otherwise
        
end


