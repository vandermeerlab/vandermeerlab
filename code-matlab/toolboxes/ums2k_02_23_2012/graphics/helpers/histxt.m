function [counts,t_inds,x_inds] = histxt(x, varargin)
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% histxt - 2D histogram of set of waveforms
%
% Usage:
%    [counts,t_inds,x_inds] = histxt(x, varargin)
%
% Description:  
%   COUNTS = HISTXT(X,D), where X is an M x T matrix returns a D x T
%   matrix COUNTS in which each column contains the histogrammed (with D
%   bins) values from the corresponding column in X.  If D is not
%   specified (or is the empty matrix), it defaults to 100.
%
%   [COUNTS,T_INDS,X_INDS] = HISTXT(X,D) returns the column indices and
%   bin centers so that the density can be visualized with
%   IMAGESC(T_INDS,X_INDS,COUNTS).
%
%   [...] = HISTXT(...,'log') uses the log of the counts (0's yield -Inf).
%
%   HISTXT(...) without output arguments produces an image of the counts.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Parse Inputs %%%%%%%%%%%%%%%%%%%%%%%%%%%%
D = 100;  logflag = 0;                        % defaults

[M,T] = size(x);
if (~isnumeric(x) || ndims(x)~=2)
	error('First input argument must be a numeric matrix.');
end

if (length(varargin) > 0)
	tail = varargin{end};
	if (ischar(tail) && strcmpi(tail,'log'))  % If the last arg was 'log' ...
		varargin = varargin(1:(end-1));       % ... chomp it and set a flag.
		logflag = 1;
	end
	if (length(varargin) > 0)
		tail = varargin{end};
		if (isnumeric(tail) && length(tail)==1)   % If next to last arg was a scalar, ...
			varargin = varargin(1:end-1);         % ... chomp it and set the bin count.
			D = tail;
		end
	end
	if (length(varargin) > 0),  error('Unknown syntax.');  end;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%% Rescale Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scale the data
[x,oldmin,oldmax] = rescale(x,1,D);  x = round(x);

% Make bin centers/column indices 
x_inds = linspace(oldmin,oldmax,D);
t_inds = 1:T;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Histogram %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
counts = hist(x,1:D);


if (logflag)
    old = warning('off');   counts = log(counts);   warning(old);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Graphics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (nargout == 0)
     imagesc(t_inds, x_inds, counts);   axis xy;
     clear counts t_inds x_inds  % clear these so nothing is dumped to output
end
