function WriteFig(varargin)
% function WriteFig(fn): writes current figure to fn
% function WriteFig(fh, fn): writes figure with handle fh to fn
% function WriteFig(fh, fn, pngOnly): as above but only png output

pngOnly = 0;
switch nargin
    
    case 1
        fn = varargin{1};
        fh = gcf;
    case 2
        fh = varargin{1};
        fn = varargin{2};
    case 3
        fh = varargin{1};
        fn = varargin{2};
        pngOnly = 1;
    otherwise
        error('RTFM');
end
        

figure(fh);
print('-dpng','-r300', cat(2,fn,'.png'));
if ~pngOnly
    %print('-dill', cat(2,fn,'.ai'));
    print('-depsc', cat(2,fn,'.eps'));
    saveas(gcf, cat(2, fn, '.fig'), 'fig');
end
