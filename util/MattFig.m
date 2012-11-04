function MattFig(figno,fn,varargin)
% function MattFig(figno,fn)

pngOnly = 1;
extract_varargin;

figure(figno);
print('-dpng','-r300',cat(2,fn,'.png'));
if ~pngOnly
    print('-dill',cat(2,fn,'.ai'));
    print('-depsc',cat(2,fn,'.eps'));
    saveas(gcf,cat(2,fn,'.fig'),'fig');
end
