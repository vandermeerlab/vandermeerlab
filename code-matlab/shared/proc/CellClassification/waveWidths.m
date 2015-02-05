function [pw,vw] = waveWidths(mWV)
% function [pw,vw] = waveWidths(mWV)
%
% compute peak and valley widths (in ms) at half amplitude (relative to
% zero) of highest amplitude waveform
%
% assumes (32x4) mWV input from wv file (use Create_CQ_File.m)
%
% MvdM 09

% example wv file
% load C:\data\R117\R117_neurophys\R117-2007-06-01\R117-2007-06-01-TT01-1-wv.mat

debug = 0;

% find biggest wv
bigAmpl = 0;
bigCh = 0;
for iCh = 1:4
   
    ampl = max(mWV(:,iCh))-min(mWV(:,iCh));
    if  ampl > bigAmpl
        bigAmpl = ampl;
        bigCh = iCh;
    end
    
end
wv = mWV(:,bigCh);

% find peak & valley; assume one peak & take valley that follows it

[pval,pind] = max(wv);
[vval,vind] = min(wv(pind:end)); vind = vind + pind - 1;

if vval == 32
    plot(mWV(:,bigCh));
    title('problematic one');
    pause;
end

% find half-widths, peak
half = pval / 2;
temp = wv(1:pind) - half; temp(temp < 0) = NaN;
[vpre_peakH,ipre_peakH] = min(temp); % first sample above half, now check which one is closer
if ipre_peakH < 1
    if abs(wv(ipre_peakH) - half) > abs(wv(ipre_peakH - 1) - half)
        ipre_peakH = ipre_peakH - 1;
    end
end

temp = wv(pind:end) - half; temp(temp < 0) = NaN;
[vpost_peakH,ipost_peakH] = min(temp); % first sample above half, now check which one is closer
ipost_peakH = ipost_peakH + pind - 1;
if abs(wv(ipost_peakH) - half) > abs(wv(ipost_peakH + 1) - half)
   ipost_peakH = ipost_peakH + 1; 
end

% find half-widths, valley
wvi = -wv;
half = -vval / 2;
temp = wvi(pind:vind) - half; temp(temp < 0) = NaN;
[vpre_valH,ipre_valH] = min(temp); % first sample below half, now check which one is closer
ipre_valH = ipre_valH + pind - 1;
if abs(wvi(ipre_valH) - half) > abs(wvi(ipre_valH - 1) - half)
   ipre_valH = ipre_valH - 1; 
end

temp = wvi(vind:end) - half; temp(temp < 0) = NaN;
[vpost_valH,ipost_valH] = min(temp); % first sample above half, now check which one is closer
ipost_valH = ipost_valH + vind - 1;

if ipost_valH < 32
    if abs(wvi(ipost_valH) - half) > abs(wvi(ipost_valH + 1) - half)
        ipost_valH = ipost_valH + 1;
    end
end

% check it looks ok
if debug
    plot(mWV(:,bigCh));
    hold on;
    plot(pind,pval,'k.');
    plot(vind,vval,'k.');
    plot(ipre_peakH,wv(ipre_peakH),'r.');
    plot(ipost_peakH,wv(ipost_peakH),'r.');
    plot(ipre_valH,wv(ipre_valH),'g.');
    plot(ipost_valH,wv(ipost_valH),'g.');
    pause; clf;
end

pw = (ipost_peakH - ipre_peakH)*(1/32);
vw = (ipost_valH - ipre_valH)*(1/32);

if pw < 0 | vw < 0
   error('Negative widths!'); 
end



