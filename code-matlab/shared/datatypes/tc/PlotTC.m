function [] = tcplot(cfg_in,tc_in)
%% TCPLOT Tuning Curve Plot
%
%   Plots the tuning curves from a tc object
%
%   INPUTS:
%       cfg_in: input cfg file
%       tc_in: input tc
%
%   CFG OPTIONS:
%       cfg.display = 'separate' % 'together'
%
%   NOTES:
%       Currently only plots 1D tuning curves.
%
% youkitan 2014-11-05
%% set defaults and process cfg
cfg_def.display = 'separate';

mfun = mfilename;
cfg = ProcessConfig(cfg_def,cfg_in,mfun);

%% create plot layout
numPF = length(tc_in.firing_rates);
[nsp,~] = numSubplots(numPF);

%% plot
cases = {'separate','together'};
assert(any(strcmp(cfg.display,cases)),...
    '"%s" is an invalid plot type!',cfg.display)

switch cfg.display
    case 'separate'
        indPF = 1;
        hold on;
        for i = 1:(nsp(1)*nsp(2))
            if indPF <= numPF
                subplot(nsp(1),nsp(2),indPF);
                plot(tc_in.stimulus{indPF},tc_in.firing_rates{indPF});
                title(['Cell ',num2str(i)],'Fontsize',8);
                xlabel('Position (cm)','Fontsize',8);
                ylabel('Firing Rate (Hz)','Fontsize',8);
                set(gca,'Fontsize',6)
                indPF = indPF + 1;
            end
        end%iterate rows

    case 'together'
        plotcolors = {'c','r','g','b','k'};
        hold on;
        for i = 1:numPF
            cval = mod(i,5) + 1;
            plot(tc_in.stimulus{i},tc_in.firing_rates{i},'Color',plotcolors{cval});
            xlabel('Position (cm)','Fontsize',8);
            ylabel('Firing Rate (Hz)','Fontsize',8);
            set(gca,'Fontsize',6)
        end %iterate rows
end %end switch


end