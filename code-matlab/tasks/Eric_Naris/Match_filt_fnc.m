function filt_out = Match_filt_fnc(data_in, span)
%% Match_smooth_fnc: a fast alternative to the built in 'smooth' funciton


% a_s = smooth(data_in,span);
fb = ones(span,1)./span;
fa = zeros(span,1)./span; fa(1) = 1;
filt_out = filtfilt(fb,fa,data_in);

% figure;
% plot(data_in); hold on;
% plot(a_s,'r');
% plot(filt_out,'g');
