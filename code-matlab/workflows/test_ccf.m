%%ccf testing script
st1 = rand(100,1);
st2 = st1 + 0.15;
s1 = ts({st1});
s2 = ts({st2});
cfg1.smooth = 0;
cfg1.binsize = 0.01;
cfg2.smooth = 1;
cfg2.xcorr = 'coeff';
cfg2.binsize = 0.01;

[q1,t1] = ccf1(cfg1,s1.t{1},s1.t{1});
[q2,t2] = ccf1(cfg2,s1.t{1},s1.t{1});
subplot(2,1,1)
plot(t1,q1);
title('CCF with raw spike count version', 'FontSize', 15)
subplot(2,1,2)
plot(t2,q2);
title('CCF with spike density but no smoothing', 'FontSize', 15)