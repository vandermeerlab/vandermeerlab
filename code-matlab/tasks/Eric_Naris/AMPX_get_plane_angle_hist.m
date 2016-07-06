function AMPX_get_plane_angle_hist(data_out)


l_angle = data_out.lg.power.stats.theta;
lr_angle = data_out.lg_ran.power.stats.theta;
h_angle = data_out.hg.power.stats.theta;
hr_angle = data_out.hg_ran.power.stats.theta;


[n1, p1] = hist(l_angle, 25);
[n2, p2] = hist(lr_angle, 25);
[n3, p3] = hist(h_angle, 25);
[n4, p4] = hist(hr_angle, 25);

%%
figure(10000)
subplot(211)
hold on
h2 = bar(p2, n2, 'r', 'barwidth', .9, 'edgecolor', 'none');
h = bar(p1,n1,  'b', 'barwidth', .9, 'edgecolor', 'none');

legend([h h2], {'low', 'low ran'})

subplot(212)
hold on
h4 = bar(p4, n4, 'r', 'barwidth', .9, 'edgecolor', 'none');
h3 = bar(p3, n3, 'b', 'barwidth', .9, 'edgecolor', 'none');

legend([h3 h4], {'high', 'high ran'})

%%
[Hl, Pl] = ttest2(l_angle, lr_angle)

[Hh, Ph] = ttest2(h_angle, hr_angle)
