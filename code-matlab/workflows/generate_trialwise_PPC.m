% label = 'R117-2007-06-09-TT03_4'; type = 1;
label = 'R131-2007-09-04-TT10_5'; type = 2;
% 1 for MSN, 2 for FSI
cd('E:\Dropbox (Dartmouth College)\AnalysisResults\FieldTripResults\ft_trialwise_ppc');
fname = strcat(extractBefore(label,'-TT'),'_ft_spec.mat');
load(fname);
%%
if type == 2 % FSI Case
    fsi_labels  = od.label(od.cell_type == 2);
    fsi_labels = cellfun(@(x) extractBefore(x, '.t'), fsi_labels, 'UniformOutput', false);
    iC = find(cellfun(@(x) strcmp(label, x), fsi_labels));
    tw_ppc = od.fsi_res.near_spec{iC}.trial_wise_ppc;
else
    msn_labels  = od.label(od.cell_type == 1);
    msn_labels = cellfun(@(x) extractBefore(x, '.t'), msn_labels, 'UniformOutput', false);
    iC = find(cellfun(@(x) strcmp(label, x), msn_labels));
    tw_ppc = od.msn_res.near_spec{iC}.trial_wise_ppc;
end
% %Remove 0s
% for iR = 1:size(tw_ppc,1)
%    tw_ppc(iR,:) = normdata(tw_ppc(iR,:));
% end
figure;
subplot(5,1,1)
imagesc(tw_ppc(1,:);  
subplot(5,1,2)
imagesc((tw_ppc(2,:));
subplot(tw_ppc(1,:
 %%
 function ndata = normdata(data)
    if (sum(isnan(data)) == length(data))
        ndata = data;
    else
        maxd = max(data);
        mind = min(data);
        ndata = (data - mind)/(maxd-mind);
    end
end
