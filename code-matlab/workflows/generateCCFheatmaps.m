%% Script to plot ccfs as heatmaps ordered by earliest occurring peak
close all;
cd('/Users/manishm/Work/vanDerMeerLab/RandomVStrDataAnalysis/ccf_results/Smooth');
rats = {'R117','R119','R131','R132'};
for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat('*',curRat,'*od.mat');
    ofiles = dir(searchString);
    for jdx = 1:length(ofiles)
        load(ofiles(jdx).name);
        fig = figure;
        d1 = od.cx1;
        d2 = od.cx2;
        d3 = od.cx3;
        msn_count = length(od.l1);
        fsi_count = length(od.l2);
        %Table for MSN-MSN
        if ~isempty(d1)
            subplot(3,1,1)
            ne = (msn_count*(msn_count-1))/2;
            m2m_label1 = cell(ne,1);
            m2m_label2 = cell(ne,1);
            nd1 = cellfun(@normdata, d1, 'UniformOutput', false);
            k = 1;
            for i = 1:msn_count-1
                for j = i+1:msn_count
                    m2m_label1{k} = od.l1{i};
                    m2m_label2{k} = od.l1{j};
                    k = k+1;
                end
            end
            m2m = popTable(nd1, m2m_label1, m2m_label2);
            m2m = sortrows(m2m,'peak_id');
            m2m_csc = reshape(cell2mat(m2m.ccf),[],length(m2m.ccf))';
            imagesc(m2m_csc)
            colormap('jet')
            set(gca, 'XTick', [1,10:10:100])
            set(gca, 'XTickLabel', mat2cell((-0.5:0.1:0.5)',11, 1))
            set(gca, 'YTick', [])
            set(gca, 'TickLength', [0 0])
            title('MSN-MSN')
        end
        
       %Table for FSI-FSI
        if ~isempty(d2)
            subplot(3,1,2)
            ne = (fsi_count*(fsi_count-1))/2;
            f2f_label1 = cell(ne,1);
            f2f_label2 = cell(ne,1);
            nd2 = cellfun(@normdata, d2, 'UniformOutput', false);
            k = 1;
            for i = 1:fsi_count-1
                for j = i+1:fsi_count
                    f2f_label1{k} = od.l2{i};
                    f2f_label2{k} = od.l2{j};
                    k = k+1;
                end
            end
            f2f = popTable(nd2, f2f_label1, f2f_label2);
            f2f = sortrows(f2f,'peak_id');
            f2f_csc = reshape(cell2mat(f2f.ccf),[],length(f2f.ccf))';
            imagesc(f2f_csc)
            colormap('jet')
            set(gca, 'XTick', [1,10:10:200])
            set(gca, 'XTickLabel', mat2cell((-0.5:0.05:0.5)',21, 1))
            set(gca, 'YTick', [])
            set(gca, 'TickLength', [0 0])
            title('FSI-FSI')
        end
        
        %Table for MSN-FSI
        if ~isempty(d3)
            subplot(3,1,3)
            ne = (msn_count*fsi_count);
            m2f_label1 = cell(ne,1);
            m2f_label2 = cell(ne,1);
            nd3 = cellfun(@normdata, d3, 'UniformOutput', false);
            k = 1;
            for i = 1:msn_count
                for j = i:fsi_count
                    m2m_label1{k} = od.l1{i};
                    m2m_label2{k} = od.l2{j};
                    k = k+1;
                end
            end
            m2f = popTable(nd3, m2f_label1, m2f_label2);
            m2f = sortrows(m2f,'peak_id');
            m2f_csc = reshape(cell2mat(m2f.ccf),[],length(m2f.ccf))';
            imagesc(m2f_csc)
            colormap('jet')
            set(gca, 'XTick', [1,10:10:200])
            set(gca, 'XTickLabel', mat2cell((-0.5:0.05:0.5)',21, 1))
            set(gca, 'YTick', [])
            set(gca, 'TickLength', [0 0])
            title('MSN-FSI')
        end        
        WriteFig(fig,strcat(curRat,'-Day-',num2str(jdx)),1);
        close all;
    end
end


% load('/Users/manishm/Work/vanDerMeerLab/RandomVStrDataAnalysis/ccf_results/Smooth/ccf_R131-2007-09-06_od.mat');
% msn_count = length(od.l1);
% fsi_count = length(od.l2);
% %Table for MSN-FSI
% ne = (msn_count*(msn_count-1))/2;
% m2f_label1 = cell(ne,1);
% m2f_label2 = cell(ne,1);
% d1 = od.cx1;
% nd1 = cellfun(@normdata, d1, 'UniformOutput', false);
% k = 1;
% for i = 1:msn_count-1
%     for j = i+1:msn_count
%         m2f_label1{k} = od.l1{i};
%         m2f_label2{k} = od.l1{j};
%         k = k+1;
%     end
% end
% m2f = popTable(nd1, m2f_label1, m2f_label2);
% m2f = sortrows(m2f,'peak_id');

% q0 = reshape(cell2mat(m2f.ccf),[],length(m2f.ccf))';

function res_table = popTable(ndata, label1, label2)
    n = length(label1);
    res_table = table;
    vid = [];
    for i = 1:n
        if ~isnan(ndata{i}(1)) %Omit if ndata is all nan
            vid = [vid i];
        end
    end
    res_table.uid = [1:length(vid)]';
    res_table.label1 = label1(vid);
    res_table.label2 = label2(vid);
    [~ , res_table.peak_id] = cellfun(@max, ndata(vid));
    res_table.ccf = ndata(vid);
end

% function to normalize data
function ndata = normdata(data)
    if (sum(isnan(data)) == length(data))
        ndata = data;
    else
        maxd = max(data);
        mind = min(data);
        ndata = (data - mind)/(maxd-mind);
    end
end

