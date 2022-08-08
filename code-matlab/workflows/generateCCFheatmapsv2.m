%% Script to plot ccfs as heatmaps ordered by earliest occurring peak as well maximum correlation value
close all;
cd('/Users/manishm/Work/vanDerMeerLab/RandomVStrDataAnalysis/ccf_results_v2/');
rats = {'R117','R119','R131','R132'};
for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat('*',curRat,'*od.mat');
    ofiles = dir(searchString);
    for jdx = 1:length(ofiles)
        load(ofiles(jdx).name);
        fig = figure;
        d1 = od.S2.cx1;
        d2 = od.S2.cx2;
        d3 = od.S2.cx3;
        msn_count = length(od.S2.l1);
        fsi_count = length(od.S2.l2);
        
        %Table for MSN-MSN
        if ~isempty(d1)
            ne = (msn_count*(msn_count-1))/2;
            m2m_label1 = cell(ne,1);
            m2m_label2 = cell(ne,1);
            nd1 = cellfun(@normdata, d1, 'UniformOutput', false);
            k = 1;
            for i = 1:msn_count-1
                for j = i+1:msn_count
                    m2m_label1{k} = od.S2.l1{i};
                    m2m_label2{k} = od.S2.l1{j};
                    k = k+1;
                end
            end
            m2m = popTable(d1, nd1, m2m_label1, m2m_label2);
            m2m = sortrows(m2m,'peak_id');
            m2m_ep_nccf = reshape(cell2mat(m2m.nccf),[],length(m2m.nccf))';
            m2m_ep_ccf = reshape(cell2mat(m2m.ccf),[],length(m2m.ccf))';
            m2m = sortrows(m2m, 'max_peak', 'descend');
            m2m_mp_nccf = reshape(cell2mat(m2m.nccf),[],length(m2m.nccf))';
            m2m_mp_ccf = reshape(cell2mat(m2m.ccf),[],length(m2m.ccf))';
            nrows = size(m2m_mp_ccf);
            nrows = nrows(1);
            
            subplot(6,2,1)
            imagesc(m2m_ep_nccf)
            colormap('jet')
            set(gca, 'XTick', [1,10:10:100])
            set(gca, 'XTickLabel', mat2cell((-0.5:0.1:0.5)',11, 1))
            set(gca, 'TickLength', [0 0])
            title('MSN-MSN-Earliest-Peak-Norm')
            
            subplot(6,2,2)
            imagesc(m2m_ep_ccf)
            colormap('jet')
            set(gca, 'XTick', [1,10:10:100])
            set(gca, 'XTickLabel', mat2cell((-0.5:0.1:0.5)',11, 1))
            set(gca, 'TickLength', [0 0])
            title('MSN-MSN-Earliest-Peak-Raw')
            
            subplot(6,2,3)
            imagesc(m2m_mp_nccf)
            colormap('jet')
            set(gca, 'XTick', [1,10:10:100])
            set(gca, 'XTickLabel', mat2cell((-0.5:0.1:0.5)',11, 1))           
            set(gca, 'TickLength', [0 0])
            title('MSN-MSN-Max-Peak-Norm')
            
            subplot(6,2,4)
            imagesc(m2m_mp_ccf)
            colormap('jet')
            set(gca, 'XTick', [1,10:10:100])
            set(gca, 'XTickLabel', mat2cell((-0.5:0.1:0.5)',11, 1))           
            set(gca, 'TickLength', [0 0])
            title('MSN-MSN-Max-Peak-Raw')
        end
        
        %Table for FSI-FSI
        if ~isempty(d2)
            ne = (fsi_count*(fsi_count-1))/2;
            f2f_label1 = cell(ne,1);
            f2f_label2 = cell(ne,1);
            nd2 = cellfun(@normdata, d2, 'UniformOutput', false);
            k = 1;
            for i = 1:fsi_count-1
                for j = i+1:fsi_count
                    f2f_label1{k} = od.S2.l2{i};
                    f2f_label2{k} = od.S2.l2{j};
                    k = k+1;
                end
            end
            f2f = popTable(d2, nd2, f2f_label1, f2f_label2);
            f2f = sortrows(f2f,'peak_id');
            f2f_ep_nccf = reshape(cell2mat(f2f.nccf),[],length(f2f.nccf))';
            f2f_ep_ccf = reshape(cell2mat(f2f.ccf),[],length(f2f.ccf))';
            f2f = sortrows(f2f, 'max_peak', 'descend');
            f2f_mp_nccf = reshape(cell2mat(f2f.nccf),[],length(f2f.nccf))';
            f2f_mp_ccf = reshape(cell2mat(f2f.ccf),[],length(f2f.ccf))';
            nrows = size(f2f_mp_ccf);
            nrows = nrows(1);
            
            subplot(6,2,5)
            imagesc(f2f_ep_nccf)
            colormap('jet')
            set(gca, 'XTick', [1,10:10:200])
            set(gca, 'XTickLabel', mat2cell((-0.5:0.05:0.5)',21, 1))           
            set(gca, 'TickLength', [0 0])
            title('FSI-FSI-Earliest-Peak-Norm')
            
            subplot(6,2,6)
            imagesc(f2f_ep_ccf)
            colormap('jet')
            set(gca, 'XTick', [1,10:10:200])
            set(gca, 'XTickLabel', mat2cell((-0.5:0.05:0.5)',21, 1))           
            set(gca, 'TickLength', [0 0])
            title('FSI-FSI-Earliest-Peak-Raw')
            
            subplot(6,2,7)
            imagesc(f2f_mp_nccf)
            colormap('jet')
            set(gca, 'XTick', [1,10:10:200])
            set(gca, 'XTickLabel', mat2cell((-0.5:0.05:0.5)',21, 1))           
            set(gca, 'TickLength', [0 0])
            title('FSI-FSI-Max-Peak-Norm')
            
            subplot(6,2,8)
            imagesc(f2f_mp_ccf)
            colormap('jet')
            set(gca, 'XTick', [1,10:10:200])
            set(gca, 'XTickLabel', mat2cell((-0.5:0.05:0.5)',21, 1))           
            set(gca, 'TickLength', [0 0])
            title('FSI-FSI-Max-Peak-Raw')
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
                    m2f_label1{k} = od.S2.l1{i};
                    m2f_label2{k} = od.S2.l2{j};
                    k = k+1;
                end
            end
            m2f = popTable(d3, nd3, m2f_label1, m2f_label2);
            m2f = sortrows(m2f,'peak_id');
            m2f_ep_nccf = reshape(cell2mat(m2f.nccf),[],length(m2f.nccf))';
            m2f_ep_ccf = reshape(cell2mat(m2f.ccf),[],length(m2f.ccf))';
            m2f = sortrows(m2f, 'max_peak', 'descend');
            m2f_mp_nccf = reshape(cell2mat(m2f.nccf),[],length(m2f.nccf))';
            m2f_mp_ccf = reshape(cell2mat(m2f.ccf),[],length(m2f.ccf))';
            nrows = size(m2f_mp_ccf);
            nrows = nrows(1);
            
            
            subplot(6,2,9)
            imagesc(m2f_ep_nccf)
            colormap('jet')
            set(gca, 'XTick', [1,10:10:200])
            set(gca, 'XTickLabel', mat2cell((-0.5:0.05:0.5)',21, 1))            
            set(gca, 'TickLength', [0 0])
            title('MSN-FSI-Earliest-Peak-Norm')
            
            subplot(6,2,10)
            imagesc(m2f_ep_ccf)
            colormap('jet')
            set(gca, 'XTick', [1,10:10:200])
            set(gca, 'XTickLabel', mat2cell((-0.5:0.05:0.5)',21, 1))            
            set(gca, 'TickLength', [0 0])
            title('MSN-FSI-Earliest-Peak-Raw')
            
            subplot(6,2,11)
            imagesc(m2f_mp_nccf)
            colormap('jet')
            set(gca, 'XTick', [1,10:10:200])
            set(gca, 'XTickLabel', mat2cell((-0.5:0.05:0.5)',21, 1))            
            set(gca, 'TickLength', [0 0])
            title('MSN-FSI-Max-Peak-Norm')
            
            subplot(6,2,12)
            imagesc(m2f_mp_ccf)
            colormap('jet')
            set(gca, 'XTick', [1,10:10:200])
            set(gca, 'XTickLabel', mat2cell((-0.5:0.05:0.5)',21, 1))           
            set(gca, 'TickLength', [0 0])
            title('MSN-FSI-Max-Peak-Raw')
     
        end        
        WriteFig(fig,strcat(curRat,'-Day-',num2str(jdx), '-Away-Trials'),1);
        close all;
    end
end


% load('/Users/manishm/Work/vanDerMeerLab/RandomVStrDataAnalysis/ccf_results/Smooth/ccf_R131-2007-09-06_od.S2.mat');
% msn_count = length(od.S2.l1);
% fsi_count = length(od.S2.l2);
% %Table for MSN-FSI
% ne = (msn_count*(msn_count-1))/2;
% m2f_label1 = cell(ne,1);
% m2f_label2 = cell(ne,1);
% d1 = od.S2.cx1;
% nd1 = cellfun(@normdata, d1, 'UniformOutput', false);
% k = 1;
% for i = 1:msn_count-1
%     for j = i+1:msn_count
%         m2f_label1{k} = od.S2.l1{i};
%         m2f_label2{k} = od.S2.l1{j};
%         k = k+1;
%     end
% end
% m2f = popTable(nd1, m2f_label1, m2f_label2);
% m2f = sortrows(m2f,'peak_id');

% q0 = reshape(cell2mat(m2f.ccf),[],length(m2f.ccf))';

function res_table = popTable(data, ndata, label1, label2)
    n = length(label1);
    res_table = table;
    vid = [];
    for i = 1:n
        if ~isnan(ndata{i}(1)) % Omit if ndata is all nan
            vid = [vid i];
        end
    end
    res_table.uid = [1:length(vid)]';
    res_table.label1 = label1(vid);
    res_table.label2 = label2(vid);
    [~ , res_table.peak_id] = cellfun(@max, ndata(vid));
    [res_table.max_peak, ~] = cellfun(@max, data(vid));
    res_table.nccf = ndata(vid);
    res_table.ccf = data(vid);
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

