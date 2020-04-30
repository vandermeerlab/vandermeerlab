%%Script to generate plots%%

% cd('/Users/manishm/Work/vanDerMeerLab/RandomVStrDataAnalysis/temp');
% % rats = {'R117','R119','R131','R132'};
% rats = {'R119','R131','R132'};
% for idx = 1:length(rats)
%     curRat = rats{idx};
%     searchString = strcat('*',curRat,'*od.mat');
%     ofiles = dir(searchString);
%     for jdx = 1:length(ofiles)
%         figure;
%         f_prefix = strcat(curRat,'-Day-',num2str(jdx));
%         f_1 = strcat(f_prefix,'-MSN');
%         f_2 = strcat(f_prefix,'-FSI');
%         f_3 = strcat(f_prefix,'-mixed');
%         load(ofiles(jdx).name);
%         d1 = od.cx1;
%         nc1 = length(od.l1);
%         squareFrame(d1,nc1, f_1);
%         d2 = od.cx2;
%         nc2 = length(od.l2);
%         squareFrame(d2,nc2, f_2);
%         d3 = od.cx3;
%         rectFrame(d3,nc1, nc2, f_3);
%     end
% end


load('/Users/manishm/Work/vanDerMeerLab/Common/temp/ccf_R117-2007-06-01_od.mat');
d1 = od.cx1;
nc1 = length(od.l1);
squareFrame(d1,nc1, 'test1');
d2 = od.cx2;
nc2 = length(od.l2);
squareFrame(d2,nc2, 'test2');
d3 = od.cx3;
rectFrame(d3,nc1, nc2, 'test3');

function sq = squareFrame(data, ncount, o_fname)
    n = ncount-1;
    m = 1;
    fig = figure;
    for i=1:n
            subplot(n,n,((i-1)*n+i:n:(n*n)));
            t = zeros(n+1-i,101);
            for j = m:m+n-i
                t(j+1-m,:) = data{j};
            end
            m = j+1;
            s = stackedplot((-.5:0.01:0.5),t');
            top_label = num2str(i);
            %hack for n=1
            if (n==1)
                side_label = {'2'};
            else
                if (i==1)
                    side_label = cell(n,1);
                    for sC = 1:length(side_label)
                        side_label{sC} = num2str(sC+1);
                    end
                else
                    side_label = cell(n+1-i,1);
                    for sC = 1:length(side_label)
                        side_label{sC} = '';
                    end
                end
            end
            s.set('DisplayLabels',side_label);
            s.set('Title',top_label);
            s.set('Fontsize',11);
    end
    WriteFig(fig,o_fname,1);
    close all;
end

function rt = rectFrame(data, ncount1, ncount2, o_fname)
    n1 = ncount1;
    n2 = ncount2;
    fig = figure;
    for i=1:n2
            subplot(n1,n2,(i:n2:(n1*n2)));
            t = zeros(n1,101);
            for j = i:n2:(n1*n2)
                t(((j-i)/n2)+1,:) = data{j};
            end
            s = stackedplot((-.5:0.01:0.5),t');
            top_label = num2str(i);
            side_label = cell(n1,1);
            if (i==1)
                for sC = 1:length(side_label)
                    side_label{sC} = num2str(sC);
                end
            else
                for sC = 1:length(side_label)
                    side_label{sC} = '';
                end
            end
            s.set('DisplayLabels',side_label);
            s.set('Title',top_label);
            s.set('Fontsize',11);
    end
    WriteFig(fig,o_fname,1);
    close all;
end
        
