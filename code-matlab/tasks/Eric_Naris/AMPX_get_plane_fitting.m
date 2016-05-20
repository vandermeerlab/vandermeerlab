function stats = AMPX_get_plane_fitting(power)
%% AMPX_plane_sandbox: Computes the R^2 of a best fit plane to the 3D power values from AMPX_pipeline_sfn
%
%
%
%%
[X,Y] = meshgrid(1:8,1:8);

%% some different views of a single event
for iEvt = 1:length(power.power_distrib)
    evt_data = power.power_distrib{iEvt};
%     f = figure(10000);
%     maximize
    [c,r] = find(isnan(evt_data));
%     subplot(221)
%     nan_imagesc_ec(evt_data);
%     colorbar('eastOutside')
%     
%     subplot(222);
%     plot3(X,Y,evt_data); grid on;
%     
%     subplot(223);
%     plot3(X,Y,evt_data,'.k','MarkerSize',20); grid on;
%     
%     subplot(224);
%     surf(X,Y,evt_data);
    
    %% fit plane
    one_d_evt = find(isnan(evt_data));
    evt_data(find(isnan(evt_data))) = [];
    for ii = 1:length(c)
        X(c(ii),r(ii)) = NaN;
        Y(c(ii),r(ii)) = NaN;
    end
    one_d_x = find(isnan(X));
    X(find(isnan(X))) = [];
    Y(find(isnan(Y))) = [];
    
    XYZ = [X(:) Y(:) evt_data(:)]; % note : unwraps data points
    
    [n,V,p] = affine_fit(XYZ); % n is first eigenvector (orthonormal to plane)
    
%     subplot(222); hold on;
    Z = - (n(1)/n(3)*X+n(2)/n(3)*Y-dot(n,p)/n(3)); % some 3D trigonometry...
    
    %% put the values back for ploting  
    Z_plot = Z;
    if isempty(one_d_evt) ==0
        for ii  = 1:length(one_d_evt)
            if one_d_evt(ii)> length(Z_plot)
                Z_plot = [Z_plot(1:one_d_evt(ii)-1), NaN];
            elseif one_d_evt(ii) == 1
                Z_plot = [NaN, Z_plot];
            else
                Z_plot = [Z_plot(1:one_d_evt(ii)-1), NaN, Z_plot(one_d_evt(ii):end)];
            end
        end
    else
        Z_plot = Z;
    end
    [X,Y] = meshgrid(1:8,1:8);
%     surf(X,Y,reshape(Z_plot,8,8),'facecolor','blue','facealpha',0.5);
%     subplot(223)
%     hold on
%     surf(X,Y,reshape(Z_plot,8,8),'facecolor','blue','facealpha',0.5);
    %% compute R^2 (percent of variance explained)
%     rsq = NaN(length(data_out.(type).power.power_distrib));
%     rsq2 = NaN(length(data_out.(type).power.power_distrib));
    % method 1
    all_var = var(evt_data(:));
    
    red_var = var(evt_data(:)-Z(:));
    rsq(iEvt) = ((all_var-red_var)./all_var)*100;
    
    % method 2
    cc = corrcoef(evt_data(:),Z(:));
    rsq2(iEvt) = 100*(cc(1,2).^2);
    fprintf(['\n R^2_1 = ' num2str(rsq(iEvt)) '    R^2_1 = ' num2str(rsq2(iEvt)) '\n'])
%     ax = axes('position',[0,0,1,1],'visible','off');
%     tx = text(0.35,0.5,[' R^2 var = ' num2str(round(rsq(iEvt))) '%    R^2 corrcoef = ' num2str(round(rsq2(iEvt))) '%']);
%     set(tx,'fontweight','bold', 'fontsize', 26);
%     tx = text(0.1,0.95,['Event: ' num2str(iEvt)]);
%     set(tx,'fontweight','bold', 'fontsize', 26);
    %% save the figures
%     set(f, 'color', 'w')
%     mkdir(['D:\DATA\' data_out.cfg.fname(1:4)  '\' strrep(data_out.cfg.fname(1:15), '_', '-') '\' date '\' data_out.cfg.type '\plane_stats\' type '\'])
%     print(f,'-dpng','-r300',['D:\DATA\' data_out.cfg.fname(1:4)  '\' strrep(data_out.cfg.fname(1:15), '_', '-') '\' date '\' data_out.cfg.type '\plane_stats\' type '\' data_out.cfg.fname '_' num2str(iEvt) '_plane.png'])
%     saveas(f, ['D:\DATA\' data_out.cfg.fname(1:4)  '\' strrep(data_out.cfg.fname(1:15), '_', '-') '\' date '\' data_out.cfg.type '\plane_stats\' type '\' data_out.cfg.fname '_' num2str(iEvt) '_plane.fig'])
%     close(f)
end
%% save the data_out
stats.rsq = rsq;
stats.rsq2 = rsq2;
% save(['D:\DATA\' data_out.cfg.fname(1:4)  '\' strrep(data_out.cfg.fname(1:15), '_', '-') '\' data_out.cfg.type '_data_out.mat'], 'data_out','-v7.3');

