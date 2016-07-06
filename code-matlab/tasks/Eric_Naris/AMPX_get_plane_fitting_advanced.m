function stats = AMPX_get_plane_fitting_advanced(cfg_in, power)
%% AMPX_get_plane_fitting_advanced: creates a plane of best fit to the power
%   across the sites on the probe.  The maximum value is then used to
%   generate
%
%          Inputs:
%           - cfg_in: [struct] contains parameters to use such as cfg.x;
%           cfg.y which are the location of the origin in "8x8 probe space"
%           [default is (4.5,4.5) since this is center of the probe though
%           this is not a physical location.
%           - power: [struct] output from AMPX_get_power
%           -
%          Outputs:
%           - mx:   [1xnEvent] x location of the highest magnitude site per event
%           - my:   [1xnEvent] y location of the highest magnitude site per event
%           - mval: [1xnEvent]   value of the highest magnitude per event
%           used as the z value for later plotting in 3D
%
%
% EC - 2016-06-07
%% set up paramters
cfg_def.x = 4.5;  % default for low gamma
cfg_def.y = 4.5;
cfg_def.plot = 1;
cfg = ProcessConfig(cfg_def, cfg_in);


rsq = NaN(1,length(power.power_distrib));
rsq2 = NaN(1,length(power.power_distrib));
mx = NaN(1,length(power.power_distrib));
my = NaN(1,length(power.power_distrib));
mval = NaN(1,length(power.power_distrib));
%%
for iEvt = 1:length(power.power_distrib)
    evt_data = power.power_distrib{iEvt};
    evt_data = evt_data/(min(min(evt_data)));
    evt_data_in = evt_data; % hold the original data in with NaNs
    [X,Y] = meshgrid(1:8,1:8);
    
    % with NaNs
    if cfg.plot ==1
        subplot(221);
        nan_imagesc_ec(evt_data);
        colorbar('eastOutside')
        
        subplot(222);
        plot3(X,Y,evt_data); grid on;
        
        subplot(223);
        plot3(X,Y,evt_data,'.k','MarkerSize',20); grid on;
        
        subplot(224);
        surf(X,Y,evt_data);
    end
    % remove NaNs
    [c,r] = find(isnan(evt_data));
    one_d_evt = find(isnan(evt_data));
    evt_data(find(isnan(evt_data))) = [];
    for ii = 1:length(c)
        X(c(ii),r(ii)) = NaN;
        Y(c(ii),r(ii)) = NaN;
    end
    one_d_x = find(isnan(X));
    X(find(isnan(X))) = [];
    Y(find(isnan(Y))) = [];
    
    
    %%
    XYZ = [X(:) Y(:) evt_data(:)]; % note : unwraps data points
    
    [n,~,p] = affine_fit(XYZ); % n is first eigenvector (orthonormal to plane)
    
    
    %% get the vector from the center point (4.5, 4.5) amnd the max point
    % annon fnc for getting X,Y from the plane data
    plane_x_y = @(X,Y)-(n(1)/n(3)*X+n(2)/n(3)*Y-dot(n,p)/n(3));
    
    % recreate the 8x8
    Z_plane = NaN*meshgrid(1:8,1:8);
    for ii = 1:8
        for jj  = 1:8
            Z_plane(ii, jj) = plane_x_y(ii, jj);
        end
    end
    
    Z_plane_no_nan = Z_plane;
    
    for ii = length(c):-1:1  % remove the NaNs for plane fitting varience metrics
        Z_plane_no_nan(c(ii), r(ii)) = NaN;
    end
    Z_plane_no_nan(find(isnan(Z_plane_no_nan))) = [];
    
    % get the center value
    plane_center = plane_x_y(cfg.x,cfg.y);
    
    % subtract from the plane values to normalize
    plane_norm = Z_plane - plane_center;
    
    % find the max in the plane
    mval(iEvt) = max(max(plane_norm));
    [mx(iEvt), my(iEvt)] = find(plane_norm == mval(iEvt));
    
    
    %% add the plane
    % figure
    % maximize
    % hold on
    [X,Y] = meshgrid(1:8);
    if cfg.plot ==1
        subplot(224)
        hold on
        surf(X,Y,reshape(Z_plane,8,8),'facecolor','blue','facealpha',0.5);
    end
    
    %% trig to get the angle at from the center to the V/L pole
    subplot(223)
    hold on
    plot3([cfg.x 8], [cfg.y 8], [plane_center Z_plane(8,8)], 'r', 'linewidth', 2)
    plot3([cfg.x 8], [cfg.y 8], [plane_center plane_center], 'r', 'linewidth', 2)
    view(90, 0)
    %
    %basic trig
    %     B = 8 - cfg.x;
    %     A = Z_plane(8,8) - plane_center;
    %     theta(iEvt) = atand(A/B);
    %     if Z_plane <= plane_center
    %         theta(iEvt) = -theta(iEvt);
    %     end
    %% location of the max pole of the plane.  used for classification later. 
    if mx(iEvt) ==8 && my(iEvt) == 8
        pole = 88;
    elseif mx(iEvt) ==8 && my(iEvt) == 1
        pole  = 81;
    elseif mx(iEvt) ==1 && my(iEvt) == 8
        pole = 18;
    elseif mx(iEvt) ==1 && my(iEvt) == 1
        pole = 11;
    end
    
    if plane_center <0
        plane_offset = abs(plane_center);
    else
        plane_offset = 0;
    end
    
    if Z_plane(8,8) <= plane_center == 1
        u = [8,8,(abs(Z_plane(8,8))+ plane_center + plane_offset)]; % convert to vectors with origin  = 0,0,0
        v = [8,8, plane_center+plane_offset];
        cos_theta = (dot(u,v)/dot(norm(u), norm(v)));
        theta(1,iEvt) =-1*acosd(cos_theta);
    else
        u = [8, 8, Z_plane(8,8)+plane_offset]; % convert to vectors with origin  = 0,0,0
        v = [8, 8, plane_center+plane_offset];
        cos_theta = (dot(u,v)/dot(norm(u), norm(v)));
        theta(1,iEvt) =acosd(cos_theta);
    end
    if theta(iEvt) > 90;
        disp('this is strange')
    end
    % get the max pole for categorizing them later
    theta(2,iEvt) = pole;
    %% get the orthogonal value to the plane to characterize the direction and magnitude of the plane
    subplot(224)
    %     if cfg.plot ==1
    %         [X,Y] = meshgrid(1:8,1:8);
    %         surf(X,Y,evt_data_in);
    plot3(p(1),p(2),p(3),'bo','markersize',15,'markerfacecolor','blue');
    %         quiver3(p(1),p(2),p(3),n(1)/3,n(2)/3,n(3)/3,'k','linewidth',2)
    %         % plot3([0 p(1)],[0 p(2)],[0 p(3)],'b','markersize',15,'markerfacecolor','blue');
    %         plot3([p(1) p(1) + n(1)/3], [p(2) p(2) + n(2)/3],[p(3) p(3) + n(3)/3], 'm')
    %         p_2 = p; n_2 = n;
    %         plot3([p_2(1) p_2(1)+n_2(1)/1],[p_2(2) p_2(2)+n_2(2)/1], [p_2(3) p_2(3)+n_2(3)/1], 'k')
    %         view(0,0)
    %     end
    point_val(iEvt,:) = p;
    orth_val(iEvt,:) =n*3;
    %% compute R^2 (percent of variance explained)
    % method 1
    all_var = var(evt_data(:));
    red_var = var(evt_data(:)-Z_plane_no_nan(:));
    rsq(iEvt) = ((all_var-red_var)./all_var)*100;
    
    % method 2
    cc = corrcoef(evt_data(:),Z_plane_no_nan(:));
    rsq2(iEvt) = 100*(cc(1,2).^2);
    Z_planes(:,:,iEvt) = Z_plane;
    close all
end
% subtract the center X Y values so that it becomes the origin
mx = mx-cfg.x;
my = my-cfg.y;
stats.rsq = rsq;
stats.rsq2 = rsq2;
stats.mx = mx;
stats.my = my;
stats.point_val = point_val;
stats.orth_val = orth_val;
stats.mval = mval;
stats.z_plane = Z_planes;
stats.theta = theta;
