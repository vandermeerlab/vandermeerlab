%%
bin_c = -0.95:0.05:0.95;
figure;       

for iT = 1:2
    
    % find all scored events
    
    what = {'score1','score3'}; whatlabel = {'identity shuffle','time shuffle'};
    for iW = 1:length(what)
    
        % get data to process
        temp = eval(cat(2,'out.',what{iW},'(iT)'));
    
        % for observed data, only use scored events (e.g. with at least 5
        % cells active, etc)
        scored_idx = find(~isinf(temp.WIN_rho_obs_pval));
        nScored = length(scored_idx); nEvt = length(temp.WIN_rho_obs_pval);
        
        observed_rhos = temp.WIN_rho_obs(scored_idx);
        
        shuffled_rhos = temp.WIN_rho_shuf(scored_idx,:);
        shuffled_rhos_p = temp.WIN_pval_shuf(scored_idx,:);
        
        shuffled_rhos = shuffled_rhos(~isinf(shuffled_rhos_p));
        
        % histogram
        observed_rhos_h = hist(observed_rhos,bin_c);
        shuffled_rhos_h = hist(shuffled_rhos(:),bin_c);
        
        % plot (v1)
%         subplot(2,2,2*(iT-1)+iW);
%         h(1) = plot(bin_c,observed_rhos_h./nScored,'LineWidth',1);
%         hold on;
%         plot(bin_c,observed_rhos_h./nScored,'b.','MarkerSize',15);
%         
%         h(2) = plot(bin_c,shuffled_rhos_h./length(shuffled_rhos),'r','LineWidth',1);
%         plot(bin_c,shuffled_rhos_h./length(shuffled_rhos),'r.','MarkerSize',15);
%         
%         set(gca,'XLim',[-1 1],'FontSize',18,'YLim',[0 0.1]); box off;
%         legend(h,{'observed',whatlabel{iW}}); legend boxoff; grid on;
%         
        % plot (v2) -- bars seem to be necessary to draw in grouped mode..
%         subplot(2,2,2*(iT-1)+iW);
%         h(1) = bar(bin_c,observed_rhos_h./nScored); 
%         set(h(1),'FaceColor','b','EdgeColor','none','BarWidth',0.4);
%         hold on;
%         h(2) = bar(bin_c,shuffled_rhos_h./length(shuffled_rhos)); 
%         set(h(2),'FaceColor','r','EdgeColor','none','BarWidth',0.4);
%         
%         set(gca,'XLim',[-1 1],'FontSize',18,'YLim',[0 0.1]); box off;
%         legend(h,{'observed',whatlabel{iW}}); legend boxoff; grid on;
%         
        
        subplot(2,2,2*(iT-1)+iW);
        h = bar(bin_c,[observed_rhos_h./nScored;shuffled_rhos_h/length(shuffled_rhos)]'); 
        set(h(1),'FaceColor','b','EdgeColor','none','BarWidth',0.5);
        set(h(2),'FaceColor','r','EdgeColor','none','BarWidth',0.5);
        
        set(gca,'XLim',[-1 1],'FontSize',18,'YLim',[0 0.1]); box off;
        legend(h,{'observed',whatlabel{iW}}); legend boxoff; grid on;
        
        
        
        
        
        
        switch iT
            case 1
                title(sprintf('LEFT nScored %d/%d total',nScored,nEvt));
            case 2
                title(sprintf('RIGHT nScored %d/%d total',nScored,nEvt));
        end
        
    end
    
end