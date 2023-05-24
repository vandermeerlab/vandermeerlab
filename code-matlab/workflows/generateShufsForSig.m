% Boilerplate code for selecting 'clean cells'

cd('D:\RandomVstrAnalysis\temp2\'); % Change this to your local machine location for results
odir = ('D:\RandomVstrAnalysis\temp_with_shuf\');
rats = {'R117','R119','R131','R132'};
clean_msn = 0;
clean_fsi = 0;
nshufs = 1000;
for idx = 1:length(rats)
    curRat = rats{idx};
    searchString = strcat(curRat,'*ft_spec.mat');
    ofiles = dir(searchString);
    for jdx = 1:length(ofiles)
        load(ofiles(jdx).name); % Load a particular session
        % do fsi stuff
        fsi_labels  = od.label(od.cell_type == 2);
        fsi_labels = cellfun(@(x) extractBefore(x, '.t'), fsi_labels, 'UniformOutput', false);
        for iC = 1:length(fsi_labels)
            if isfield(od.fsi_res.near_spec{iC}, 'flag_no_control_split') && ~od.fsi_res.near_spec{iC}.flag_no_control_split
                % do fsi_stuff
                clean_fsi = clean_fsi + 1;
                shuf_sts = zeros(nshufs, length(od.pool_sts.freq));
                shuf_ppc = zeros(nshufs, length(od.pool_sts.freq));
                for iShuf = 1:nshufs
                    pool_count = length(od.pool_sts.time{1});
                    keep = randperm(pool_count);
                    keep = keep(1:od.fsi_res.near_spec{iC}.spk_count);
                    this_shuf = od.pool_sts;
                    this_shuf.label{1} = fsi_labels{iC}; % Because the pooled STS might have a differnt spike channel label
                    this_shuf.fourierspctrm{1} = od.pool_sts.fourierspctrm{1}(keep,:,:);
                    this_shuf.time{1} = od.pool_sts.time{1}(keep,:);
                    this_shuf.trial{1} = od.pool_sts.trial{1}(keep,:);
                    shuf_sts(iShuf,:) = nanmean(sq(abs(this_shuf.fourierspctrm{1})));
    
                    cfg_ppc               = [];
                    cfg_ppc.method        = 'ppc0'; % compute the Pairwise Phase Consistency
                    cfg_ppc.spikechannel  = this_shuf.label;
                    cfg_ppc.channel       = this_shuf.lfplabel; % selected LFP channels
                    cfg_ppc.avgoverchan   = 'weighted';
                    this_shuf_ppc = ft_spiketriggeredspectrum_stat(cfg_ppc, this_shuf);
                    shuf_ppc(iShuf,:) = this_shuf_ppc.ppc0;
                end
                od.fsi_res.near_spec{iC}.shuf_sts = shuf_sts;
                od.fsi_res.near_spec{iC}.shuf_ppc = shuf_ppc;
                clear shuf_ppc shuf_sts
            end    
        end
        % do msn stuff
        msn_labels  = od.label(od.cell_type == 1);
        msn_labels = cellfun(@(x) extractBefore(x, '.t'), msn_labels, 'UniformOutput', false);
        for iC = 1:length(msn_labels)
            if isfield(od.msn_res.near_spec{iC}, 'flag_no_control_split') && ~od.msn_res.near_spec{iC}.flag_no_control_split
                % do msn_stuff
                clean_msn = clean_msn + 1;
                shuf_sts = zeros(nshufs, length(od.pool_sts.freq));
                shuf_ppc = zeros(nshufs, length(od.pool_sts.freq));
                for iShuf = 1:nshufs
                    pool_count = length(od.pool_sts.time{1});
                    keep = randperm(pool_count);
                    keep = keep(1:od.msn_res.near_spec{iC}.spk_count);
                    this_shuf = od.pool_sts;
                    this_shuf.label{1} = msn_labels{iC}; % Because the pooled STS might have a differnt spike channel label
                    this_shuf.fourierspctrm{1} = od.pool_sts.fourierspctrm{1}(keep,:,:);
                    this_shuf.time{1} = od.pool_sts.time{1}(keep,:);
                    this_shuf.trial{1} = od.pool_sts.trial{1}(keep,:);
                    shuf_sts(iShuf,:) = nanmean(sq(abs(this_shuf.fourierspctrm{1})));
    
                    cfg_ppc               = [];
                    cfg_ppc.method        = 'ppc0'; % compute the Pairwise Phase Consistency
                    cfg_ppc.spikechannel  = this_shuf.label;
                    cfg_ppc.channel       = this_shuf.lfplabel; % selected LFP channels
                    cfg_ppc.avgoverchan   = 'weighted';
                    this_shuf_ppc = ft_spiketriggeredspectrum_stat(cfg_ppc, this_shuf);
                    shuf_ppc(iShuf,:) = this_shuf_ppc.ppc0;
                end
                od.msn_res.near_spec{iC}.shuf_sts = shuf_sts;
                od.msn_res.near_spec{iC}.shuf_ppc = shuf_ppc;
                clear shuf_ppc shuf_sts
            end
        end
        fn_out = strcat(odir, ofiles(jdx).name);
        save(fn_out,'od');
    end
end
fprintf("Total number of clean MSNs are %d.\n", clean_msn);
fprintf("Total number of clean FSIs are %d.\n", clean_fsi);
