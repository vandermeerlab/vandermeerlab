function spikes = ss_align( spikes )
% UltraMegaSort2000 by Hill DN, Mehta SB, & Kleinfeld D  - 07/12/2010
%
% ss_align - Multi-channel alignment of waveforms on negative peaks
%
% Usage:
%       spikes = ss_align( spikes )
%
% Description:  
%   This function uses spline interpolation in order to find the true peak
% of the event waveforms in a spikes object.  The threshold crossing used
% in spike detection may not be an optimal location for aligning spikes 
% because (1) finite sampling by a data acquisition card may not capture
% the true moment of threshold crossing and (2) the threshold crossing
% is susceptible to noise. So we use alignment in order to compare spikes
% aligned on their maximum negative deviation from 0.
%
% Alignment is performed on the channel with the largest deviation and then
% this offset is applied to all other channels.
%

    % Check for previous alignment
    if isfield( spikes.info, 'align' )
        disp('Spike data is already aligned.')
    else

        % make helper variables
        w = spikes.waveforms;
        num_spikes     = size(w,1);
        num_samples    = size(w,2);
        num_channels   = size(w,3);
        max_samples    = round(spikes.params.max_jitter * spikes.params.Fs/1000);
        keep_samples = num_samples-max_samples;
        thresh_sample  = spikes.info.detect.align_sample;
        thresh_channel = spikes.info.detect.event_channel;

        % build a matrix of only the waveforms that triggered spike detection
        w2  = zeros( [num_spikes num_samples] );    
        for j = 1:num_spikes
            w2( j, : ) = w(j,:,thresh_channel(j));
        end

        % find alignment point on the triggered channel
        shifts = get_offsets( w2, thresh_sample, max_samples );  

        % apply alignment to all channels
        for j = 1:num_channels
            spikes.waveforms(:,1:keep_samples,j) = respline( w(:,:,j), shifts, max_samples );
        end
        spikes.waveforms(:,keep_samples+1:end,:) = [];

        % update spike times
        spikes.spiketimes = spikes.spiketimes + shifts' / spikes.params.Fs;

        % update alignment flagged
        spikes.info.align.aligned = 1;
    end

    [pca.u,pca.s,pca.v] = svd(detrend(spikes.waveforms(:,:),'constant'), 0);             % SVD the data matrix
    spikes.info.pca = pca;
        
%
%  HELPER FUNCTIONS
%
  
 function shifts = get_offsets( w2, thresh_sample, max_samples )
 % find the alignment point in a matrix of data
 
    num_spikes = size(w2,1);
    num_samples = size(w2,2);
    
    % get spline coefficients for region of interest
    pp = spline( 1:num_samples, w2 );
    coefs = reshape( pp.coefs, [num_spikes,num_samples-1,4] );
    ind  = thresh_sample + -2 + [1:max_samples];
    a = coefs(:,ind,1)';  b = coefs(:,ind,2)';  c = coefs(:,ind,3)'; d =coefs(:,ind,4)';
    
    % find value at positive peak, negative peak, and edge
    p = ( -b + sqrt(b.^2 - 3*a.*c) ) ./ (3*a);
    n = ( -b - sqrt(b.^2 - 3*a.*c) ) ./ (3*a);
    val1 = a.*(p.^3) + b.*(p.^2) + c.*p + d;    
    val2 = a.*(n.^3) + b.*(n.^2) + c.*n + d;
    val3 = w2(:,ind+1)';
    
    % strike out values found out of range or that are complex
    val1( p<0 | p>1 | imag(p)~=0 ) = inf;
    val2( n<0 | n>1 | imag(n)~=0 ) = inf;
     
    % find best in each category
    [val11, pos11] = min(val1,[],1);
    [val22, pos22] = min(val2,[],1);
    [val33, pos33] = min(val3,[],1);
    
    % find best overall
    [val,pos] = min( [val11' val22' val33' ]',[],1 );
    
    % save the peak locations
    peak_loc = zeros([num_spikes 1 ]);
    for j = 1:length(pos)
        if      pos(j)==1, peak_loc(j) = ind(pos11(j)) + p(pos11(j),j);
        elseif  pos(j)==2, peak_loc(j) = ind(pos22(j)) + n(pos22(j),j);
        elseif  pos(j)==3, peak_loc(j) = pos33(j) + thresh_sample - 2;
        end
    end
    shifts = peak_loc - thresh_sample;
       
    
function new_w = respline(w,  shifts, max_s)  
% generate values around new alignment point using spline interpolation

   num_spikes = size(w,1);
   num_samples = size(w,2);
   total_samples = num_samples-max_s;
   
   pp = spline(1:num_samples, w);
   
   % the efficient way to call spline is on a single vector rather than
   % on a stack.  so we are going to concatenate all the waveforms together
   % with zeros in between  
   pp.coefs = reshape(pp.coefs, num_spikes, num_samples-1, []);
   pp.coefs = permute(pp.coefs, [2 1 3]);
   padzeros = zeros(1,num_spikes, 4);
   pp.coefs = cat(1, pp.coefs, padzeros);
   pp.coefs(num_samples,:,4) = w(:,end)';
   pp.coefs = reshape(pp.coefs, [], 4);
   pp.pieces = num_spikes*num_samples;   
   pp.dim = 1;   
   pp.breaks = [1:(pp.pieces+1)];
  
   % get  indices for new waveforms
  
   shift_mat = repmat( shifts,1, total_samples );
   offset    = repmat(([1:num_spikes]-1)' * num_samples, 1,total_samples); 
   ind_mat   = repmat( [1:total_samples], num_spikes, 1 );
   new_inds  = offset + shift_mat + ind_mat;
 
   % evaluate spline at the locations of the new waveforms
   new_w = ppval(pp, new_inds);
