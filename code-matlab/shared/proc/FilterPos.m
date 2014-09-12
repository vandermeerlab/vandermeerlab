function pos = FilterPos(cfg_in,pos)
% pos = FilterPos(cfg)
%
% Basic Kalman filter for position tracking data
%
% inputs:
%
% pos: position tsd from LoadPos(), DO NOT REMOVE ZEROS (cfg.removeZeros =
% 0)
%
% input cfg fields:
%
% cfg.dt = 0.001; % time step (should be << dt of input for accuracy)
% cfg.sigma_model = 4; % higher: less weight to model
% cfg.P_model = 10^-5; % model noise
% cfg.sigma_meas = 0.01; % measurement noise
% cfg.RB = 1; % run filter forward and backward, take weighted average
%
% output:
%
% pos: filtered position tsd
%
% MvdM 2014-08-09
% based on example by Alex Blekhman

cfg.dt = 0.001;
cfg.sigma_model = 4; % higher: less weight to model
cfg.P_model = 10^-5; % model noise
cfg.sigma_meas = 0.01; % measurement noise
cfg.RB = 1;

ProcessConfig; % note: this converts cfg fields into workspace variables!

mfun = mfilename;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% set up kalman time and observation vectors %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ktvec = pos.tvec(1):cfg.dt:pos.tvec(end);

% find nearest idx (into kalman time vector) for observed samples
near_idx = nearest_idx3(pos.tvec,ktvec);

x_obs = nan(size(ktvec));
x_obs(near_idx) = getd(pos,'x');
x_obs(x_obs == 0) = NaN;

y_obs = nan(size(ktvec));
y_obs(near_idx) = getd(pos,'y');
y_obs(y_obs == 0) = NaN;

%%%%%%%%%%%%%%%%%%%%%%%
%%% set up matrices %%%
%%%%%%%%%%%%%%%%%%%%%%%

% state vector is [position; velocity]
Xk_prev = [nanmean(pos.data(1,1:10)); 5]; % initial estimate
Xk=[]; % current state -- unknown

Yk_prev = [nanmean(pos.data(2,1:10)); 5]; % initial estimate
Yk=[]; % current state -- unknown

% Motion equation: Xk = Phi*Xk_prev + Noise, that is Xk(n) = Xk(n-1) + Vk(n-1) * dt
% Of course, V is not measured, but it is estimated
% Phi represents the dynamics of the system: it is the motion equation
Phi = [1 cfg.dt;
       0  1];

% The error matrix (or the confidence matrix): P states whether we should 
% give more weight to the new measurement or to the model estimate 

% P = sigma^2*G*G';
Q = [(cfg.dt^4)/4  (cfg.dt^3)/3;
                 (cfg.dt^3)/3 cfg.dt^2] * cfg.sigma_model;
             
% G is [(dt.^2)/2; dt] for acceleration

% Q is the process noise covariance. It represents the amount of
% uncertainty in the model. In our case, we arbitrarily assume that the model is perfect (no
% acceleration allowed for the train, or in other words - any acceleration is considered to be a noise)
P = [cfg.P_model 0;
     0 cfg.P_model];

% M is the measurement matrix. 
% We measure X, so M(1) = 1
% We do not measure V, so M(2)= 0
M = [1 0];

% R is the measurement noise covariance. Generally R and sigma_meas can
% vary between samples. 
R = cfg.sigma_meas^2;

%%%%%%%%%%%%%%%%%%
%%% run filter %%%
%%%%%%%%%%%%%%%%%%

fprintf('Forward pass for x: ');
[x_out,S_out1] = runKalman(x_obs,Xk_prev,Xk,Phi,Q,P,M,R);

fprintf('Forward pass for y: ');
[y_out,S_out2] = runKalman(y_obs,Yk_prev,Yk,Phi,Q,P,M,R);

x_out = x_out(1,:); y_out = y_out(1,:);

if cfg.RB % also run filter backwards
    
    Xk_prev = [nanmean(pos.data(1,end-10:end)); 5]; % initial estimate
    Xk=[]; % current state -- unknown
    
    Yk_prev = [nanmean(pos.data(2,end-10:end)); 5]; % initial estimate
    Yk=[]; % current state -- unknown

    fprintf('Backward pass for x: ');
    [x_outR,S_out1R] = runKalman(x_obs(end:-1:1),Xk_prev,Xk,Phi,Q,P,M,R);
    
    fprintf('Backward pass for y: ');
    [y_outR,S_out2R] = runKalman(y_obs(end:-1:1),Yk_prev,Yk,Phi,Q,P,M,R);
    
    x_outR = x_outR(1,:); y_outR = y_outR(1,:);
    
    x_outR = x_outR(end:-1:1); S_out1R = S_out1R(end:-1:1);
    y_outR = y_outR(end:-1:1); S_out2R = S_out2R(end:-1:1);
    
    % combine
    S_out1 = 1./S_out1; S_out2 = 1./S_out2;
    S_out1R = 1./S_out1R; S_out2R = 1./S_out2R;
    x_out = (1./(S_out1+S_out1R)).*((x_out.*S_out1)+(x_outR.*S_out1R));
    y_out = (1./(S_out2+S_out2R)).*((y_out.*S_out2)+(y_outR.*S_out2R));
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% get samples for desired tvec (default: original) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pos.data(1,:) = x_out(near_idx);
pos.data(2,:) = y_out(near_idx);

% housekeeping
pos.cfg.history.mfun = cat(1,pos.cfg.history.mfun,mfun);
pos.cfg.history.cfg = cat(1,pos.cfg.history.cfg,{cfg});

function [out,S_out] = runKalman(in,Xk_prev,Xk,Phi,Q,P,M,R)

nSamples = length(in);

out = zeros(2,nSamples+1);
out(:,1) = Xk_prev;
S_out = zeros(1,nSamples+1);

prev_slen = 0; % counter for fprint output

for k=1:nSamples
    
    % Z is the measurement vector
    Z = in(k);
    
    % Kalman iteration
    P1 = Phi*P*Phi' + Q;
    S = M*P1*M' + R;
    
    % K is Kalman gain. If K is large, more weight goes to the measurement.
    % If K is low, more weight goes to the model prediction.
    K = P1*M'*inv(S);
    
    % if missing data, set K = 0
    if isnan(Z)
        K = zeros(size(K));
        Z = 0; % to prevent NaNs in output
    end
    
    P = P1 - K*M*P1;
    
    Xk = Phi*Xk_prev + K*(Z-M*Phi*Xk_prev);
    out(:,k+1) = Xk;
    
    S_out(k+1) = S;
   
    % For the next iteration
    Xk_prev = Xk;
    
    % Update status
    if mod(k,100) == 0
             
       bc = repmat('\b',[1 prev_slen]);
       fprintf(bc);
        
       pc = (k./nSamples).*100;
       pc = sprintf('%.2f%%',pc);
       fprintf('%s',pc);
       
       prev_slen = length(pc);
    end
end;

fprintf('\n');