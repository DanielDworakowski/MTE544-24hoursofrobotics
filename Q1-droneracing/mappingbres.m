% Occupancy Grid Mapping
clear; clc;

% Simulation time
Tmax = 9999;
dt = 1/5; 

T = 0:Tmax / dt;

% Initial Robot location
% x0 = [5 20 0]
% Robot motions
u = [3 0 -3 0;
     0 3 0 -3];
 ui=1;
% Robot sensor rotation command
% w = 0.3*ones(length(T));

run('racing.m')
% map = map.';

x0 = [startpos ./ dxy -pi];

% Belief map
m = 0.5*ones(M,N);
L0 = log(m./(1-m));
L=L0;

% Sensor model parameters
meas_phi = linspace(deg2rad(-69/2), deg2rad(69/2), 128);
rmax = 10 / dxy; % Max range
rmin = 0.3 / dxy;
alpha = 1; % Width of an obstacle (Distance about measurement to fill in)
beta = 0.05; % Width of a beam (Angle beyond which to exclude) 

%State Initialization
x = zeros(3,length(T)+1);
x(:,1) = x0;

goal0 = [440 642];
goal = goal0;
unseen = ones(size(map));
frontier = zeros(size(map));

spincnt = 0
nspins = 0 
gain = 10
minDist = 10
%% Main simulation
for t=2:length(T)
    % Robot motion
    [x(:,t), spincnt, nspins] =  trajRoll(x(:,t - 1), m, goal, unseen,frontier,spincnt, nspins, gain);
    % Generate a measurement data set
    meas_r = getranges(map, x(:,t),meas_phi,rmax, rmin);
    disttoGoal = norm(goal - x(1:2,t).');
    
%     MODIFIERS TO BIAS WHERE TO GO.
    if disttoGoal < minDist

      goal = startpos;
      startpos = startpos * 10;
      gain = 5;
      if minDist == 10
      minDist = 300;
      else
        minDist = 10
        gain = 0
%         goal = [250, 350]
      end
      
    end
    t
    if t == 30 
    spincnt = 15
    end
%     if mod(t,60) == 0
%       spincnt = 15;
%     end

    %% Map update
    measL = zeros(M,N);
    for i = 1:length(meas_phi)
        % Get inverse measurement model
        invmod = inversescannerbres(M,N,x(1,t),x(2,t),meas_phi(i)+x(3,t),meas_r(i),rmax);
        if (~isempty(invmod))
          
          sz = length(invmod(:,1));
          if meas_r(i) == rmax
            sz = sz - 1;
          end
          for j = 1:sz
              ix = invmod(j,1);
              iy = invmod(j,2);
              il = invmod(j,3);
              % Calculate updated log odds
              L(ix,iy) = L(ix,iy) +log(il./(1-il))-L0(ix,iy);
              measL(ix,iy)= measL(ix,iy) +log(il./(1-il))-L0(ix,iy);
              unseen(ix, iy) = 0;
              frontier(ix, iy) = 0;
          end
          if meas_r(i) == rmax
            sz = sz + 1;
            ix = invmod(sz,1);
            iy = invmod(sz,2);
            il = invmod(sz,3);
            L(ix,iy) = L(ix,iy) +log(il./(1-il))-L0(ix,iy);
            measL(ix,iy)= measL(ix,iy) +log(il./(1-il))-L0(ix,iy);
            frontier(ix, iy) = 1 * unseen(ix,iy);

          end
        end
    end
    % Calculate probabilities
    m = exp(L)./(1+exp(L));
    invmod_T = exp(measL)./(1+exp(measL));

    %% Plot results
    
    % Map and vehicle path
    figure(1);clf; hold on;
%     image(100*(1-map));
image(100*frontier);
    colormap(gray);
    plot(x(2,1:t),x(1,1:t),'bx-')
    axis([0 N 0 M])

%     % Inverse measurement model
%     figure(2);clf; hold on;
%     image(100*(invmod_T));
%     colormap(gray);
%     plot(x(2,t),x(1,t),'bx')
%     for i=1:length(meas_r)
%         plot( x(2,t)+meas_r(i)*sin(meas_phi(i) + x(3,t)),x(1,t)+meas_r(i)*cos(meas_phi(i)+ x(3,t)),'ko')
%     end
%     axis([0 N 0 M])
%     %F2(t-1) = getframe;
%     title('Measurements and inverse measurement model');

    % Belief map
    figure(3);clf; hold on;
    image(100*(m));
    colormap(gray);
    plot(x(2,max(1,t-10):t),x(1,max(1,t-10):t),'bx-')
    axis([0 N 0 M])
    %F3(t-1) = getframe;
    title('Current occupancy grid map')

end

function [x_plus] = motModel(x, v, theta, dt)
  if (abs(v) > 20)
    v = sign(v) * 20;
  end
  v = v / 0.1;
  x_plus = x; 
  x_plus(1) = x_plus(1) + v * cos(x(3)) * dt; 
  x_plus(2) = x_plus(2) + v * sin(x(3)) * dt; 
%   x_plus(3) = x_plus(3) + angleWrap(theta) * dt;
  x_plus(3) = angleWrap(theta);
end

function [x_new, spincnt, nspins] = trajRoll(x_cur, map, xF, unseen, frontier, spincnt, nspins, gain)
  dt = 1/5;
  uMin = [0.8 -1.8 + x_cur(3)]; % bounds on inputs, [velocity, rotation rate]
  uMax = [0.8 1.8 + x_cur(3)]; % bounds on inputs, [velocity, rotation rate]
%   uMin = [1. -0.3];
%   uMax = [1. 0.3];
  uR = uMax-uMin; % range of inputs
  sMin = 10; % steps to move
  sMax = 30; % steps to compute for rollout
  n_traj = 15; % number of trajectories to roll out
  score_step = inf;
  x_new = x_cur;
  checkLength = 5;
  [frontRows, frontCols] = find(frontier);
  [N M] =  size(map);
  for i = 1:n_traj
      % Constant speed, linear distribution of turn rates
      input = [uR(1)/2+uMin(1) uR(2)*(i-1)/(n_traj-1)+uMin(2)];
      steps = sMax; 

      % Propagate Dynamics
      x = x_cur;
      invalid = 0;
      for j=2:steps
%           x(:, j) = x(:,j-1)+[input(1)*cos(x(3,j-1))*dt; input(1)*sin(x(3,j-1))*dt; input(2)*dt]
          x(:, j) = motModel(x(:, j-1), input(1), input(2), dt);
          if (x(2,j) <=1||x(2,j)>=M||x(1,j)<=1||x(1,j)>=N)
            invalid = 1;
          end
      end
      if invalid
        continue
      end
        
      idx = sub2ind(size(map), uint32(x(1,:)), uint32(x(2,:)));
      keep = map(idx) < 0.4;
      if (sum(keep)==steps)

        iX = x(1,end);
        iY = x(2,end);
        dist = 0;
        if ~isempty(frontRows)
%           [frontRows - iY, frontCols - iX]
% norm([frontRows - iY, frontCols - iX])
%           norm([frontRows - iX, frontCols - iY],1);
          dist = min(sqrt(sum([frontRows - iX, frontCols - iY].^2,2)))
          if dist > 200
          dist = 0
          end
        end
         
        
%         
%         unseenScore = unseen(x(1, end), x(2, end))
        plot(x(2,:),x(1,:),'g');
        % Score the trajectory
%         thetaWeight = 

        togo_cur = norm(x(1:2,end).'-xF) + gain * dist;

        score_cur = togo_cur %- 0.1*obs_dist;
        if (score_cur < score_step)
            score_step = score_cur;

            x_new = x(:,sMin);
            x_plot = x;
        end
      else
          plot(x(2,:),x(1,:),'r');
      end
  end
  % Check if no progress is made
  if (all(x_new==x_cur) || ((spincnt > 0) && nspins < 1))
      x_new = x_cur;
      x_new(3)=x_new(3)-0.4;
      spincnt = spincnt - 1
      if spincnt < 0
        spincnt = 15;
      end
      if spincnt == 1
        nspins = nspins + 1;
      end
      display('here')
  else
      plot(x_plot(2,:),x_plot(1,:),'b');
      plot(x_plot(2,end),x_plot(1,end),'bo');
      nspins = 0;
  end
  drawnow;
end

function [wrapped_angle] = angleWrap(angle)

wrapped_angle = mod(angle + pi, 2*pi) - pi;
end