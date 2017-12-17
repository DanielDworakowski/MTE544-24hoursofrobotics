clc; clear; close all;
run('warehouse.m');

newfeature = ones(204,1);
mufeat_known = zeros(2*4,1);
S0feat_known = zeros(2*4);

for i=1:4
  mufeat_known((2*(i-1)+1):(2*i)) = known_fiducials(i,:).'; 
end

mufeat_unknown = zeros(2*200,1);
S0feat_unknown = 100*eye(2*200);


% 
% Motion disturbance.
n = 4;
R = 0.02*ones(n,n);
% R(5,5) = 0.002;
R(4,4) = 0.002;
Qi = eye(2) * 0.1;

% Fixed vehicle parameters
v = 5;

% Desired trajectory as a plot of points. The plot of points defines
% multiple line segments all joined together. The vehicle will first follow
% along the first line segment. After it is has driven past the length of
% said line segment, it defines a new line segment to follow by taking the
% next point.

traj_points = path;
traj_point_counter = 1; % Keep track of where we are in the trajectory
           
% Initial conditions in [x y heading]
x0_r = [x0; 0]; 

% Simulation time
Tmax = 80;  % End point
dt =0.05; % Time step
T = 0:dt:Tmax; % Time vector

nO = 4;
nE = 4;
obsEdges = [];
for i=9:(length(warehouse_map) -1)
    if (any(isnan(warehouse_map(i, :))) || any(isnan(warehouse_map(i+1, :))))
      continue
    end
    obsEdges = [obsEdges; warehouse_map(i, :) warehouse_map(i+1, :)];
end
%%
% Simulation setup
% xd = zeros(length(T)-1,3); % Derivative of state ([edot psidot])
x = zeros(length(T),n);  % State ([e psi] 
x(1,:) = [x0_r; 5]; % Initial condition
s_rr = eye(n) * 4;
s_rr(4,4) = 0;
% s_rr(5,5) = 0;
y = zeros(408,1);

% Construct overall mean and variances. 
mu = [x(1,:).'; mufeat_known; mufeat_unknown];
S = [s_rr zeros(n,2*4); zeros(2*4,n) S0feat_known];
S = [S zeros(n+8, 400); zeros(400, n+8) S0feat_unknown];
%%
% figure(1);clf; hold on;
hold on;
% axis([-10 10 -10 10]);
for t=1:length(T)-1
    end_point = traj_points(traj_point_counter+1, :);
    start_point = traj_points(traj_point_counter, :);
    
    [x(t+1,:), next_point, a] = motModel(x(t,:), start_point, end_point, v, dt,n);
 
    
    [y_known, y_unknown, flistknown, flistunknown] = measModel(x(t+1,:).', known_fiducials, unknown_fiducials, obsEdges);
    flist = [flistknown; flistunknown];
    %
    % Fill in measurements. 
    for i=1:length(flistknown)
      if (flistknown(i))
        y(2*(i-1)+1:2*i) = y_known(i,:).'; 
      end
    end
    for i=1:length(flistunknown)
      idx = i + 4;
      if (flistunknown(i))
        y(2*(idx-1)+1:2*idx) = y_unknown(i,:).';

      end
    end
%%
% motion cov update.
% mu(1:5) = x(t+1,:).';
[mu(1:n), next_point, a] = motModel(mu(1:n).', start_point, end_point, v, dt,n);
S(1:n,1:n) = a * S(1:n, 1:n) * a' + R;
% Measurement update. 
    for i=1:204
        if (flist(i))
          if (i < 5)
          display('found a static point');
          end
            % Feature initialization
            if (newfeature(i) == 1)
%                 display('new feature');
                dx = y(2*(i-1)+1);
                dy = y(2*i);
                Rot = rot2D(-mu(3));
                newy = Rot * [dx dy].';
                mu(n+2*(i-1)+1) = mu(1)+newy(1);
                mu(n+2*i) = mu(2)+newy(2);
                newfeature(i) = 0;
            end
            % Linearization
            % Predicted range
%             display('0000000')
%             x(t+1,3)
%             mu(3)
%             x(t+1,1)
%             mu(1)
%             x(t+1,2)
%             mu(2)
%             
%             mu(5+2*(i-1)+1)
%             mu(5+2*i)
            
              
%             S(5+2*i,5+2*i)
            
            
%             plot(mu(5+2*(i-1)+1), mu(5+2*i), 'gx')
       
            mu(1:n) - x(t+1,:).'
            dx = mu(n+2*(i-1)+1) - mu(1);
            dy = mu(n+2*i) - mu(2);
            Fi = zeros(n,n+408);
            Fi(1:n,1:n) = eye(n);
            Fi(n+1:n+2,n+2*(i-1)+1:n+2*i) = eye(2);
            th = mu(3);
            rot = rot2D(th);
            xf = mu(n+2*(i-1)+1);
            xr = mu(1);
            yf = mu(n+2*i);
            yr = mu(2);
            Ht = [-cos(th), sin(th), (-xf*sin(th)+xr*sin(th)-yf*cos(th)+yr*cos(th)), 0,   cos(th), -sin(th),;
                  -sin(th), -cos(th), (xf*cos(th)-xr*cos(th)-yf*sin(th)+yr*sin(th)), 0,   sin(th), cos(th)]*Fi;
%             Ht
            I = y(2*(i-1)+1:2*i)-rot*[dx dy].';
%             I
%             mu(3)
            x(t+1,3)
            % Measurement update
            K = S*Ht'/(Ht*S*Ht'+Qi);
%             display('before')
%             mu(3)
%             mu(5+2*(i-1)+1)
%             mu(5+2*i)
            mu = mu + K*I;
            
%             display('after');
%             mu(5+2*(i-1)+1)
%             mu(5+2*i)
            mu(3) = angleWrap(mu(3));
%             mu(3)
            S = (eye(n+408)-K*Ht)*S;
                        % In cases if S bemoes not positive definite, manually make it
            % P.D.
            if min(eig(S))<0
               S=S-eye(length(S)).*min(eig(S));
               warning('S was manually made positive definite')
            end
        end
    end

    
%     figure(1);
    %subplot(1,2,1); hold on;
    %plot(trees(:,1),trees(:,2),'go', 'MarkerSize',10,'LineWidth',2);
%     plot(xr(1,1:t),xr(3,1:t), 'ro--')
%     plot([xr(1,t) xr(1,t)+1*cos(yaw)],[xr(3,t) xr(3,t)+1*sin(yaw)], 'r-')
    plot(mu(1),mu(2), 'rx')
    %plot([mu_S(1,t) mu_S(1,t)+1*cos(yaw)],[mu_S(3,t) mu_S(3,t)+1*sin(yaw)], 'b-')
    mu_pos = [mu(1) mu(2)];
    S_pos = [S(1,1) S(1,2); S(2,1) S(2,2)];
%     error_ellipse(S_pos,mu_pos,0.75);
%     error_ellipse(S_pos,mu_pos,0.95);
%     plot( [mu(1) mu(1)+2*cos(phides)],[mu(3) mu(3)+2*sin(phides)], 'b')
%     for i=1:204
%         if (flist(i))
%             fi = 2*(i-1)+1;
%             fj = 2*i;
% %             testY = y_unknown + x0.'
% %             plot([mu(1) mu(1)+rxy*cos(y(fj,t)+yaw)], [mu(3) mu(3)+rxy*sin(y(fj,t)+yaw)], 'c');
%             plot(mu(5+fi),mu(5+fj), 'gx')
%             mu_pos = [mu(5+fi) mu(5+fj)];
%             S_pos = [S(5+fi,5+fi) S(5+fi,5+fj); S(5+fj,5+fi) S(5+fj,5+fj)];
%             error_ellipse(S_pos,mu_pos,0.95);
%         end
%     end
%     axis equal
%     axis([-4 6 -1 7])
    title('SLAM with Range & Bearing Measurements')

%%
    % Check if we have travelled the distance of the line segment. 
    % If we have, then get the next point
    if (next_point == 1)
        traj_point_counter = traj_point_counter+1;
        if (traj_point_counter == length(traj_points(:,1)))
            break;
        end
    end
    
    plot(x(1:t,1),x(1:t,2),'bo');

  drawnow
  clear y_known y_unknown flistknown flistunknown a;
end

%% Plotting

% Plotting the trajectory of the vehicle
figure(1);clf; hold on;
plot(x(1:t,1),x(1:t,2),'b-');

for t=1:30:t
      drawbox(x(t,1),x(t,2),x(t,3),.5,1);
end
xlabel('x (m)')
ylabel('y (m)')
axis equal


% Phase portrait
% [crosstrack,heading] = meshgrid(-10:.5:10,-3:.2:3); % Create a grid over values of the crosstrack error and heading error
% delta = max(-delta_max,min(delta_max,heading+atan2(k*crosstrack,velocity)));  % Calculate steering angle at each point
% ed = velocity*sin(heading-delta); % Find crosstrack derivative
% psid = -(velocity*sin(delta))/(robot_length); % Find heading derivative

% psibplus = -atan2(k*crosstrack(1,:),velocity)+delta_max; % Find border of max region
% psibminus = -atan2(k*crosstrack(1,:),velocity)-delta_max; % Find border of min region
% 
% figure(2);clf; hold on;
% quiver(crosstrack,heading, ed, psid)
% plot(crosstrack(1,:),psibplus,'r', 'LineWidth',2);
% plot(crosstrack(1,:),psibminus,'r', 'LineWidth',2);
% axis([-10 10 -3 3])
% xlabel('e (m)')
% ylabel('\psi (rad)')


%%

% Collision test
obsEdges(2,1:2)
[col, edge] = checkCollision(x0.', x0+ 1, obsEdges)
%%

% Measurement test.
[y_known, y_unknown, flistknown, flistunknown] = measModel([x0; pi/2], known_fiducials, unknown_fiducials, obsEdges)
%
hold on
R= rot2D(-pi/2);
testY = y_unknown;
for i=1:length(flistunknown)
  if (flistunknown(i))
    testing = R*testY(i,:).' + x0
    plot(testing(1), testing(2), 'b*')
  end
end
%%


function [x_plus, next_point, a] = motModel(x, start_point, end_point, v,dt,n)
    % The desired trajectory is a line segment consisting of 2 points from
    % the desired trajectory
    delta_max = 25*pi/180; % max steering angle
    k = 2.5; % Gain
    kp = 1;
    ki = 0.01;
    robot_length = 1; % Car length

    v = x(4);
    
    traj_angle = atan2(end_point(2) - start_point(2), end_point(1) - start_point(1));
    
    [crosstrack_error, next_point] = distanceToLineSegment(start_point,end_point,x(1:2));
    
    % Calculate steering angle
    delta = max(-delta_max,min(delta_max, angleWrap(traj_angle - x(3))+ atan2(-k*crosstrack_error,v)));
    % State derivatives
%     xd(1) = v*cos(x(3));
%     xd(2) = v*sin(x(3));
%     xd(3) = v*tan(delta/robot_length);
%     
    % State update
%     x_plus(1) = x(1)+dt*xd(1);
%     x_plus(2) = x(2)+dt*xd(2);
%     x_plus(3) = x(3)+dt*xd(3);
%     x_plus(4) = x(4);
%     x_plus(5) = x(4);

    a = eye(n);
    a(1,3) = -dt * x(4) * sin(x(3));
    a(1,4) = dt * cos(x(3));
    a(2,3) = dt * x(4) * cos(x(3));
    a(2,4) = dt * sin(x(3));
    a(3,4) = dt * tan(delta/robot_length);
    a(4,4) = 1;
%     a(5,5) = 1;
%     a(5,4) = 1;
%     a(5,5) = 1;
    x_plus = a * x.';
    x_plus(4) = 5 + kp*(5-x(4)); % + ki * x(5);
%     x_plus(5) = x_plus(5) + ki * (5 - x(4));

    % angle wrap the heading + noise.
    x_plus(1) = x_plus(1) + 0.02 * randn();
    x_plus(2) = x_plus(2) + 0.02 * randn();
    x_plus(3) = x_plus(3) + 0.002 * randn();
    x_plus(4) = x_plus(4) + 0.002 * randn();
%     x_plus(5) = x_plus(5) + 0.002 * randn();
    x_plus(3) = angleWrap(x_plus(3));
end

function [ inCollision, edge ] = checkCollision( ptA, ptB, obstEdges )
  %CHECKCOLLISION Checks if a line interesects with a set of edges
  %   Detailed explanation goes here

  % Check for each edge
  edge = [];
  for k = 1:size(obstEdges,1)
      % If both vertices aren't to the left, right, above, or below the
      % edge's vertices in question, then test for collision
      if ~((max([ptA(1),ptB(1)])<min([obstEdges(k,1),obstEdges(k,3)])) || ...
           (min([ptA(1),ptB(1)])>max([obstEdges(k,1),obstEdges(k,3)])) || ...
           (max([ptA(2),ptB(2)])<min([obstEdges(k,2),obstEdges(k,4)])) || ...
           (min([ptA(2),ptB(2)])>max([obstEdges(k,2),obstEdges(k,4)])))
          if (EdgeCollision([ptA, ptB], obstEdges(k,:)))
              % Eliminate end-to-end contacts from collisions list
              if (sum(abs(ptA-obstEdges(k,1:2)))>0 && ...
                  sum(abs(ptB-obstEdges(k,1:2)))>0 && ...
                  sum(abs(ptA-obstEdges(k,3:4)))>0 && ...
                  sum(abs(ptB-obstEdges(k,3:4)))>0)

                  edge = k;
                  inCollision = 1 ; % In Collision
                  return
              end
          end
      end
  end
  inCollision = 0 ; % Not in collision
end


function [y_known, y_unknown, flist_known , flist_unknown] = measModel(x, markers_known, markers_unknown, edges)
  szKnown = length(markers_known);
  szunKnown = length(markers_unknown);
  flist_known = zeros(szKnown,1);
  flist_unknown = zeros(szunKnown,1);
  y_known = zeros(szKnown,2);
  y_unknown = zeros(szunKnown,2);
  diff = markers_known - x(1:2).';
  R = rot2D(x(3));
  S = eye(2) * 0.1;
  [QiE, Qie] = eig(S);

  for i = 1:szKnown
    [collide, edge] = checkCollision(x(1:2).', markers_known(i,:), edges);
    if (collide)
      continue;
    end
    rotated =  R * diff(i,:).';
    if (abs(rotated(1)) > 10 || abs(rotated(2)) > 7.5)
      continue;
    end
    flist_known(i) = 1;
    d = QiE*sqrt(Qie)*randn(2,1);
    y_known(i, :) = rotated.' + d.';
    
  end
  diff = markers_unknown - x(1:2).';
  for i = 1:szunKnown
    [collide, edge] = checkCollision(x(1:2).', markers_unknown(i,:), edges);
    if (collide)
      continue;
    end
    rotated = R * diff(i,:).';
    if (abs(rotated(1)) > 10 || abs(rotated(2)) > 7.5)
      continue;
    end
    flist_unknown(i) = 1;
    d = QiE*sqrt(Qie)*randn(2,1);
    y_unknown(i,:) = rotated.' + d.';
  end

end




















