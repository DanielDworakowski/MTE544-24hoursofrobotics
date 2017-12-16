clc; clear; close all;
run('planningEnv.m');

xMax = posMaxBound; % State bounds
xMin = posMinBound;
xR = xMax-xMin;
x0 = [startPos]
xF = [endPos ]
nO = numObsts;
nE = 4;

% Define single freespace nonconvex polygon by its vertices, each hole
% separated by NaNs
env = [xMin(1) xMin(2);xMin(1) xMax(2);xMax(1) xMax(2);xMax(1) xMin(2); xMin(1) xMin(2)];
obsEdges = [];
figure(1); hold on;
for i=1:nO
    env = [env; NaN NaN; obsPtsStore(:,2*(i-1)+1:2*i);obsPtsStore(1,2*(i-1)+1:2*i)];
    obsEdges = [obsEdges; obsPtsStore(1:nE,2*(i-1)+1:2*i) obsPtsStore([2:nE 1],2*(i-1)+1:2*i)];
end

% Plot obstacles
fig = figure(1); clf; hold on;
plotEnvironment(obsPtsStore,xMin, xMax, x0, xF, fig);
drawnow();
figure(1); hold on;
disp('Time to create environment');

goalIdx = -1;

%% RRT, created until solution found
tic;
done = 0;
milestones = [x0];
nM = 1;
t= 0;
max_iters = 150;
collision_step_size = 0.1;
e = zeros(max_iters,max_iters);
nCon = 1;
while ((~done) && (t < max_iters))
    t=t+1;
    % Select goal location
    % Uniform with growth factor - not in obstacle
    goal_found = false;
    while (~goal_found)
        growth_factor = 2*(t+50)/(max_iters+100);
        cur_goal = (1-growth_factor)*x0 + growth_factor * xR.*rand(1,2);
        if (inpolygon(cur_goal(1), cur_goal(2), env(:,1), env(:,2)))
            goal_found = true;
        end
    end
   % plot(cur_goal(1),cur_goal(2), 'kx', 'MarkerSize', 6)
   
    % Find closest node 
    dist = zeros(1,length(milestones(:,1)));
    for i = 1:length(milestones(:,1))
        dist(i) = norm(cur_goal-milestones(i,1:2));
    end
    [maxdist, curstone] = min(dist);
    
    cur_edge = [milestones(curstone,1:2); cur_goal];
    %plot(cur_edge(:,1), cur_edge(:,2),'k','LineWidth',2);

    steps = floor(norm(cur_edge(1,:)-cur_edge(2,:))/collision_step_size);
    samples = milestones(curstone,1:2);
    for i=2:steps
        samples(i,:) = ((steps-i)/steps)*milestones(curstone,1:2) + (i/steps)*cur_goal;
    end
    %plot(samples(:,1), samples(:,2),'g','LineWidth',2);
    
    keep = inpolygon(samples(:,1), samples(:,2), env(:,1),env(:,2));

    if (sum(keep)==steps)
        milestones = [milestones; samples(end,:)];
        nCon = nCon + 1;
        e(curstone, nCon) = 1;
        e(nCon,curstone) = 1;
        nM = nM+1;
        plot(samples(:,1),samples(:,2),'m');
        plot(milestones(end,1),milestones(end,2),'mo');
    end
    
    % Check if a path from start to end is found
    last_edge = [xF;milestones(end,1:2)];
    steps = floor(norm(last_edge(1,:)-last_edge(2,:))/collision_step_size);
    samples = milestones(end,1:2);
    for i=2:steps
        samples(i,:) = ((steps-i)/steps)*milestones(end,1:2) + (i/steps)*xF;
    end
    keep = inpolygon(samples(:,1), samples(:,2), env(:,1),env(:,2));

    if (sum(keep)==steps)
        milestones = [milestones; samples(end,:)];
        if goalIdx == -1
          nCon = nCon + 1;
          goalIdx = nCon;
        end
        e(nCon - 1, goalIdx) = 1;
        e(goalIdx,nCon - 1) = 1;
        nM = nM+1;
        plot(samples(:,1),samples(:,2),'m');
        plot(milestones(end,1),milestones(end,2),'mo');
%         done = 1;
    end
    if mod(t, 10) == 0
      [sp, sd] = shortestpath_new(milestones, e, 1, goalIdx);
      if ((sd - 61.55) / 61.55) < 0.05
        done = 1;
      end
    end
end
toc;
% % Find and plot final path through back tracing
% if (done)
  [sp, sd] = shortestpath_new(milestones, e, 1, goalIdx);
  sd
  for i=1:length(sp)-1
      plot(milestones(sp(i:i+1),1),milestones(sp(i:i+1),2), 'go-', 'LineWidth',3);
  end
  
%     done = 0;
%     cur = nM
%     curC = milestones(nM,:);
%     prev = curC(3);
%     i=2;
%     p=1;
%     dtot= 0;
%     nMiles = 0;
%     while (~done)
%         if (prev == 1)
%             done = 1;
%         end
%         plot([milestones(prev,1) milestones(cur,1)], [milestones(prev,2) milestones(cur,2)],'go','MarkerSize',6, 'LineWidth',2)
%         dtot = dtot + norm(milestones(prev,1:2)-milestones(cur,1:2));
%         nMiles = nMiles+1;
%         cur = prev;
%         curC = milestones(cur,:);
%         prev = curC(3);
%         p=p+1;
%     end
%     disp('Time to find a path');
%     toc;
% else
%     disp('No path found');
% end
