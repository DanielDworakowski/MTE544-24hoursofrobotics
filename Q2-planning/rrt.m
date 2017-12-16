clc; clear; close all;
run('planningEnv.m');

xMax = posMaxBound; % State bounds
xMin = posMinBound;
xR = xMax-xMin;
x0 = [startPos 0]
xF = [endPos ]

% Define single freespace nonconvex polygon by its vertices, each hole
% separated by NaNs
env = [xMin(1) xMin(2);xMin(1) xMax(2);xMax(1) xMax(2);xMax(1) xMin(2); xMin(1) xMin(2)];
obsEdges = [];
figure(1); hold on;
nO = numObsts;
nE = 4;
for i=1:nO
    env = [env; NaN NaN; obsPtsStore(:,2*(i-1)+1:2*i);obsPtsStore(1,2*(i-1)+1:2*i)];
    obsEdges = [obsEdges; obsPtsStore(1:nE,2*(i-1)+1:2*i) obsPtsStore([2:nE 1],2*(i-1)+1:2*i)];
end
%% Vehicle
dt = 0.1;
uMin = [0 -2]*dt;
uMax = [5 2]*dt;
uR = uMax-uMin;
sMin = 30;
sMax = 100;
sR = sMax-sMin;
%% Multi-query PRM, created until solution found
tic;
done = 0;
milestones = [x0 0];
nM = 1;
t= 0;
while ((~done) && (t < 1000))
    t=t+1;
    % Select node to expand
    % Uniform
    % curstone = max(1,min(nM,round(nM*rand(1,1))))
    % Weighted on distance to goal
    for i=1:nM
        d(i) = norm(milestones(i,1:2)-xF);
    end
    [ds,ind] = sort(d);
    w(ind) = exp(-0.1*[1:nM]);
    W = cumsum(w);
    seed = W(end)*rand(1);
    curstone = find(W>seed,1);

    % Get new control input and trajectory
    newstone = 0;
    s = 0;
    while (~newstone && (s < 10))
        s=s+1;
        input = [uR(1)*rand(1,1)+uMin(1) uR(2)*rand(1,1)+uMin(2)];
        steps = sR*rand(1,1)+sMin;
        samples = milestones(curstone,1:3);
        % Dynamics
        for i=2:steps
            samples(i,:) = samples(i-1,:)+[input(1)*cos(samples(i-1,3))*dt input(1)*sin(samples(i-1,3))*dt input(2)*dt]; 
        end
        keep = inpolygon(samples(:,1), samples(:,2), env(:,1),env(:,2));
        
        if (sum(keep)==length(samples(:,1)))
            milestones = [milestones; samples(end,:) curstone];
            newstone = 1;
            nM = nM+1;
            plot(samples(:,1),samples(:,2),'m');
            plot(milestones(end,1),milestones(end,2),'mo');
            drawnow
        end
    end
    % Check if a path from start to end is found
    if (norm(milestones(end,1:2)-xF)<1)
        done = 1;
    end
end

% Find and plot final path through back tracing
done = 0;
cur = nM
curC = milestones(nM,:);
prev = curC(4);
i=2;
p=1;
dtot= 0;
nMiles = 0;
while (~done)
    if (prev == 1)
        done = 1;
    end
    plot([milestones(prev,1) milestones(cur,1)], [milestones(prev,2) milestones(cur,2)],'go','MarkerSize',6, 'LineWidth',2)
    dtot = dtot + norm(milestones(prev,1:2)-milestones(cur,1:2));
    nMiles = nMiles+1;
    cur = prev;
    curC = milestones(cur,:);
    prev = curC(4);
    p=p+1;
end
disp('Time to find a path');
toc;
