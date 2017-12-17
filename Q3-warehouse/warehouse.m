% Warehouse
warehouse_map = [0 0; 0 50; 100 50; 100 -20; 70 -20; 70 0; 0 0; NaN NaN;
          20 10; 25 10; 25 15; 20 15; 20 10; NaN NaN;
          20 35; 25 35; 25 40; 20 40; 20 35; NaN NaN;
          80 10; 85 10; 85 15; 80 15; 80 10; NaN NaN;
          80 35; 85 35; 85 40; 80 40; 80 35; NaN NaN];
      

% Initial quadrotor position
s = 4;
x0 = [10; 40] + s*randn(2,1);

% Plot initial configuration
figure(1); clf; hold on;
plot(warehouse_map(:,1),warehouse_map(:,2), 'LineWidth',2);
plot(x0(1),x0(2), 'go', 'MarkerSize',8,'LineWidth',2)
axis equal

% Path
path = [5 45; 95 45; 95 -15; 75 -15; 75 25; 10 25; 10 5; 50 5; 50 45; 5 45];
plot(path(:,1), path(:,2), 'g');

% Fiducials
known_fiducials = [26 12.5; 26 37.5; 79 12.5; 79 37.5];
nF = 200;
count = 0;
unknown_fiducials = [];
while (count < nF)
    new_fiducial = [100 70].*rand(1,2) + [0 -20];
    if (inpolygon(new_fiducial(:,1), new_fiducial(:,2), warehouse_map(:,1), warehouse_map(:,2)))
        unknown_fiducials = [unknown_fiducials; new_fiducial];
        count = count + 1;
    end
end
plot(known_fiducials(:,1),known_fiducials(:,2), 'ro', 'MarkerSize',6,'LineWidth',2);
plot(unknown_fiducials(:,1),unknown_fiducials(:,2), 'mx', 'MarkerSize',4,'LineWidth',1);
