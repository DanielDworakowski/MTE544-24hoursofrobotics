% Aerial racing

I = imread('Racecourse.png');
map = im2bw(I, 0.4); % Convert to 0-1 image
map = flipud(1-map)'; % Convert to 0 free, 1 occupied and flip.
[M,N]= size(map); % Map size

% Robot start position
dxy = 0.1;
startpos = dxy*[350 250];
checkpoints = dxy*[440 620; 440 665];

% Plotting
figure(1); clf; hold on;
colormap('gray');
imagesc(1-map');
plot(startpos(1)/dxy, startpos(2)/dxy, 'ro', 'MarkerSize',10, 'LineWidth', 3);
plot(checkpoints(:,1)/dxy, checkpoints(:,2)/dxy, 'g-x', 'MarkerSize',10, 'LineWidth', 3 );
xlabel('North (decimeters)')
ylabel('East (decimeters)')
axis equal
