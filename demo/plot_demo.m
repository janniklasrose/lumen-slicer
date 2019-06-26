% this script reads and plots the rather cryptic output file of mesh_slicer

%% process

% read centreline first
fid_line = fopen('centreline.dat', 'r');
nPlanes = fscanf(fid_line, '%d\n', 1); % read number of planes
centreline = fscanf(fid_line, '%g %g %g %g %g %g\n', [6, nPlanes]);
orig = centreline(1:3, :);
fclose(fid_line);

% open file
fid_data = fopen('slices.dat', 'r');
assert(fscanf(fid_data, '%d\n', 1) == nPlanes, 'inconsistent number of points');

% read data
planes = cell(1, nPlanes);
for i = 1:nPlanes

    % read curves (may be more than 1 per slice)
    nCurves = fscanf(fid_data, '%d\n', 1); % read number of curves
    curves = cell(1, nCurves);
    areas = zeros(1, nCurves);
    for j = 1:nCurves
        areas(j) = fscanf(fid_data, '%g\n', 1); % area of this curve
        nPoints = fscanf(fid_data, '%d\n', 1); % read number of points
        curves{j} = fscanf(fid_data, '%g %g %g\n', [3, nPoints]); % read points
    end

    % store
    planes{i} = struct('curves', {curves}, 'areas', areas);

end
fclose(fid_data);

%% areas

% compute
t = linspace(0, 1, nPlanes); % along centreline
areas = cellfun(@(p) sum(p.areas), planes); % sum in case of multiple curves

% plot
figure();
plot(t*100, areas, 'x-', 'DisplayName', 'Area (mesh_slicer)');
xlabel('t [%]', 'Interpreter', 'tex');
ylabel('A [units^2]', 'Interpreter', 'tex');
legend('Interpreter', 'none');

%% geometry

% prepare figure
figure();
hold on;
axis equal;
xlabel('x [units]');
ylabel('y [units]');
zlabel('z [units]');
legend();

colour = t; % usually, use areas, but this is not very interesting here

% plotting
plot3(orig(1,:), orig(2,:), orig(3,:), 'ro-', 'DisplayName', 'centreline');
for i = 1:nPlanes
    plane = planes{i};
    curves = plane.curves;
    for j = 1:numel(curves)
        xyz = curves{j};
        fill3(xyz(1,:), xyz(2,:), xyz(3,:), colour(i), ...
            'DisplayName', sprintf('plane %d',i));
    end
end
view(3);
