clc
clear
%% first dataset: cube, cylinder, sphere all separated, points only on the surfaces
% Random seed in order to always obtain the same points
rng(23); 

% Total number of points
num_points = 1800;

% Creation of the points on the surface of the sphere
% Ray of the sphere
r_sphere = 1.5;

% Azimuthal angles of the points of the sphere
theta = 2 * pi * rand(num_points / 3, 1); 

% Polar angles of the points of the sphere
phi = acos(2 * rand(num_points / 3, 1) - 1); 

% Coordinates of the points, centre : (0, 0, -5)
x_sphere = r_sphere * sin(phi) .* cos(theta);
y_sphere = r_sphere * sin(phi) .* sin(theta);
z_sphere = r_sphere * cos(phi) - 5;

% Points of the sphere
sphere_points = [x_sphere, y_sphere, z_sphere];

% Creation of the points on the surface of the cube
% Cube edge
cube_edge = 2.5; 
cube_points = zeros(num_points / 3, 3);

for i = 1:num_points / 3
    % Select a face 
    face = randi(6); 

    % Random point on the interval [-cube_edge/2, cube_edge/2]
    new_point = cube_edge * (rand(1, 3) - 0.5); 
    
    switch face
        case 1, new_point(1) = cube_edge / 2; 
        case 2, new_point(1) = -cube_edge / 2; 
        case 3, new_point(2) = cube_edge / 2;  
        case 4, new_point(2) = -cube_edge / 2; 
        case 5, new_point(3) = cube_edge / 2;  
        case 6, new_point(3) = -cube_edge / 2; 
    end

    % Move the points
    cube_points(i, :) = new_point - [2, 2.5, 2.5]; 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Creation of the points on the surface of the cylinder
% Ray of the cylinder
cylinder_radius = cube_edge / 2;

% Height of the cylinder
cylinder_height = cube_edge;
cylinder_points = zeros(num_points / 3, 3);

for i = 1:num_points / 3

    % Azimuthal angles of the points of the cylinder
    angle = 2 * pi * rand;

    % Positions x and y
    x = cylinder_radius * cos(angle);
    y = cylinder_radius * sin(angle);

    % Height of the points (z coordinates)
    z = rand * cylinder_height - cylinder_height / 2;

    % Move the points
    cylinder_points(i, :) = [x - 3, y + 0.25, z - 3]; 
end

% Save the dataset
X = [sphere_points; cube_points; cylinder_points];
save('3d_objects.mat', 'X');

%%
clc 
clear
%% second dataset: cube inside a sphere, points only on the surface
% Fixed random seed
rng(23); 

% Total number of pointsi
num_points = 1000;

% Numeber of points assigned to each object
sphere_ratio = 0.8;  
cube_ratio = 0.2;   
num_sphere_points = round(num_points * sphere_ratio);
num_cube_points = num_points - num_sphere_points;

% Creation of the points on the surface of the sphere
% Ray of the sphere
r_sphere = 3;

% Azimuthal angles of the points of the sphere
theta = 2 * pi * rand(num_sphere_points , 1); 

% Polar angles of the points of the sphere
phi = acos(2 * rand(num_sphere_points , 1) - 1); 

% Coordinates of the points, centre : (0, 0, 0)
x_sphere = r_sphere * sin(phi) .* cos(theta);
y_sphere = r_sphere * sin(phi) .* sin(theta);
z_sphere = r_sphere * cos(phi) ;

% Points of the sphere
sphere_points = [x_sphere, y_sphere, z_sphere];

% Creation of the points on the surface of the cube
% Cube edge
cube_edge = 1; 
cube_points = zeros(num_cube_points, 3);

for i = 1:num_cube_points
    % Select a face 
    face = randi(6); 

    % Random point on the interval [-cube_edge/2, cube_edge/2]
    new_point = cube_edge * (rand(1, 3) - 0.5); 
    
    switch face
        case 1, new_point(1) = cube_edge / 2; 
        case 2, new_point(1) = -cube_edge / 2; 
        case 3, new_point(2) = cube_edge / 2;  
        case 4, new_point(2) = -cube_edge / 2; 
        case 5, new_point(3) = cube_edge / 2;  
        case 6, new_point(3) = -cube_edge / 2; 
    end

    % Move the points
    cube_points(i, :) = new_point; 
end

% Save the dataset
X = [sphere_points; cube_points; cube_points];
save('3d_sphere_cube.mat', 'X');
