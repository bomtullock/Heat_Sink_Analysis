 %FINITE ELEMENT ANALYSIS OF THERMAL DIFFUSION IN A 3D DOMAIN
%
% This script performs finite element analysis on a 3D domain to solve
% a transient heat conduction problem. The domain consists of multiple regions,
% each with its unique material properties.
%
% The script involves the following main steps:
% 1. Assign material properties based on domain regions.
% 2. Compute local stiffness, mass, and load matrices.
% 3. Assemble global matrices and vectors.
% 4. Perform convective boundary condition adjustments.
% 5. Utilize a transient solver to iterate over time steps.
% 6. Visualize results with a 3D mesh plot.
% 7. Extract temperature along a specific edge and visualize it over time.
% 8. Compare results with external ANSYS data.
% 9. Analyze differences in terms of percentage differences over time.
%
% Dependencies/Related Functions:
% - Requires the MATLAB PDE Toolbox for the `pdeplot3D` function.
% - Requires imported ANSYS data in 'ANSYSdata.csv'.
%
% Assumptions:
% - The Convective_Faces data structure, containing boundary information, is available.
% - Material properties for each region (Fr, Si, Al) are predefined.
% - Initial temperature is assumed to be ambient temperature everywhere.
% - Boundary conditions, such as convection, are applied on specific faces.
%
% Output:
% - Visual plots illustrating the temperature distribution over time.
% - Comparisons between MATLAB's results and ANSYS data.
% - Saves output figures as PNG files.
%
% Notes:
% - Ensure that all necessary data and parameters are correctly initialized before running the script.
% - For larger domains or finer meshes, computation might take longer.
%
% Created by: Tom Bullock and Iestyn Hughes
%
% Tom Bullock Contributions:
% Steady State Solver, Convection Boundary Condition 
% ANSYS Section of report etc 
%
% Iestyn Hughes Contributions: Formatting, Transient Solver, Plots,
% Verfication Techniques, ANSYS Data comparison.


clear
close all

%% ... Material Properties
kappa_Fr = 0.2;  %W/m/K
kappa_Al = 200;  %W/m/K
kappa_Si = 120;  %W/m/K
qv_chip = 5e7;   %W/m^3
T_atm = 20;      %Celc

c_Fr = 1100; % Specific heat of Fr in J/kg/K
c_Al = 903; % Specific heat of Al in J/kg/K
c_Si = 700; % Specific heat of Si in J/kg/K

rho_Fr = 1900; % Density of Fr in kg/m^3
rho_Al = 2710; % Density of Al in kg/m^3
rho_Si = 2330; % Density of Si in kg/m^3


%% ... Boundary Condition Constants
h_convective = 120; %W/m^2/K
h_insulated = 0; %W/m^2/K

%% Time settings
T_end = 150;       % End time in seconds (arbitrary, adjust as needed)
dt = 1;           % Time step in seconds (arbitrary, adjust as needed)
num_timesteps = floor(T_end/dt);

%% Model Import and Mesh Generation
model = createpde;
importGeometry(model, "3D_FIN.stl");
% MESH GENERATION SET HMAX TO 0.5 UNLESS YOU WANT ULTRA LONG RUN TIME EVEN
% 0.5 STILL TAKES ~MINS
mesh = generateMesh(model, "GeometricOrder", "linear", "Hmax", 0.5);

% Plot 3D mesh
pdeplot3D(mesh, NodeLabels="off", ElementLabels="off", FaceAlpha=1)

% Plot geometry with face labels
%subplot(1,3,2)
%gm = fegeometry("3D_Fin_v2.stl");
%pdegplot(gm,FaceLabels="on",EdgeLabels="off",FaceAlpha=0.6,CellLabels="off")

%% Nodes and Edges of Interest

% Convective Faces
face_ids = [5, 6, 3, 1, 4, 7];
for i = 1:length(face_ids)
    eval(sprintf('Nf%d = findNodes(mesh, "region", "Face", %d);', face_ids(i), face_ids(i)));
end

% Edge for analysis
Nedge_of_interest = findNodes(mesh, "region", "Edge", 24);


%% INITIALISE GLOBAL STIFFNESS MATRIX, LOAD VECTOR

tt = mesh.Elements';            % Elements with 4 nodes
pt = mesh.Nodes' / 1000  ;      % Node coordinates (divide by 1000 to convert to m from mm)

nE = size(tt,1);      % Number of elements
nN = size(pt,1);      % Number of nodes
RegionId = (zeros(nE,1));

% Global Stiffness Matrix, Load Vector, and C Matrix
K = sparse(nN,nN);
F = sparse(nN,1);
C = sparse(nN,nN);

%everything but the sides and front and back
Convective_Faces = {Nf5, Nf6, Nf4, Nf7, Nf3, Nf1};
elementsOnFaces = cell(length(Convective_Faces), 1);

T_old = sparse(nN,1);

%% FINDING SUBDOMAIN
% only way to do subdomains with an .stl import is take the avg. of the z coordinate, and then
% use that to determine what material it is. not ideal as some elements
% will be half in half out but with small enough mesh fairly negligble

for ee = 1:nE
    gnodes = tt(ee,:);
    z_avg  = sum(pt(gnodes,3)) / 4;

    %MOBO 
    if (0 <= z_avg) && (1e-3 >= z_avg)
        RegionId(ee,1) = 1;
    %CHIP

    elseif (1e-3 < z_avg) && (4e-3 >= z_avg)
        RegionId(ee,1) = 2;
        
    %FIN
    elseif (4e-3 < z_avg) && (16e-3 >= z_avg)
        RegionId(ee,1) = 3;
    end
end

%% Loop over each element
for ee = 1:nE
    % Global Node IDs for each element
    gnodes = tt(ee,:);

    % GLOBAL X,Y,Z COORDINATES OF EACH NODE IN ELEMENT
    x = pt(gnodes, 1);
    y = pt(gnodes, 2);
    z = pt(gnodes, 3);

    % Points Matrix for volume calculation
    Pm = [1, x(1), y(1), z(1)
          1, x(2), y(2), z(2)
          1, x(3), y(3), z(3)
          1, x(4), y(4), z(4)];

    % Element volume
    V = (1 / 6) * abs(det(Pm));

    %% Gradient Shape Function Matrix
    
    P_inv = inv(Pm)*(6*V);
    
    B = (P_inv([2 3 4],[1 2 3 4])) * (1/(6*V));

    %% MATERIAL PROPERTIES DEPENDING ON DOMAIN
    if RegionId(ee) == 1
        kappa = kappa_Fr;
        qv = 0; %W/m^3
        rho = rho_Fr;
        c = c_Fr;
    
    elseif RegionId(ee) == 2
        kappa = kappa_Si;
        qv = qv_chip; %W/m^3
        rho = rho_Si;
        c = c_Si;
    
    else 
        kappa = kappa_Al;
        qv = 0;       %W/m^3
        rho = rho_Al;
        c = c_Al;
    end


    %% Local Matricies
    lk = kappa * (B' * B) * V;

    % LOCAL LOAD VECTOR 
    lf = (6 * qv * V / 24) * [1;1;1;1];

    % LOCAL MASS MATRIX 
    lc = (rho * c)*(V /20)*(eye(4) + ones(4));

    % Assembly
    K(gnodes,gnodes) = K(gnodes,gnodes) + lk;
    C(gnodes,gnodes) = C(gnodes,gnodes) + lc;
    F(gnodes,1) = F(gnodes,1) + lf; 

    % Loop through all convection faces
    for idx = 1:length(Convective_Faces)

        node_check = intersect(gnodes, Convective_Faces{idx});
        % if true then that element has a face on the convective surface
        if length(node_check) == 3
            % To see what elements are on the face not actually important
            %elementsOnFaces{idx} = [elementsOnFaces{idx}, ee];
            
            % Coordinates of J,K,L nodes for each element of face
            coords = pt(node_check,:);
            xk = coords(1,1); yk = coords(1,2); zk = coords(1,3);
            xl = coords(2,1); yl = coords(2,2); zl = coords(2,3);
            xj = coords(3,1); yj = coords(3,2); zj = coords(3,3);

            % Edge lengths
            alpha = ((xk-xj)^2 + (yk - yj)^2 + (zk-zj)^2)^0.5;
            beta = ((xl-xk)^2 + (yl - yk)^2 + (zl-zk)^2)^0.5;
            phi = ((xj-xl)^2 + (yj - yl)^2 + (zj-zl)^2)^0.5;

            s = 0.5 * (alpha + beta + phi);

            % Area of face
            A_face = sqrt(s * (s-alpha) * (s-beta) * (s-phi)); 

            % Local load and stiffness matrices due to convection
            load = ((h_convective * T_atm * A_face) / 3 ) * [1;1;1];
            stiffness = (h_convective * A_face / 12) * [2,1,1; 1,2,1; 1,1,2];

            % Add the components to global matrix and vector
            F(node_check,1) = F(node_check,1) + load;
            K(node_check,node_check) = K(node_check,node_check) + stiffness;
        end
    end
end

%% Transient Solver

% Intitialise arrays
T_current = T_atm * ones(nN,1); % Assuming initial temperature is zero everywhere
edge_temperatures = zeros(length(Nedge_of_interest), num_timesteps);
max_temperature = zeros(1, num_timesteps);
min_temperature = zeros(1, num_timesteps);

% Initialise contant M matrix
M = C + dt * K;

%% Time loop
hFig = figure('Name','Woah Look its a rotating dinosaur isnt that great'); % Store handle to the figure
hFig.Position = [100, 50, 1200, 700]; % Set position and size: [left, bottom, width, height]
colorRange = [150, 240]; % Set color range to fixed limits
azimuth = 20+90;

for n = 1:num_timesteps
    % Compute the right-hand side
    rhs = C * T_current + dt * F;

    % Solve for T_next using the combined matrix
    T_next = M \ rhs;

    edge_temperatures(:,n) = T_current(Nedge_of_interest);
    max_temperature(:,n) = max(T_current);
    min_temperature(:,n) = min(T_current);

    % Update T_current for the next iteration
    T_current = T_next;

    % Plot results with mesh
    pdeplot3D(model, 'ColorMapData', T_current, 'FaceAlpha', 0.6); % plot the temperature data
    hold on; % hold the current plot
    pdeplot3D(model, 'FaceAlpha', 0, 'EdgeColor', 'none'); % plot the wireframe mesh
    hold off; % release the hold
    title(['Time: ', num2str(n*dt), ' seconds']);
    clim(colorRange); % Fix colorbar limits to 0 to 100 degrees
    cb = colorbar;
    cb.Label.String = 'Temperature (°C)';
    cb.Label.FontWeight = 'bold';
    view(azimuth, 15); % Set view angle (azimuth and elevation)
    %azimuth = azimuth + 5; % Increment angle for next rotation
    drawnow; % This will update the figure window with each plot


end



%% Edge Surface Plot
% Extract coordinates of the edge nodes
edge_coordinates = pt(Nedge_of_interest,:);

% Compute cumulative distances along the edge
cumulative_distances = zeros(size(Nedge_of_interest));
for i = 2:length(Nedge_of_interest)
    dist = norm(edge_coordinates(i,:) - edge_coordinates(i-1,:));
    cumulative_distances(i) = cumulative_distances(i-1) + dist;
    cumulative_distances_mm(i) = cumulative_distances(i) * 1000;
end

% Generate the time values for the x-axis
time_values = (0:dt:T_end-dt);

% Create a new figure
figure;
sgtitle('Temperature Distribution along the Edge over Time', 'FontWeight', 'bold');

% First subplot with the first viewing angle (3D view)
subplot(1, 2, 1);  % 1 row, 2 columns, first plot
surf(time_values, cumulative_distances_mm, edge_temperatures);
xlabel('Time (seconds)', 'FontWeight', 'bold');
ylabel('Distance along the edge (mm)', 'FontWeight', 'bold');
zlabel('Temperature (°C)', 'FontWeight', 'bold');
title('3D View', 'FontWeight', 'bold');
cb = colorbar;
cb.Label.String = 'Temperature (°C)';
cb.Label.FontWeight = 'bold';
view(-45, 30);  % Set viewing angle (azimuth, elevation); adjust these values for best viewing


% Second subplot with the top-down contour view
subplot(1, 2, 2);  % 1 row, 2 columns, second plot
contourf(time_values,cumulative_distances_mm, edge_temperatures);
xlabel('Time (s)', 'FontWeight', 'bold');
ylabel('Distance along edge from the base (mm)' ,'FontWeight', 'bold');
title('Top-Down Contour View', 'FontWeight', 'bold');
colorbar;

% Adjusting the figure properties for better spacing
set(gcf, 'Position', [100, 100, 1500, 600]);  % Adjusts the figure size for better presentation
%print('-dpng', '-r500', 'Final Results.png'); % Save the current figure as 'filename.png' at 300 DPI


%% Verification Plots 

% Import ANSYS data
data = readtable('ANSYS_Data.csv');
time_csv = data.("Time");
max_temp_csv = data.("MaxTemp");
min_temp_csv = data.("MinTemp");


% Linear Interpolation

resolution = 1; %s
common_start_time = max(min(time_csv), min(time_values));
common_end_time = min(max(time_csv), max(time_values));
common_time_range = common_start_time:resolution:common_end_time;  

max_temp_csv_interp = interp1(time_csv, max_temp_csv, common_time_range, 'linear');
min_temp_csv_interp = interp1(time_csv, min_temp_csv, common_time_range, 'linear');

max_temp_matlab_interp = interp1(time_values, max_temperature, common_time_range, 'linear');
min_temp_matlab_interp = interp1(time_values, min_temperature, common_time_range, 'linear');

percentage_diff_max = ((max_temp_csv_interp - max_temp_matlab_interp) ./ max_temp_csv_interp) * 100;
percentage_diff_min = ((min_temp_csv_interp - min_temp_matlab_interp) ./ min_temp_csv_interp) * 100;


% Raw Data Comparison
figure;

subplot(1,2,1); % Interpolated Ansys Data
plot(common_time_range, max_temp_csv_interp, '-.r', 'LineWidth', 2);
hold on;
plot(common_time_range, min_temp_csv_interp, '-.b', 'LineWidth', 2);
title('ANSYS Data');
xlabel('Time (s)');
ylabel('Temperature (°C)');
legend('Max Temp', 'Min Temp');
grid on;
set(gca, 'FontWeight', 'bold');


subplot(1,2,2); % Interpolated MATLAB Data
plot(common_time_range, max_temp_matlab_interp, '-r', 'LineWidth', 2);
hold on;
plot(common_time_range, min_temp_matlab_interp, '-b', 'LineWidth', 2);
title('MATLAB Data');
xlabel('Time (s)');
ylabel('Temperature (°C)');
legend('Max Temp', 'Min Temp');
grid on;
set(gca, 'FontWeight', 'bold');
set(gcf, 'Position', [100, 100, 1200, 600]);  % Adjusts the figure size for better presentation
%print('-dpng', '-r500', 'DataPlots.png'); % Save the current figure as 'filename.png' at 300 DPI

% Percentage Difference
figure;
plot(common_time_range, percentage_diff_max, '-r', 'LineWidth', 2);
hold on;
plot(common_time_range, percentage_diff_min, '-b', 'LineWidth', 2);
title('Percentage Difference Over Time');
xlabel('Time (s)');
ylabel('Percentage Difference (%)');
legend('Max Temp', 'Min Temp');
grid on;
set(gca, 'FontWeight', 'bold');
set(gcf, 'Position', [100, 100, 1200, 600]);  % Adjusts the figure size for better presentation
%print('-dpng', '-r500', 'PercentageDifferences.png'); % Save the current figure as 'filename.png' at 300 DPI


