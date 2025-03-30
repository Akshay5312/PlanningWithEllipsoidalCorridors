experiment_number = 1; %
 
%Read the collision points from sample_point_cloud_data/point_cloud_{experiment_number}.txt
%the point cloud data is stored such that each row is a point in the format [x y z]. If in 2d, it's [x y]
point_cloud = [];
fileID = fopen(['sample_point_cloud_data/point_cloud_' num2str(experiment_number) '.txt'], 'r');
while ~feof(fileID)
    line = fgetl(fileID);
    if ischar(line)
        point = str2num(line); %#ok<ST2NM>
        point_cloud = [point_cloud; point];
    end
end
fclose(fileID);

point_cloud


%$Also read the ellipsoid cloud data from sample_result_data/ellipsoid_cloud_{experiment_number}.txt
%the ellipsoid cloud data is stored such that each row is a point in the format [x y z]. If in 2d, it's [x y]
ellipsoid_cloud = [];
fileID = fopen(['sample_result_data/ellipsoid_' num2str(experiment_number) '.txt'], 'r');
while ~feof(fileID)
    line = fgetl(fileID);
    if ischar(line)
        point = str2num(line); %#ok<ST2NM>
        ellipsoid_cloud = [ellipsoid_cloud; point];
    end
end
fclose(fileID);
%Visualize the ellipsoid cloud

%Visualize the point cloud
figure;



%Visualize the ellipsoid cloud
if size(ellipsoid_cloud, 2) == 2
    plot(ellipsoid_cloud(:,1), ellipsoid_cloud(:,2), 'b.'); % Blue color for 2D ellipsoids
elseif size(ellipsoid_cloud, 2) == 3
    plot3(ellipsoid_cloud(:,1), ellipsoid_cloud(:,2), ellipsoid_cloud(:,3), 'b.'); % Blue color for 3D ellipsoids
else
    error('Ellipsoid cloud data must be 2D or 3D.');
end

hold on;

if size(point_cloud, 2) == 2
    plot(point_cloud(:,1), point_cloud(:,2), 'r.'); % Red color for 2D points
    xlabel('X-axis');
    ylabel('Y-axis');
    title('2D Point Cloud');
elseif size(point_cloud, 2) == 3
    plot3(point_cloud(:,1), point_cloud(:,2), point_cloud(:,3), 'r.'); % Red color for 3D points
    xlabel('X-axis');
    ylabel('Y-axis');
    zlabel('Z-axis');
    title('3D Point Cloud');
else
    error('Point cloud data must be 2D or 3D.');
end

axis equal;
grid on;
%Add a title and labels to the axes
%Save the figure
saveas(gcf, ['sample_result_data/ellipsoid_experiment_' num2str(experiment_number) '.png']);
%Close the figure
close(gcf);