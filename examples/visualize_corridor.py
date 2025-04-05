import numpy as np
import matplotlib.pyplot as plt
def read_points(file_path):
    points = []
    with open(file_path, 'r') as file:
        for line in file:
            values = list(map(float, line.split()))
            points.append(np.array(values))
    return points


experiment_number = 1  # Replace with your experiment number

point_cloud_file = f'examples/sample_point_cloud_data/point_cloud_{experiment_number}.txt'
point_cloud = read_points(point_cloud_file)

ellipsoid_file = f'examples/sample_result_data/ellipsoid_corridor_{experiment_number}.txt'
ellipsoid = read_points(ellipsoid_file)

N_q = len(ellipsoid[0])

if(N_q == 3):
    #plot the point cloud
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    #plot the ellipsoid
    ellipsoid = np.array(ellipsoid)
    ax.scatter(*zip(*ellipsoid), c='r', marker='-', label='Ellipsoid')
    ax.scatter(*zip(*point_cloud), c='b', marker='x', label='Point Cloud')
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_zlabel('Z axis')
    ax.set_title('Point Cloud and Ellipsoid Visualization')
    ax.legend()
    plt.show()
else:
    #plot the point cloud
    fig, ax = plt.subplots()
    #plot the ellipsoid
    ellipsoid = np.array(ellipsoid)
    ax.scatter(*zip(*ellipsoid), c='r', marker='o', label='Ellipsoid')
    ax.scatter(*zip(*point_cloud), c='b', marker='x', label='Point Cloud')
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_title('Point Cloud and Ellipsoid Visualization')
    ax.legend()
    plt.show()