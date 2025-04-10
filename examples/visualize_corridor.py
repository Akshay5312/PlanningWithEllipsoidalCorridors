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

ellipsoid_file = f'examples/sample_result_data/corridor_{experiment_number}.txt'
ellipsoid = read_points(ellipsoid_file)

path_file = f'examples/sample_result_data/path_{experiment_number}.txt'
path = read_points(path_file)

path = np.array(path)

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
    x0 = [0.1, 0.1]
    xf = [0.9, 0.9]
    # create a line between x0 and xf
    line = np.linspace(x0, xf, 100)

    #plot the point cloud
    fig, ax = plt.subplots()
    #plot the ellipsoid
    ellipsoid = np.array(ellipsoid)
    ax.scatter(*zip(*ellipsoid), c='r', marker='o', label='Corridor', s=0.01)
    ax.scatter(*zip(*point_cloud), c='b', marker='x', label='Collision Point Cloud')
    #plot the line
    ax.plot(line[:, 0], line[:, 1], c='g', label='reference', linewidth=2)
    ax.plot(path[:,0], path[:,1], c='k', label='Resulting Path', linewidth=2)
    ax.set_xlabel('X axis')
    ax.set_ylabel('Y axis')
    ax.set_title('Point Cloud and Ellipsoid Visualization')
    ax.legend()
    plt.show()

# Save the figure
fig.savefig(f'examples/sample_result_data/corridor_{experiment_number}.png')