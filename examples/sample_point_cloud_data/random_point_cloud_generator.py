import numpy as np
import random
import os
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def generate_random_spheres(Ns, Nq, num_points_per_sphere):
    spheres = []
    for _ in range(Ns):
        center = np.random.rand(Nq)  # Random center in Nq-dimensional space
        radius = random.uniform(0.05, 0.1)  # Random radius between 1 and 10
        spheres.append((center, radius))
    return spheres

def sample_points_in_sphere(center, radius, num_points):
    points = []
    for _ in range(num_points):
        direction = np.random.randn(len(center))
        direction /= np.linalg.norm(direction)
        distance = random.uniform(0, radius)
        point = center + direction * distance
        points.append(point)
    return points

def save_point_cloud(points, filename):
    # os.makedirs(os.path.dirname(filename), exist_ok=True)
    with open(filename, 'w') as f:
        for point in points:
            f.write(' '.join(map(str, point)) + '\n')

def visualize_point_cloud(points):
    points = np.array(points)
    if points.shape[1] == 2:
        plt.scatter(points[:, 0], points[:, 1])
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.title('2D Point Cloud')
        plt.show()
    elif points.shape[1] == 3:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(points[:, 0], points[:, 1], points[:, 2])
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_zlabel('Z')
        ax.set_title('3D Point Cloud')
        plt.show()
    else:
        print("Visualization is only supported for 2D or 3D point clouds.")

def main():
    Ns = 10  # Number of spheres
    Nq = 2   # Dimensionality of space
    pc_index = 1  # Index of the Point Cloud
    num_points_per_sphere = 100  # Number of points to sample per sphere

    spheres = generate_random_spheres(Ns, Nq, num_points_per_sphere)
    all_points = []

    for center, radius in spheres:
        points = sample_points_in_sphere(center, radius, num_points_per_sphere)
        all_points.extend(points)

    save_point_cloud(np.array(all_points), f'examples/sample_point_cloud_data/point_cloud_{pc_index}.txt')
    visualize_point_cloud(all_points)

if __name__ == "__main__":
    main()