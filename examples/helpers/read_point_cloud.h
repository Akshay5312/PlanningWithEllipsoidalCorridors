#pragma once

#include <iostream>

#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <stdexcept>
#include <Eigen/Dense>


Eigen::MatrixXd readPointCloud(const std::string& filename, int& dimensions) {
    std::vector<Eigen::VectorXd> points;
    std::ifstream file(filename);
    std::string line;

    std::cout << "filename: " << filename << std::endl;

    if (!file.is_open()) {
        throw std::runtime_error("Could not open file");
    }

    dimensions = 0;
    bool firstLine = true;

    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::vector<double> values;
        double value;

        while (iss >> value) {
            values.push_back(value);
        }

        if (firstLine) {
            dimensions = values.size();
            firstLine = false;
        }

        if (values.size() != dimensions) {
            throw std::runtime_error("Inconsistent point dimensions");
        }

        Eigen::VectorXd point(dimensions);
        for (int i = 0; i < dimensions; ++i) {
            point[i] = values[i];
        }
        points.push_back(point);
    }

    file.close();

    Eigen::MatrixXd pc_matrix(points.size(), dimensions);
    for (size_t i = 0; i < points.size(); ++i) {
        for (int j = 0; j < dimensions; ++j) {
            pc_matrix(i, j) = points[i][j];
        }
    }

    drake::log()->info("Read {} points with {} dimensions from file: {}", points.size(), dimensions, filename);
    drake::log()->info("Point cloud matrix:\n{}", pc_matrix);

    return pc_matrix.transpose();
}


void savePointCloud(const std::string& filename, const std::vector<Eigen::VectorXd>& point_cloud) {
    std::ofstream file(filename);
    if (file.is_open()) {
        for (const auto& point : point_cloud) {
            file << point.transpose();
            file << "\n";
        }
        file.close();
    } else {
        std::cerr << "Unable to open file: " << filename << std::endl;
    }
}