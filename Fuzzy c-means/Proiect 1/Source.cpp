#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <ctime>
#include <fstream>


double euclideanDistance(const std::vector<double>& point1, const std::vector<double>& point2) {
    double distance = 0.0;
    for (size_t i = 0; i < point1.size(); ++i) {
        distance += std::pow(point1[i] - point2[i], 2);
    }
    return std::sqrt(distance);
}


void updateMembershipMatrix(const std::vector<std::vector<double>>& data, const std::vector<std::vector<double>>& centroids, std::vector<std::vector<double>>& membershipMatrix, double fuzziness) {
    size_t clusterCount = centroids.size();
    size_t pointCount = data.size();

    for (size_t i = 0; i < pointCount; ++i) {
        for (size_t j = 0; j < clusterCount; ++j) {
            double distanceToCluster = euclideanDistance(data[i], centroids[j]);
            double distanceSum = 0.0;

            for (size_t k = 0; k < clusterCount; ++k) {
                distanceSum += std::pow(distanceToCluster / euclideanDistance(data[i], centroids[k]), 2.0 / (fuzziness - 1.0));
            }

            membershipMatrix[i][j] = 1.0 / distanceSum;
        }
    }
}


void updateCentroids(const std::vector<std::vector<double>>& data, const std::vector<std::vector<double>>& membershipMatrix, std::vector<std::vector<double>>& centroids) {
    size_t clusterCount = centroids.size();
    size_t pointCount = data.size();
    size_t dimensionCount = data[0].size();

    for (size_t j = 0; j < clusterCount; ++j) {
        for (size_t d = 0; d < dimensionCount; ++d) {
            double numerator = 0.0;
            double denominator = 0.0;

            for (size_t i = 0; i < pointCount; ++i) {
                numerator += std::pow(membershipMatrix[i][j], 2.0) * data[i][d];
                denominator += std::pow(membershipMatrix[i][j], 2.0);
            }

            centroids[j][d] = numerator / denominator;
        }
    }
}


bool hasConverged(const std::vector<std::vector<double>>& oldMembershipMatrix, const std::vector<std::vector<double>>& newMembershipMatrix, double epsilon) {
    size_t pointCount = oldMembershipMatrix.size();
    size_t clusterCount = oldMembershipMatrix[0].size();

    for (size_t i = 0; i < pointCount; ++i) {
        for (size_t j = 0; j < clusterCount; ++j) {
            if (std::abs(oldMembershipMatrix[i][j] - newMembershipMatrix[i][j]) > epsilon) {
                return false;
            }
        }
    }

    return true;
}

void fuzzyCMeans(const std::vector<std::vector<double>>& data, size_t clusterCount, double fuzziness, double epsilon, std::vector<std::vector<double>>& centroids, std::vector<std::vector<double>>& membershipMatrix) {
    size_t pointCount = data.size();
    size_t dimensionCount = data[0].size();

    
    centroids.resize(clusterCount, std::vector<double>(dimensionCount));
    for (size_t i = 0; i < clusterCount; ++i) {
        for (size_t d = 0; d < dimensionCount; ++d) {
            centroids[i][d] = static_cast<double>(std::rand()) / RAND_MAX;
        }
    }

    
    membershipMatrix.resize(pointCount, std::vector<double>(clusterCount));
    for (size_t i = 0; i < pointCount; ++i) {
        double membershipSum = 0.0;
        for (size_t j = 0; j < clusterCount; ++j) {
            membershipMatrix[i][j] = static_cast<double>(std::rand()) / RAND_MAX;
            membershipSum += membershipMatrix[i][j];
        }

        
        for (size_t j = 0; j < clusterCount; ++j) {
            membershipMatrix[i][j] /= membershipSum;
        }
    }

    
    std::vector<std::vector<double>> oldMembershipMatrix(pointCount, std::vector<double>(clusterCount, 0.0));
    do {
        
        oldMembershipMatrix = membershipMatrix;

        
        updateMembershipMatrix(data, centroids, membershipMatrix, fuzziness);

        
        updateCentroids(data, membershipMatrix, centroids);

    } while (!hasConverged(oldMembershipMatrix, membershipMatrix, epsilon));
}

 
void input(std::vector<std::vector<double>>& data, size_t& clusterCount, double& fuzziness)
{
    size_t dimensionCount, pointCount;
    std::cout << "Enter the number of dimensions (n):\n";
    std::cin >> dimensionCount;
    std::cout << "Enter the number of randomly generated points:\n";
    std::cin >> pointCount;
    std::cout << "Enter the number of clusters:\n";
    std::cin >> clusterCount;
    if (clusterCount >= pointCount) std::cout << "The values entered may cause errors!\n";
    do
    {
        std::cout << "Enter the \"fuzziness\" parameter (>=1). Recommended values are in the range [1.5, 2.5]:\n";
        std::cin >> fuzziness;
    } while (fuzziness < 1);

    std::cout << "The upper asymptotic complexity of this Fuzzy C-means algorithm is approximately O(K*N*D*C), where:\n";
    std::cout << "K is the number of iterations until convergence,\n";
    std::cout << "N is the number of points,\n";
    std::cout << "C is the number of clusters,\n";
    std::cout << "D is the number of dimensions.\n";
    std::cout << "For large input values, the data processing time will be longer.\n\n";

    std::ofstream generatedPointsFile("generatedPoints.txt");
    for (int i = 0; i < pointCount; i++)
    {
        std::vector<double> point;
        for (int j = 0; j < dimensionCount; j++)
        {
            double coordinate = ((std::rand()) % 1000);
            coordinate = double(coordinate / 100);
            point.push_back(coordinate);
            generatedPointsFile << coordinate << " ";
        }
        data.push_back(point);
        generatedPointsFile << '\n';
    }
    generatedPointsFile.close();
}

void printResult(std::vector<std::vector<double>> data, std::vector<std::vector<double>> centroids, size_t clusterCount, std::vector<std::vector<double>> membershipMatrix)
{
    std::ofstream resultFile("results.txt");
    std::cout << "Centroids:\n";
    resultFile << "Centroids:\n";
    for (size_t i = 0; i < clusterCount; ++i) {
        std::cout << "Cluster " << i + 1 << ": ";
        resultFile << "Cluster " << i + 1 << ": ";
        for (size_t d = 0; d < data[0].size(); ++d) {
            std::cout << centroids[i][d] << " ";
            resultFile << centroids[i][d] << " ";
        }
        std::cout << std::endl;
        resultFile << std::endl;
    }

    std::cout << "\n\"Membership\" Matrix:\n";
    resultFile << "\n\"Membership\" Matrix:\n";
    for (size_t i = 0; i < data.size(); ++i) {
        std::cout << "Point " << i + 1 << ": ";
        resultFile << "Point " << i + 1 << ": ";
        for (size_t j = 0; j < clusterCount; ++j) {
            std::cout << membershipMatrix[i][j] << " ";
            resultFile << membershipMatrix[i][j] << " ";
        }
        std::cout << std::endl;
        resultFile << std::endl;
    }
    resultFile.close();
}


int main() {
    std::srand(std::time(0)); 

    std::vector<std::vector<double>> data;
    size_t clusterCount;
    double fuzziness; 
    double epsilon = 1e-6;  
    input(data, clusterCount, fuzziness);

    std::vector<std::vector<double>> centroids;
    std::vector<std::vector<double>> membershipMatrix;

    fuzzyCMeans(data, clusterCount, fuzziness, epsilon, centroids, membershipMatrix);

    
    printResult(data, centroids, clusterCount, membershipMatrix);

    return 0;
}
