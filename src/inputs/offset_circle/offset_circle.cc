#include "inputs/offset_circle/offset_circle.h"
#include "helpers/geometry_objects.h"

#include <vector>
#include <iostream>

/// A circle input designed to have a boundary that enters and exits on the same cell edge
/// Designed to test refinement criteria

namespace boundary {

namespace inputs {

std::vector<double> OffsetCircle::BoundaryFunction(double x_value){
  std::vector<double> y_values;
  y_values.push_back(std::sqrt(1.005 - std::pow(x_value - 0.125, 2)));
  y_values.push_back(-std::sqrt(1.005 - std::pow(x_value - 0.125, 2)));
  return y_values;
};

double OffsetCircle::BoundaryDerivatives(helpers::Point a_point,
                                         std::vector<int> degree){
  if (degree[0] == 0 && degree[1] == 0){
    return std::pow(a_point.y_val, 2) + std::pow(a_point.x_val - 0.125, 2) - 1.005;
  }

  else if (degree[0] + degree[1] == 1){
    return 2*(degree[0]*(a_point.x_val - 0.125) + degree[1]*a_point.y_val);
  }

  else if ((degree[0] == 2 && degree[1] == 0) || (degree[0] == 0 && degree[1] == 2)){
    return 2;
  }

  else {
    return 0;
  }
};


std::vector<double> OffsetCircle::BoundaryInverse(double y_value){
  std::vector<double> x_values;
  x_values.push_back(std::sqrt(1.005 - std::pow(y_value, 2)) + 0.125);
  x_values.push_back(-std::sqrt(1.005 - std::pow(y_value, 2)) + 0.125);
  return x_values;
};


int OffsetCircle::Inside(helpers::Point point){
  if ((std::pow(point.x_val - 0.125, 2) + std::pow(point.y_val, 2)) > 1.005){
    return 0;
  }
  else if ((std::pow(point.x_val - 0.125, 2) + std::pow(point.y_val, 2)) < 1.005){
    return 1;
  }
  else {
    return 2;
  }
};


double OffsetCircle::XMin(){
  return -1.25; // domain min
};


double OffsetCircle::XMax(){
  return 1.25; // domain max
};


double OffsetCircle::YMin(){
  return -1.25; // domain min
};


double OffsetCircle::YMax(){
  return 1.25; // domain max
};


double OffsetCircle::CellSize(){
  return .25;
};

int OffsetCircle::QOrder(){
  return 1; // Q order
}

} // namespace inputs

} // namespace boundary
