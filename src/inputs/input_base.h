#ifndef EMBEDDED_BOUNDARY_INPUT_INPUT_BASE_H
#define EMBEDDED_BOUNDARY_INPUT_INPUT_BASE_H

#include "helpers/geometry_objects.h"

#include <vector>
#include <array>

namespace boundary {

namespace inputs {

/**
The input file base class. All input files are based on this class which
contains all relevant input geometry information.
*/
class InputBase {
public:
  virtual ~InputBase() = default;
  /// The function describing the boundary.
  virtual double BoundaryFunction(double x_value) = 0;
  /// The derivatives of the boundary function.
  virtual double BoundaryDerivatives(helpers::Point a_point, std::vector<int> degree) = 0;
  /// The inverse of the boundary function.
  virtual double BoundaryInverse(double y_value) = 0;
  /// A function returning whether or not a point is inside or outside the
  /// boundary.
  /* It returns 0 for outside, 1 for inside, 2 for on boundary.
  */
  virtual int Inside(std::array<double, 2> point) = 0;
  /// Domain minimum in the x direction.
  virtual double XMin() = 0;
  /// Domain maximum in the x direction.
  virtual double XMax() = 0;
  /// Domain minimum in the y direction.
  virtual double YMin() = 0;
  /// Domain maximum in the y direction.
  virtual double YMax() = 0;
  /// Size of cells (assumes square cells).
  virtual double CellSize() = 0;
  /// Desired Q order.
  virtual int QOrder() = 0;

};

} // namespace inputs

} // namespace boundary

#endif // EMBEDDED_BOUNDARY_INPUT_INPUT_BASE_H
