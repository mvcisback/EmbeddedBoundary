#ifndef EMBEDDED_BOUNDARY_SRC_INPUTS_OFFSET_CIRCLE_OFFSET_CIRCLE_H
#define EMBEDDED_BOUNDARY_SRC_INPUTS_OFFSET_CIRCLE_OFFSET_CIRCLE_H

#include "inputs/input_base.h"
#include "helpers/geometry_objects.h"

#include <vector>
#include <array>

namespace boundary {

namespace inputs {


/// An input file representing the Circle 0 = (x-0.125)^2 + y^2 - 1.005
class OffsetCircle : public InputBase{
  public:
    ~OffsetCircle() = default;
    std::vector<double> BoundaryFunction(double x_value) override;
    double BoundaryDerivatives(helpers::Point a_point, std::vector<int> degree) override;
    std::vector<double> BoundaryInverse(double y_value) override;
    int Inside(helpers::Point point) override;
    double XMin() override;
    double XMax() override;
    double YMin() override;
    double YMax() override;
    double CellSize() override;
    int QOrder() override;

};

} // namespace inputs

} // namespace boundary

#endif // EMBEDDED_BOUNDARY_SRC_INPUTS_OFFSET_CIRCLE_OFFSET_CIRCLE_H
