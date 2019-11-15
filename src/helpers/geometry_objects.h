#ifndef EMBEDDED_BOUNDARY_HELPERS_GEOMETRY_OBJECTS_H
#define EMBEDDED_BOUNDARY_HELPERS_GEOMETRY_OBJECTS_H

#include <vector>
#include <cmath>

namespace boundary {

namespace helpers {

/// A Class describing a point
/** This class represents a two dimensional point
*/
  class Point{
    public:
      Point();
      Point(double first_dim, double second_dim);
      ~Point() = default;
      Point operator + (const Point &a_point);
      /// Lexicographic Ordering
      bool operator < (const Point &a_point) const;
      double x_val;
      double y_val;
      //int MortonOrder(); will implement later using Phil's code
  };

} // namespace helpers

} // namespace boundary

#endif // EMBEDDED_BOUNDARY_HELPERS_GEOMETRY_OBJECTS_H
