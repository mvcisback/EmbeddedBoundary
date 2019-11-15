#ifndef EMBEDDED_BOUNDARY_GEOMETRY_BOUNDARY_H
#define EMBEDDED_BOUNDARY_GEOMETRY_BOUNDARY_H

#include "inputs/input_base.h"

#include <array>
#include <map>
#include <vector>

namespace boundary {

namespace geometry {

  struct geo_info {
    /// tells if cell is irregular or not
    bool irregular;
    /// id to index into vector with more cell information
    int id;
    /// derivatives of the normal to the boundary
    std::vector<std::vector<std::vector<double>>> normal_derivatives;
    /// 1d volume fraction for boundary cell edges
    std::array<double, 4> vol_frac_1d;
    /// volume moments
    std::vector<std::vector<double>> volume_moments;
    /// boundary moments
    std::vector<std::vector<double>> boundary_moments;
    // children
    std::map<std::array<double, 2>, geo_info> children;
  };

  struct cell_info{
    /// tells if cell is inside (0), outside (1), or boundary (2)
    int cell_type;
    // children
    std::map<std::array<double, 2>, geo_info> children;
  };

/// A class describing the boundary geometry.
/**
This class stores a map of all boundary cells with necessary geometry information
*/
  class Boundary{
    public:
      /// Constructor
      /**
      Iterates through all cells in geometry and adds boundary cells to map
      */
      Boundary(boundary::inputs::InputBase* input);
      /// Determines if a cell is a boundary cell. Deprecated -- go through and replace with inside outside
      static bool IsBoundaryCell(helpers::Point lower_left,
                                 helpers::Point lower_right,
                                 helpers::Point upper_right,
                                 helpers::Point upper_left,
                                 boundary::inputs::InputBase* input);
      /// Determines if cell corners are all inside (0), all outside (1), or mixed (2)
      static int InsideOutside(helpers::Point lower_left,
                               helpers::Point lower_right,
                               helpers::Point upper_right,
                               helpers::Point upper_left,
                               boundary::inputs::InputBase* input);
      std::map<std::array<double, 2>, geo_info> BoundaryCells();
      static double WhichValue(std::vector<double> values,
                         double first_bound,
                         double second_bound);
      void MakeCellTree();
      void RecursiveRefine(double x_min, double x_max, double y_min,
                           double y_max);
      std::map<helpers::Point, cell_info> CellTree();                    
    private:

      void CalculateMoments_(std::array<double, 2> cell_center);
      double DIntegral_(double beginning,
                        double end,
                        std::array<int, 2> q,
                        int index,
                        double fixed_value);
      double CalcD_(double bd_length,
                    double fixed_value,
                    std::array<double, 2> cell_center,
                    std::array<int, 2> q,
                    int d,
                    std::array<int, 2> which_d);
      std::map<std::array<double, 2>, geo_info> boundary_cells_;
      boundary::inputs::InputBase* input_;
      int Q_;
      double cell_size_;
      std::map<helpers::Point, cell_info> cell_tree_;
  };
}

}
#endif // EMBEDDED_BOUNDARY_GEOMETRY_BOUNDARY_H
