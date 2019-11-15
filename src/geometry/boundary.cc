#include "geometry/boundary.h"
#include "helpers/geometry_objects.h"
#include "normals/normals.h"

#include <array>
#include <iostream>
#include <stdexcept>
#include <Eigen/Dense>

namespace boundary {

namespace geometry {

Boundary::Boundary(boundary::inputs::InputBase* input){
  input_ = input;
  cell_size_ = input->CellSize();
  Q_ = input->QOrder();
  double y_min = input->YMin();
  double y_max = y_min + cell_size_;
  int id_count = 0;
  while (y_max <= input->YMax()){
    double x_min = input->XMin();
    double x_max = x_min + cell_size_;
    while (x_max <= input->XMax()){
      // corners of cell
      std::array<helpers::Point, 4> corners;
      corners[0] = helpers::Point(x_min, y_min); // lower left
      corners[1] = helpers::Point(x_min, y_max); // upper left
      corners[2] = helpers::Point(x_max, y_max); // upper right
      corners[3] = helpers::Point(x_max, y_min); // lower right

      bool is_boundary = IsBoundaryCell(corners[0], corners[1],
                                        corners[2], corners[3], input_);
      if (is_boundary){
        // make struct with tag and id
        geo_info cell;
        cell.irregular = true;
        cell.id = id_count;
        // update id count
        id_count += 1;
        helpers::Point cell_center = helpers::Point(x_max - cell_size_/2,
                                                    y_max - cell_size_/2);
        // resize vectors
        cell.normal_derivatives.resize(Q_+1);
        cell.volume_moments.resize(Q_+1);
        cell.boundary_moments.resize(Q_+1);
        for (int i=0; i < Q_+1; ++i){
          cell.normal_derivatives[i].resize(Q_+1);
          cell.volume_moments[i].resize(Q_+1);
          cell.boundary_moments[i].resize(Q_+1);
          for (int j = 0; j < Q_ + 1; ++j){
            cell.normal_derivatives[i][j].resize(2);
          }
        }
        // Compute derivatives of the normal and store in cell struct
        for (int q_mag=Q_; q_mag >= 0 ; q_mag--){
          for (int q1 = 0; q1 < q_mag + 1; q1++){
            std::vector<int> q = {q1, q_mag - q1};
            cell.normal_derivatives[q1][q_mag - q1][0] = normals::Normal::NormalDerivative(q,
                                                            1,
                                                            cell_center,
                                                            input);
            cell.normal_derivatives[q1][q_mag - q1][1] = normals::Normal::NormalDerivative(q,
                                                            2,
                                                            cell_center,
                                                            input);
          }
        }
        // Compute 1d volume fractions and store in cell
        // Check if four corners are inside or outside boundary
        std::vector<int> inside{input_->Inside(corners[0]),
                                input_->Inside(corners[1]),
                                input_->Inside(corners[2]),
                                input_->Inside(corners[3])};
        for (int corner=0; corner < 4; corner++){
          // iterate four edges to determine which ones intersect boundary
          if ((inside[corner] == inside[(corner+1)%4]) || inside[corner] == 2){ // no change from corner to corner, TODO: come up with more robust way to do this
            cell.vol_frac_1d[corner] = cell_size_*inside[(corner+1)%4];
          }
          else if (inside[(corner+1)%4] == 2){ // next corner is on the boundary
            cell.vol_frac_1d[corner] = cell_size_*inside[corner];
          }
          else { // change from corner to corner, indicating boundary cuts through
            // check if edge is horizontal or vertical
            if ((corners[corner].x_val - corners[(corner+1)%4].x_val) == 0){ // horizontal
              // find the intersection of x=corners[i, 0]
              std::vector<double> y_values = input->BoundaryFunction(corners[corner].x_val);
              double y;
              if (y_values.size() == 1){
                y = input->BoundaryFunction(corners[corner].x_val)[0];
              }
              else if (y_values.size() == 2){
                // choose the value that is in the cell
                // Robustness TODO: fix this in the case that both values are in cell
                y = WhichValue(y_values, corners[corner].y_val, corners[(corner+1)%4].y_val);
              }

              double edge_inside = abs(y - corners[corner].y_val);
              if (inside[corner]){
                cell.vol_frac_1d[corner] = edge_inside;
              }
              else{
                cell.vol_frac_1d[corner] = cell_size_ - edge_inside;
              }
            }
            else { // vertical
              // find intersection of y=corners[i, 1]
              std::vector<double> x_values = input->BoundaryInverse(corners[corner].y_val);
              double x;
              if (x_values.size() == 1){
                x = input->BoundaryFunction(corners[corner].x_val)[0];
              }
              else if (x_values.size() == 2){
                // choose the value that is in the cell
                // Robustness TODO: fix this in the case that both values are in cell
                x = WhichValue(x_values, corners[corner].x_val, corners[(corner+1)%4].x_val);
              }
              double edge_inside = abs(x - corners[corner].x_val);
              if (inside[corner]){
                cell.vol_frac_1d[corner] = edge_inside;
              }
              else {
                cell.vol_frac_1d[corner] = cell_size_ - edge_inside;
              }
            }
          }
        }
        // add point to map
        std::array<double, 2> point = {x_min + cell_size_/2, y_min + cell_size_/2};
        boundary_cells_.insert(
          std::pair<std::array<double, 2>, geo_info>(point, cell));
      }
      // Go to next cell
      x_max += cell_size_;
      x_min += cell_size_;
    }
    // Go to next row
    y_max += cell_size_;
    y_min += cell_size_;
  }

  for (std::map<std::array<double, 2>, geo_info>::iterator it=boundary_cells_.begin();
       it != boundary_cells_.end(); it++){
          CalculateMoments_(it->first);
       }

};


void Boundary::MakeCellTree(){
  double y_min = input_->YMin();
  double y_max = y_min + cell_size_;
  while (y_max <= input_->YMax()){
    double x_min = input_->XMin();
    double x_max = x_min + cell_size_;
    while (x_max <= input_->XMax()){
      RecursiveRefine(x_min, x_max, y_min, y_max);
      }
      // Go to next cell
      x_max += cell_size_;
      x_min += cell_size_;
    }
    // Go to next row
    y_max += cell_size_;
    y_min += cell_size_;
};

void Boundary::RecursiveRefine(double x_min, double x_max, double y_min,
  double y_max){
  double cell_size = std::abs(x_max - x_min);
  // corners of cell
  std::array<helpers::Point, 4> corners;
  corners[0] = helpers::Point(x_min, y_min); // lower left
  corners[1] = helpers::Point(x_min, y_max); // upper left
  corners[2] = helpers::Point(x_max, y_max); // upper right
  corners[3] = helpers::Point(x_max, y_min); // lower right

  helpers::Point cell_center = helpers::Point(x_min + cell_size/2,
                                              y_min + cell_size/2);

  // boundary check
  int cell_inside = InsideOutside(corners[0], corners[1], corners[2],
                                  corners[3], input_);
  cell_info cell;

  // if not boundary
  if (cell_inside == 0 || cell_inside == 1){
    // gradient test at each corner
    double corner_value;
    for (int i = 0; i < 4; i++){
      corner_value = input_->BoundaryDerivatives(corners[i],
                                                 std::vector<int> {0, 0});
      if (corner_value < std::sqrt(2)*cell_size_){
        // refine recursively by splitting cell into four
        RecursiveRefine(x_min, cell_center.x_val, y_min, cell_center.y_val); // lower left
        RecursiveRefine(x_min, cell_center.x_val, cell_center.y_val, y_max); // upper left
        RecursiveRefine(cell_center.x_val, x_max, cell_center.y_val, y_max); // upper right
        RecursiveRefine(cell_center.x_val, x_max, y_min, cell_center.y_val); // lower right
        // don't bother checking the rest of the corners
      }
      else {
        cell.cell_type = cell_inside;
        // add point to map
        cell_tree_.insert(std::pair<helpers::Point, cell_info> (cell_center, cell));
      }
    }
  }

  // if boundary, mark as boundary, go to next cell
  else {
    cell.cell_type = 2;
    // add point to map
    cell_tree_.insert(std::pair<helpers::Point, cell_info> (cell_center, cell));
 }
};

double Boundary::WhichValue(std::vector<double> values, double first_bound, double second_bound){
  if (first_bound < second_bound){
    for (std::vector<double>::iterator it = values.begin(); it != values.end(); it++){
      if (*it >= first_bound && *it <= second_bound){
        return *it;
      }
    }
  }

  else if (first_bound > second_bound){
    for (std::vector<double>::iterator it = values.begin(); it != values.end(); it++){
      if (*it <= first_bound && *it >= second_bound){
        return *it;
      }
    }
  }

  else {
    throw std::invalid_argument("No change between bounds.");
  }
};


bool Boundary::IsBoundaryCell(helpers::Point lower_left,
                              helpers::Point lower_right,
                              helpers::Point upper_right,
                              helpers::Point upper_left,
                              boundary::inputs::InputBase* input){
  // put into array with inside values
  std::vector<int> inside{input->Inside(lower_left),
                          input->Inside(upper_left),
                          input->Inside(upper_right),
                          input->Inside(lower_right)};
  // try changing all 2s to 0s and see if they are all the same
  int it_prev = inside[0];
  bool same;
  for (int it = 1; it < 4; it++){
    int current = inside[it];
    if (current == 2){
      current = 0;
    }
    if (it_prev == current){
      if (it == 3){
        // they're all the same, not a boundary cell
        same = true;
      }
      continue;
    }
    else{
      same = false;
      break;
    }
  }
  // if they're not try changing all the 2s to 1s and see if they're all the same
  if (!same){
    it_prev = inside[0];
    for (int it = 1; it < 4; it++){
      int current = inside[it];
      if (current == 2){
        current = 1;
      }
      if (it_prev == current){
        if (it == 3){
          same = true;
        }
        continue;
      }
      else{
        same = false;
        break;
      }
    }
  }
  return !same;
}

int Boundary::InsideOutside(helpers::Point lower_left,
                            helpers::Point lower_right,
                            helpers::Point upper_right,
                            helpers::Point upper_left,
                            boundary::inputs::InputBase* input){
  // put into array with inside values
  std::vector<int> inside{input->Inside(lower_left),
                          input->Inside(upper_left),
                          input->Inside(upper_right),
                          input->Inside(lower_right)};
  // try changing all 2s to 0s and see if they are all the same
  int it_prev = inside[0];
  int status;
  for (int it = 1; it < 4; it++){
    int current = inside[it];
    if (current == 2){
      current = 0;
    }
    if (it_prev == current){
      if (it == 3){
        // they're all the same, not a boundary cell
        status = current;
      }
      continue;
    }
    else{
      status = 2;
      break;
    }
  }
  // if they're not try changing all the 2s to 1s and see if they're all the same
  if (status == 2){
    it_prev = inside[0];
    for (int it = 1; it < 4; it++){
      int current = inside[it];
      if (current == 2){
        current = 1;
      }
      if (it_prev == current){
        if (it == 3){
          status = current;
        }
        continue;
      }
      else{
        status = 2;
        break;
      }
    }
  }
  return status;
}


std::map<helpers::Point, cell_info> Boundary::CellTree(){
  return cell_tree_;
}


std::map<std::array<double, 2>, geo_info> Boundary::BoundaryCells(){
  return boundary_cells_;
}


double Boundary::DIntegral_(double beginning, double end, std::array<int, 2> q,
                            int index, double fixed_value){
  double next_plus = q[(index + 1) % 2] + 1;
  return std::pow(fixed_value, q[index])*(1/(next_plus))*(std::pow(end, next_plus) - std::pow(beginning, next_plus));
}


double Boundary::CalcD_(double bd_length,
              double fixed_value,
              std::array<double, 2> cell_center,
              std::array<int, 2> q,
              int d,
              std::array<int, 2> which_d){
  int d_op = (d + 1) % 2;
  double d_pm;
  if (bd_length == cell_size_){
    d_pm = DIntegral_(cell_center[d_op] - cell_size_/2,
                            cell_center[d_op] + cell_size_/2,
                            q, d, fixed_value);
  }
  else if (bd_length > 0){
    // TODO: fix this to be more general
    helpers::Point corner = helpers::Point(cell_center[0] + cell_size_/2*which_d[0],
                                         cell_center[1] + cell_size_/2*which_d[1]);
    std::vector<double> int_values = input_->BoundaryInverse(fixed_value);
    double intersection;
    if (int_values.size() == 1){
      intersection = int_values[0];
    }
    else if (int_values.size() == 2){
      // choose the value that is in the cell
      // Robustness TODO: fix this in the case that both values are in cell
      intersection = WhichValue(int_values, cell_center[d_op] - cell_size_/2, cell_center[d_op] + cell_size_/2);
    }
    if (input_->Inside(corner)){
          if (d_op == 0){
            d_pm = DIntegral_(intersection, corner.x_val, q, d, fixed_value);
          }
          else {
            d_pm = DIntegral_(intersection, corner.y_val, q, d, fixed_value);
          }

    }
    else {
      if (d_op == 0){
        d_pm = DIntegral_(corner.x_val - cell_size_, intersection, q, d, fixed_value);
      }
      else {
        d_pm = DIntegral_(corner.y_val - cell_size_, intersection, q, d, fixed_value);
      }
    }
  }
  else {
    d_pm = 0;
  }
  return d_pm;
}


void Boundary::CalculateMoments_(std::array<double, 2> cell_center){
  for (int q_mag = Q_; q_mag >= 0; q_mag--){
    Eigen::VectorXf rho = Eigen::VectorXf::Zero(2*(q_mag + 1));
    Eigen::MatrixXf lhs = Eigen::MatrixXf::Zero(2*(q_mag + 1), 2*q_mag + 1);

    int it = -1;
    for (int q1 = 0; q1 < (q_mag + 1); q1++){
      std::array<int, 2> q = {q1, q_mag - q1};
      for (int d = 0; d < 2; d++){
        it += 1;
        std::array<int, 2> e = {0, 0};
        e[d] = 1;

        std::array<int, 2> plus_vec = {1, 1};
        std::array<int, 2> first_minus_vec = {-1, 1};
        std::array<int, 2> second_minus_vec = {1, -1};

        double ptv_bd_length;
        double d_plus;
        double ntve_bd_length;
        double d_minus;
        // set d_plus and d_minus
        if (d == 0){
          // calculate d_plus
          ptv_bd_length = boundary_cells_[cell_center].vol_frac_1d[2];
          double fixed_ptve_x = cell_center[0] + cell_size_/2;
          d_plus = CalcD_(ptv_bd_length, fixed_ptve_x, cell_center, q, d, plus_vec);

          // calculate d_minus
          ntve_bd_length = boundary_cells_[cell_center].vol_frac_1d[0];
          double fixed_ntve_x = cell_center[0] - cell_size_/2;
          d_minus = CalcD_(ntve_bd_length, fixed_ntve_x, cell_center, q,
                           d, first_minus_vec);
        }
        else if (d == 1){
          // calculate d_plus
          ptv_bd_length = boundary_cells_[cell_center].vol_frac_1d[1];
          double fixed_ptve_y = cell_center[1] + cell_size_/2;
          d_plus = CalcD_(ptv_bd_length, fixed_ptve_y, cell_center, q,
                                 d, plus_vec);

          // calculate d_minus
          ntve_bd_length = boundary_cells_[cell_center].vol_frac_1d[3];
          double fixed_ntve_y = cell_center[1] - cell_size_/2;
          d_minus = CalcD_(ntve_bd_length, fixed_ntve_y, cell_center, q,
                           d, second_minus_vec);
        }
        int s_sum = 0;
        int S = Q_ - q_mag;
        for (int s1 = 0; s1 < S+1; s1++){
          for (int s2 = 0; s2 < S+1; s2++){
            if (s1 + s2 >= S || s1 + s2 < 1){
              continue;
            }
            s_sum += (boundary_cells_[cell_center].normal_derivatives[s1][s2][d]
                      *boundary_cells_[cell_center].boundary_moments[q[0] + s1][q[1] + s2]);

          }
        }
        std::array<int, 2> q_minus_e = {q[0] - e[0], q[1] - e[1]};
        if (q_minus_e[0] >= 0 && q_minus_e[1] >= 0){
          lhs(it, q_minus_e[0]) = q[d];
        }
        // M_B Unknowns (stored based on value of q[0
        lhs(it, q_mag + q[0]) = -boundary_cells_[cell_center].normal_derivatives[0][0][d];
        rho(it) = d_plus - d_minus + s_sum;
      }
    }
    if (cell_center[0] == .625 && cell_center[1] == -.875){
    }
    // Solve
    Eigen::VectorXf v_and_b = lhs.colPivHouseholderQr().solve(rho);

    // Unpack V and B
    for (int i = 0; i < q_mag; i++){
      boundary_cells_[cell_center].volume_moments[i][q_mag - 1 - i] = v_and_b(i);
    }
    for (int i = 0; i < q_mag + 1; i++){
      boundary_cells_[cell_center].boundary_moments[i][q_mag - i] = v_and_b(q_mag + i);
    }
  }
}


} // namespace geometry

} // namespace boundary
