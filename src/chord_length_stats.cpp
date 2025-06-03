#include "openmc/chord_length_stats.h"

#include "openmc/capi.h"
#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/geometry.h"
#include "openmc/hdf5_interface.h"
#include "openmc/material.h"
#include "openmc/message_passing.h"
#include "openmc/openmp_interface.h"
#include "openmc/output.h"
#include "openmc/random_lcg.h"
#include "openmc/settings.h"
#include "openmc/timer.h"
#include "openmc/xml_interface.h"

#include "xtensor/xadapt.hpp"
#include "xtensor/xview.hpp"
#include <fmt/core.h>

#include <algorithm> // for copy
#include <cmath>     // for pow, sqrt
#include <unordered_set>

namespace openmc {
//==============================================================================
// Global variables
//==============================================================================
namespace model {
vector<ChordLengthStats> chordl_stats;
}
//==============================================================================
// Chord Length Statistic implementation
//==============================================================================
ChordLengthStats::ChordLengthStats(pugi::xml_node node)
{
  // Read domain type (material)
  std::string domain_type = get_node_value(node, "domain_type");
  if (domain_type == "material") {
    domain_type_ = TallyDomain::MATERIAL;
  } else {
    fatal_error(std::string("Unrecognized domain type for chord length "
                            "statistics: " +
                            domain_type));
  }

  // Read domain IDs and number of samples
  matrix_domain_ids_ = get_node_array<int>(node, "matrix_domain_ids");
  stochastic_media_domain_ids_ =
    get_node_array<int>(node, "stochastic_media_domain_ids");
  n_samples_ = std::stoull(get_node_value(node, "samples"));
  tally_bins = get_node_array<double>(node, "tally_bins");

  // Ensure that the matrix material and the stochastic media region
  // correspond to each other and that they are not identical
  if (matrix_domain_ids_.size() != stochastic_media_domain_ids_.size()) {
    throw std::runtime_error(
      "Matrix domain IDs and stochastic media domain IDs must have the same "
      "length.");
  }
  for (size_t i = 0; i < matrix_domain_ids_.size(); ++i) {
    if (matrix_domain_ids_[i] == stochastic_media_domain_ids_[i]) {
      throw std::runtime_error(
        fmt::format("Matrix domain ID {} and stochastic media domain ID "
                    "{} cannot be the same.",
          matrix_domain_ids_[i], stochastic_media_domain_ids_[i]));
    }
  }
  // Ensure that the domain IDs are unique
  std::unordered_set<int> unique_matrix_ids(
    matrix_domain_ids_.cbegin(), matrix_domain_ids_.cend());
  if (unique_matrix_ids.size() != matrix_domain_ids_.size()) {
    throw std::runtime_error {"Matrix domain IDs for chord length statistics "
                              "must be unique."};
  }
  std::unordered_set<int> unique_stochastic_ids(
    stochastic_media_domain_ids_.cbegin(), stochastic_media_domain_ids_.cend());
  if (unique_stochastic_ids.size() != stochastic_media_domain_ids_.size()) {
    throw std::runtime_error {"Stochastic media domain IDs for chord length "
                              "statistics must be unique."};
  }
}
vector<ChordLengthStats::Result> ChordLengthStats::execute() const
{
  // Check to make sure domain IDs are valid
  for (int id : matrix_domain_ids_) {
    if (model::material_map.find(id) == model::material_map.end()) {
      throw std::runtime_error {fmt::format(
        "Matrix material {} in chord length statistics does not exist in "
        "geometry.",
        id)};
    }
  }
  for (int id : stochastic_media_domain_ids_) {
    if (model::material_map.find(id) == model::material_map.end()) {
      throw std::runtime_error {fmt::format(
        "Stochastic media material {} in chord length statistics does not "
        "exist in geometry.",
        id)};
    }
  }

  while (true) {
    Particle p;
    // Initialize particle with a random position and direction
    int64_t id = p.id();
    uint64_t seed = init_seed(id, STREAM_VOLUME);
    p.n_coord() = 1;
    Position xi {prn(&seed), prn(&seed), prn(&seed)};
    p.r() = lower_left_ + xi * (upper_right_ - lower_left_);
    p.u() = {prn(&seed), prn(&seed), prn(&seed)};

    // Initialization result for chord length statistics
    Result result;
    result.chord_length.resize(tally_bins.size());
    for (auto& chord : result.chord_length) {
      chord = {0.0, 0.0}; // Initialize mean and std deviation
    }

    if (!exhaustive_find_cell(p))
      continue;
    double total_chord_length = 0.0;
    while (p.alive()) {
      p.boundary() = distance_to_boundary(p);
      double distance = p.boundary().distance;

      // Advance particle in space
      // Short-term solution until the surface source is revised and we can use
      // this->move_distance(distance)
      for (int j = 0; j < p.n_coord(); ++j) {
        p.coord(j).r += distance * p.coord(j).u;
      }
      p.event_cross_surface();

      if (p.material() == p.material_last()) {
        total_chord_length += distance;
      } else {
        // If the material has changed, we need to record the chord length
        // statistics for the chord length in matrix to stochastic media
        if (p.material_last() != C_NONE) {
          // Find the index of the matrix material
          auto it = std::find(matrix_domain_ids_.begin(),
            matrix_domain_ids_.end(), p.material_last());
          if (it != matrix_domain_ids_.end()) {
            size_t index = std::distance(matrix_domain_ids_.begin(), it);
            // Update the chord length statistics for this domain
            double chord_length = total_chord_length;
            double mean = chord_length / n_samples_;
            double variance =
              chord_length * chord_length / n_samples_ - mean * mean;
            result.chord_length[index][0] += mean;
            result.chord_length[index][1] += variance;
          }
        }
        // Reset the total chord length for the next material
        total_chord_length = 0.0;
      }

      // If the particle is outside the bounding box, kill it
      if (p.r().x < lower_left_.x || p.r().y < lower_left_.y ||
          p.r().z < lower_left_.z || p.r().x > upper_right_.x ||
          p.r().y > upper_right_.y || p.r().z > upper_right_.z) {
        p.wgt() = 0;
      }
    }
    // If the material is not changed, cumulative chord length
  }
}
} // namespace openmc