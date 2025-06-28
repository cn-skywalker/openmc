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
  matrix_domain_id_ = std::stoull(get_node_value(node, "matrix_domain_id"));
  stochastic_media_domain_id_ =
    std::stoull(get_node_value(node, "stochastic_media_domain_id"));
  n_samples_ = std::stoull(get_node_value(node, "samples"));
  tally_bins_ = get_node_array<double>(node, "tally_bins");
  lower_left_ = get_node_array<double>(node, "lower_left");
  upper_right_ = get_node_array<double>(node, "upper_right");
}
ChordLengthStats::Result ChordLengthStats::execute() const
{
  // Check to make sure domain IDs are valid

  if (model::material_map.find(matrix_domain_id_) ==
      model::material_map.end()) {
    throw std::runtime_error {fmt::format(
      "Matrix material {} in chord length statistics does not exist in "
      "geometry.",
      matrix_domain_id_)};
  }

  if (model::material_map.find(stochastic_media_domain_id_) ==
      model::material_map.end()) {
    throw std::runtime_error {fmt::format(
      "Stochastic media material {} in chord length statistics does not "
      "exist in geometry.",
      stochastic_media_domain_id_)};
  }

  // Initialize the result structure
  Result result;
  result.chord_length.resize(tally_bins_.size() - 1);
  for (auto& chord : result.chord_length) {
    chord = {0.0}; //  Initialize frequency
  }

  // Create a local result for each thread to avoid race conditions
  std::vector<std::vector<double>> local_results(
    omp_get_max_threads(), std::vector<double>(tally_bins_.size() - 1, 0.0));
// Parallelize the particle loop
#pragma omp parallel for
  for (int64_t i = 0; i < n_samples_; ++i) {
    // Get the thread ID for local result access
    int thread_id = omp_get_thread_num();
    //  Generate a particle with random position and direction, initialize
    //  current chord length to 0
    // Initialize a particle with a unique ID
    Particle p;
    p.id() = i; // Set particle ID
    uint64_t seed = init_seed(i, STREAM_VOLUME);
    bool if_first = true;
    Position xi {prn(&seed), prn(&seed), prn(&seed)};
    p.r() = lower_left_ + xi * (upper_right_ - lower_left_);
    // Generate a random isotropic direction and set it to p.u()
    double theta = 2.0 * M_PI * prn(&seed); // Random azimuthal angle [0, 2π)
    double phi = std::acos(2.0 * prn(&seed) - 1.0); // Random polar angle [0, π]

    // Convert spherical coordinates to Cartesian coordinates
    p.u() = {std::sin(phi) * std::cos(theta), std::sin(phi) * std::sin(theta),
      std::cos(phi)};

    p.r_born() = p.r();            // Set the initial position of the particle
    double total_chord_length = 0; // Initialize total chord length

    // If the particle is not within the geometry, skip it
    if (!exhaustive_find_cell(p))
      continue;

    while (p.alive()) {
      //  Find the position  of the next face along the initial flight direction
      //  and pass through it, and determine whether it exceeds the boundary of
      //  the area.
      p.boundary() = distance_to_boundary(p);
      double distance = p.boundary().distance;

      //  Move the particle to the next position
      for (int j = 0; j < p.n_coord(); ++j) {
        p.coord(j).r += distance * p.coord(j).u;
      }
      p.event_cross_surface();
      total_chord_length += distance;

      // Determine whether it exceeds the regional boundary
      if (check_hit_boundary(p)) {
        p.wgt() = 0; //  Kill the particle
        break;       // Exit the current particle loop and regenerate particles.
      }

      if (model::material_map[stochastic_media_domain_id_] ==
          p.material_last()) {
        // Set if_first to false after the first material switch
        if_first = false;
      }

      // Determine whether the material has been switched
      if (p.material() != p.material_last()) {

        // The material is switched, and it is determined whether it matches
        // both the base material and the random medium material.

        if (check_material_match(p.material(), p.material_last()) &&
            !if_first) {
          // If matched, use match_index.value() to obtain the index.
          // Find the index of the tally bin for the current chord length
          auto index = find_tally_bin_index(total_chord_length);

          // Add the current chord length to the corresponding frequency
          // statistics.
          if (index.has_value()) {
            size_t tally_index = index.value();
            local_results[thread_id][tally_index] += 1;
          }
        }
        // Clear current chord length
        total_chord_length = 0.0;
      }
    }
  }

  // Combine results from all threads
  for (const auto& local_result : local_results) {
    for (size_t i = 0; i < result.chord_length.size(); ++i) {
      result.chord_length[i] += local_result[i];
    }
  }
  return result; // Return the result containing chord length statistics
}

bool ChordLengthStats::check_hit_boundary(const Particle& p) const
{
  //! \brief Check if the position is within the defined bounding box
  // Check if the position is within the defined boundaries
  if (p.r().x < lower_left_.x || p.r().y < lower_left_.y ||
      p.r().z < lower_left_.z || p.r().x > upper_right_.x ||
      p.r().y > upper_right_.y || p.r().z > upper_right_.z)
    return true; // If the position is not within bounds, return true
  if (p.surface() != SURFACE_NONE) {
    // If the particle is on a surface, check if the surface has a boundary
    // condition
    const auto& surf {model::surfaces[p.surface_index()].get()};
    if (surf->bc_) {
      return true; // If the surface has a boundary condition, return true
    }
  }
  return false; // If the position is within bounds, return false
}

// Determine if the current material ID and the previous material ID match
// both the base material and the stochastic medium material
bool ChordLengthStats::check_material_match(
  int32_t index1, int32_t index2) const
{
  auto current_material_id {model::materials[index1]->id()};
  auto last_material_id {model::materials[index2]->id()};

  if (last_material_id == matrix_domain_id_ &&
      current_material_id == stochastic_media_domain_id_) {
    return true; // Return true if the current material ID matches the
                 // stochastic media domain ID
  }

  return false; // If no match is found,return std::nullopt
}

std::optional<size_t> ChordLengthStats::find_tally_bin_index(
  double total_chord_length) const
{
  if (tally_bins_.empty()) {
    throw std::runtime_error("tally_bins_ is empty.");
  }

  // Find the interval index
  for (size_t i = 0; i < tally_bins_.size() - 1; ++i) {
    if (total_chord_length >= tally_bins_[i] &&
        total_chord_length < tally_bins_[i + 1]) {
      return i;
    }
  }

  // If total_chord_length matches the last bin boundary
  if (total_chord_length == tally_bins_.back()) {
    return tally_bins_.size() - 2;
  }

  return std::nullopt; //  If no match is found,return std::nullopt
}

void ChordLengthStats::to_hdf5(
  const std::string& filename, const Result& result) const
{
  //  Create an HDF5 file with write access
  hid_t file_id = file_open(filename, 'w');

  //  Write attributes to the HDF5 file
  write_attribute(file_id, "filetype", "chord_length_stats");
  write_attribute(file_id, "version", VERSION_VOLUME);
  write_attribute(file_id, "openmc_version", VERSION);
#ifdef GIT_SHA1
  write_attribute(file_id, "git_sha1", GIT_SHA1);
#endif

  //  Write attribute with current date and time
  write_attribute(file_id, "date_and_time", time_stamp());

  //  Write basic metadata
  write_attribute(file_id, "samples", n_samples_);
  write_attribute(file_id, "lower_left", lower_left_);
  write_attribute(file_id, "upper_right", upper_right_);

  //  Write attribute for domain type
  if (domain_type_ == TallyDomain::MATERIAL) {
    write_attribute(file_id, "domain_type", "material");
  } else {
    throw std::runtime_error(
      "Unsupported domain type for chord length statistics.");
  }
  write_dataset(file_id, "tally_bins", tally_bins_);
  //  Write results
  write_dataset(file_id, "chord_length", result.chord_length);

  //  Close the HDF5 file
  file_close(file_id);
}

void free_memory_chordl()
{
  // Free memory for chord length statistics
  model::chordl_stats.clear();
}
} // namespace openmc

int openmc_chord_length_stats()
{
  using namespace openmc;
  // Execute chord length statistics for each domain
  for (auto& stats : openmc::model::chordl_stats) {
    ChordLengthStats::Result result;
    result = stats.execute();
    // Write results to HDF5 file
    stats.to_hdf5(settings::path_output + "chord_length_stats.h5", result);
  }
  return 0;
}