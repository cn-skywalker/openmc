#ifndef OPENMC_CHORD_LENGTH_STATS_H
#define OPENMC_CHORD_LENGTH_STATS_H

#include <algorithm> // for find
#include <cstdint>
#include <string>
#include <vector>

#include "openmc/array.h"
#include "openmc/openmp_interface.h"
#include "openmc/particle.h"
#include "openmc/position.h"
#include "openmc/tallies/trigger.h"
#include "openmc/vector.h"

#include "pugixml.hpp"
#include "xtensor/xtensor.hpp"
#include <optional>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace openmc {
//==============================================================================
// Chord length statistic class
//==============================================================================

class ChordLengthStats {
public:
  // Aliases, types
  struct Result {
    vector<double> chord_length; //!< Frequency of chord length
  };                             // Results for a single domain

  // Constructors
  ChordLengthStats(pugi::xml_node node);

  ChordLengthStats() = default;

  // Methods

  //! \brief Stochastically determine the chord length statistics of a set of
  //! domains along with the average number densities of nuclides within the
  //! domain
  //!
  //! \return Vector of results for each user-specified domain
  Result execute() const;

  //! \brief Write chord length statistics results to HDF5 file
  //!
  //! \param[in] filename Path to HDF5 file to write
  //! \param[in] results Vector of results for each domain
  void to_hdf5(const std::string& filename, const Result& result) const;
  bool check_material_match(int32_t index1, int32_t index2) const;
  bool check_hit_boundary(const Particle& p) const;

  //! \brief Determine the index of the interval in tally_bins_ where
  //! total_chord_length falls
  //! \param[in] total_chord_length The chord length to be checked
  //! \return Index of the interval in tally_bins_
  std::optional<size_t> find_tally_bin_index(double total_chord_length) const;
  // Tally filter and map types
  enum class TallyDomain { MATERIAL };

  // Data members
  TallyDomain domain_type_;  //!< Type of domain (cell, material, etc.)
  size_t n_samples_;         //!< Number of samples to use
  Position lower_left_;      //!< Lower-left position of bounding box
  Position upper_right_;     //!< Upper-right position of bounding box
  int32_t matrix_domain_id_; //!< IDs of matrix domains
  int32_t stochastic_media_domain_id_; //!< IDs of stachastic media domains
  vector<double> tally_bins_;          //!< Bins for the chord length statistics
};

//==============================================================================
// Global variables
//==============================================================================

namespace model {
extern vector<ChordLengthStats> chordl_stats; //!< Chord length statistics
}

void free_memory_chordl();

} // namespace openmc

#endif // OPENMC_CHORD_LENGTH_STATS_H