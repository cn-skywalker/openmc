#ifndef OPENMC_CHORD_LENGTH_STATS_H
#define OPENMC_CHORD_LENGTH_STATS_H

#include <algorithm> // for find
#include <cstdint>
#include <string>
#include <vector>

#include "openmc/array.h"
#include "openmc/openmp_interface.h"
#include "openmc/position.h"
#include "openmc/tallies/trigger.h"
#include "openmc/vector.h"

#include "pugixml.hpp"
#include "xtensor/xtensor.hpp"
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
    vector<array<double, 2>>
      chord_length; //!< Mean/standard deviation of chord length
  };                // Results for a single domain

  // Constructors
  ChordLengthStats(pugi::xml_node node);

  ChordLengthStats() = default;

  // Methods

  //! \brief Stochastically determine the chord length statistics of a set of
  //! domains along with the average number densities of nuclides within the
  //! domain
  //!
  //! \return Vector of results for each user-specified domain
  vector<Result> execute() const;

  //! \brief Write chord length statistics results to HDF5 file
  //!
  //! \param[in] filename Path to HDF5 file to write
  //! \param[in] results Vector of results for each domain
  void to_hdf5(
    const std::string& filename, const vector<Result>& results) const;

  // Tally filter and map types
  enum class TallyDomain { MATERIAL };

  // Data members
  TallyDomain domain_type_;       //!< Type of domain (cell, material, etc.)
  size_t n_samples_;              //!< Number of samples to use
  Position lower_left_;           //!< Lower-left position of bounding box
  Position upper_right_;          //!< Upper-right position of bounding box
  vector<int> matrix_domain_ids_; //!< IDs of matrix domains
  vector<int> stochastic_media_domain_ids_; //!< IDs of stachastic media domains
  vector<double> tally_bins; //!< Bins for the chord length statistics
};

//==============================================================================
// Global variables
//==============================================================================

namespace model {
extern vector<ChordLengthStats> chordl_stats; //!< Chord length statistics
}

void free_memory_ChordL();

} // namespace openmc

#endif // OPENMC_CHORD_LENGTH_STATS_H