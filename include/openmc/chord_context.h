#ifndef OPENMC_CHORD_CONTEXT_H
#define OPENMC_CHORD_CONTEXT_H

#include <cstdint>

namespace openmc {

// Chord sampling context for context-aware chord length sampling.
// CLS implementation ignores ctx (uniform exponential distribution).
enum class ChordContext : uint8_t {
  BOUNDARY_ENTER, // First entry into matrix from external boundary or cross-media
  PARTICLE_EXIT,  // Entry into matrix from particle exit
  SCATTER,        // Scatter in matrix changes direction
};

} // namespace openmc

#endif // OPENMC_CHORD_CONTEXT_H
