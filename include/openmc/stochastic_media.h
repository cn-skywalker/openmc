#ifndef OPENMC_STOCHASTIC_MEDIA_H
#define OPENMC_STOCHASTIC_MEDIA_H

#include "openmc/chord_context.h"
#include "openmc/constants.h"
#include "openmc/position.h"
#include "openmc/random_lcg.h"
#include "pugixml.hpp"
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>
#include <memory>
#include <cstring>
#include <cmath>

namespace openmc {

// Forward declarations
class Surface;

//==============================================================================
// Stochastic media RNG utilities
//==============================================================================

// Different operations use different salts to avoid RNG correlation
enum class StochasticRNGRole : uint64_t {
  CHORD_LENGTH      = 0x9E3779B97F4A7C15ULL,
  SPHERE_CENTER_XI  = 0x517CC1B727220A95ULL,
  SPHERE_CENTER_PHI = 0x6C62272E07BB0142ULL,
  PHASE_ALLOCATION  = 0x3F2E1A0B7C5D9E8FULL,
  SOURCE_RADIUS     = 0x4A7B3C9D1E2F5068ULL,
  SOURCE_POLAR      = 0x8D6E5F4A3B2C1D0EULL,
  SOURCE_AZIMUTH    = 0x2F1E3D4C5B6A7980ULL,
};

// Quantization grid: eliminate floating-point path dependency
inline uint64_t quantize(double v, double grid_size)
{
  return static_cast<uint64_t>(std::floor(v / grid_size));
}

// hash_combine (boost-style)
inline void hash_combine(uint64_t& h, uint64_t v)
{
  h ^= v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
}

inline void hash_combine(uint64_t& h, double v)
{
  uint64_t tmp;
  std::memcpy(&tmp, &v, sizeof(double));
  hash_combine(h, tmp);
}

// Salted seed computation (with zero-direction defense)
uint64_t compute_stochastic_seed(const Position& r, const Direction& u,
  double grid_size, StochasticRNGRole role, uint64_t user_seed = 0);

//==============================================================================
// StochasticMedia abstract base class
//==============================================================================

class StochasticMedia {
public:
  // Basic info
  int32_t id_ {C_NONE};
  std::string name_;
  double pf_ {0.0};  // Packing fraction

  // Indices (converted from IDs during initialization)
  int32_t particle_cell_ {C_NONE};     // Particle cell index
  int32_t remainder_cell_ {C_NONE};    // Remainder cell index (auto-detected)
  int32_t boundary_surface_ {C_NONE};  // Boundary sphere surface index
  int32_t matrix_fill_ {C_NONE};       // Matrix material index
  double matrix_sqrtkT_ {0.0};         // Matrix temperature sqrt(kT)
  double grid_size_ {0.0};             // Quantization grid size
  int32_t particle_universe_ {C_NONE}; // Particle universe index (user-provided)
  uint64_t seed_ {0};                   // User-configurable RNG seed for different realizations

  // Get particle radius (from boundary_surface_ sphere)
  double radius() const;

  // Chord sampling interface
  virtual double sample_matrix_chord(
    Position r, Direction u, ChordContext ctx) const = 0;
  virtual Position compute_sphere_center(Position r, Direction u) const = 0;

  // Source particle interior sampling: uniform volume sampling inside sphere.
  // u is used only for seed diversification and orthogonal basis, not physical direction.
  Position compute_source_sphere_center(Position r, Direction u) const;

  // XML reading
  virtual void from_xml(pugi::xml_node node);

  // ID-to-index conversion (called from adjust_indices)
  void adjust_indices();

  // Validation and derived-field computation (called from finalize_geometry)
  virtual void initialize();

  virtual ~StochasticMedia() = default;
};

//==============================================================================
// CLSMedia
//==============================================================================

class CLSMedia : public StochasticMedia {
public:
  double sample_matrix_chord(
    Position r, Direction u, ChordContext ctx) const override;
  Position compute_sphere_center(Position r, Direction u) const override;
};

//==============================================================================
// model namespace declarations
//==============================================================================

namespace model {
extern std::unordered_map<int32_t, int32_t> stochastic_media_map;
extern std::vector<std::unique_ptr<StochasticMedia>> stochastic_media;
}

} // namespace openmc

#endif // OPENMC_STOCHASTIC_MEDIA_H
