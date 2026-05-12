#include "openmc/stochastic_media.h"

#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/geometry.h"
#include "openmc/material.h"
#include "openmc/surface.h"
#include "openmc/universe.h"
#include "openmc/xml_interface.h"

#include <cmath>
#include <fmt/core.h>

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace model {
std::unordered_map<int32_t, int32_t> stochastic_media_map;
std::vector<std::unique_ptr<StochasticMedia>> stochastic_media;
} // namespace model

//==============================================================================
// RNG utilities
//==============================================================================

uint64_t compute_stochastic_seed(const Position& r, const Direction& u,
  double grid_size, StochasticRNGRole role, uint64_t user_seed)
{
  uint64_t h = static_cast<uint64_t>(role);
  hash_combine(h, user_seed);
  hash_combine(h, quantize(r.x, grid_size));
  hash_combine(h, quantize(r.y, grid_size));
  hash_combine(h, quantize(r.z, grid_size));

  // Zero-direction defense: if direction is zero vector (abnormal path),
  // use default direction {0,0,1} to avoid all particles producing same seed.
  // Normal source particles have direction initialized before find_cell.
  if (u.x == 0.0 && u.y == 0.0 && u.z == 0.0) {
    hash_combine(h, 0.0);
    hash_combine(h, 0.0);
    hash_combine(h, 1.0);
  } else {
    hash_combine(h, u.x);
    hash_combine(h, u.y);
    hash_combine(h, u.z);
  }

  return h;
}

//==============================================================================
// StochasticMedia
//==============================================================================

double StochasticMedia::radius() const
{
  auto* sphere =
    dynamic_cast<SurfaceSphere*>(model::surfaces[boundary_surface_].get());
  if (!sphere) {
    fatal_error(fmt::format(
      "Boundary surface for stochastic media {} is not a sphere.", id_));
  }
  return sphere->radius_;
}

void StochasticMedia::from_xml(pugi::xml_node node)
{
  if (check_for_node(node, "id"))
    id_ = std::stoi(get_node_value(node, "id"));
  if (check_for_node(node, "name"))
    name_ = get_node_value(node, "name");
  if (check_for_node(node, "pack_fraction"))
    pf_ = std::stod(get_node_value(node, "pack_fraction"));
  // Read particle universe ID, particle cell ID, and matrix material ID
  // (converted to indices later by adjust_indices)
  if (check_for_node(node, "particle_universe"))
    particle_universe_ = std::stoi(get_node_value(node, "particle_universe"));
  if (check_for_node(node, "particle_cell"))
    particle_cell_ = std::stoi(get_node_value(node, "particle_cell"));
  if (check_for_node(node, "matrix_material"))
    matrix_fill_ = std::stoi(get_node_value(node, "matrix_material"));
  if (check_for_node(node, "seed"))
    seed_ = std::stoull(get_node_value(node, "seed"));
}

void StochasticMedia::adjust_indices()
{
  // particle_universe_ from universe ID to universe index
  auto search_univ = model::universe_map.find(particle_universe_);
  if (search_univ != model::universe_map.end()) {
    particle_universe_ = search_univ->second;
  } else {
    fatal_error(fmt::format("Could not find particle universe {} for "
                            "stochastic media {}.",
      particle_universe_, id_));
  }

  // particle_cell_ from cell ID to cell index
  auto search_cell = model::cell_map.find(particle_cell_);
  if (search_cell != model::cell_map.end()) {
    particle_cell_ = search_cell->second;
  } else {
    fatal_error(fmt::format("Could not find particle cell {} for "
                            "stochastic media {}.",
      particle_cell_, id_));
  }

  // matrix_fill_ from material ID to material index
  auto search_mat = model::material_map.find(matrix_fill_);
  if (search_mat != model::material_map.end()) {
    matrix_fill_ = search_mat->second;
  } else {
    fatal_error(fmt::format("Could not find matrix material {} for "
                            "stochastic media {}.",
      matrix_fill_, id_));
  }
}

void StochasticMedia::initialize()
{
  // Validate particle_universe contains exactly 2 cells
  auto& univ = *model::universes[particle_universe_];
  if (univ.cells_.size() != 2) {
    fatal_error(
      fmt::format("Stochastic media {}: particle_universe must "
                  "contain exactly 2 cells (particle cell + remainder cell), "
                  "but found {} cells.",
        id_, univ.cells_.size()));
  }

  // Validate particle_cell has sphere boundary
  Cell& pcell = *model::cells[particle_cell_];
  auto surf_tokens = pcell.surfaces();
  double max_radius = -1.0;
  int32_t boundary_surf_idx = C_NONE;
  for (auto token : surf_tokens) {
    int32_t surf_idx = std::abs(token) - 1;
    auto* sphere =
      dynamic_cast<SurfaceSphere*>(model::surfaces[surf_idx].get());
    if (sphere && sphere->radius_ > max_radius) {
      max_radius = sphere->radius_;
      boundary_surf_idx = surf_idx;
    }
  }
  if (boundary_surf_idx == C_NONE) {
    fatal_error(fmt::format("Particle cell {} for stochastic media {} "
                            "has no sphere surface as boundary.",
      pcell.id_, id_));
  }
  boundary_surface_ = boundary_surf_idx;

  // Detect remainder_cell (the cell that is not particle_cell)
  remainder_cell_ = C_NONE;
  for (auto cell_idx : univ.cells_) {
    if (cell_idx != particle_cell_) {
      remainder_cell_ = cell_idx;
      break;
    }
  }

  // Set derived fields
  grid_size_ = radius() / 100.0;
  double T = model::materials[matrix_fill_]->temperature();
  matrix_sqrtkT_ = std::sqrt(K_BOLTZMANN * T);
}

//==============================================================================
// CLSMedia
//==============================================================================

double CLSMedia::sample_matrix_chord(
  Position r, Direction u, ChordContext /*ctx*/) const
{
  // CLS: uniform exponential distribution for all contexts.
  double R = radius();
  double mean = 4.0 * R * (1.0 - pf_) / (3.0 * pf_);
  uint64_t seed =
    compute_stochastic_seed(r, u, grid_size_, StochasticRNGRole::CHORD_LENGTH, seed_);
  double xi = prn(&seed);
  return -mean * std::log(xi);
}

Position CLSMedia::compute_sphere_center(Position r, Direction u) const
{
  double R = radius();

  // Independent seeds to avoid correlation with chord_length
  uint64_t seed_xi = compute_stochastic_seed(
    r, u, grid_size_, StochasticRNGRole::SPHERE_CENTER_XI, seed_);
  uint64_t seed_phi = compute_stochastic_seed(
    r, u, grid_size_, StochasticRNGRole::SPHERE_CENTER_PHI, seed_);

  double xi = prn(&seed_xi);
  double phi = 2.0 * M_PI * prn(&seed_phi);

  // Perpendicular basis vectors (u.cross(ref) is non-zero when ref is not
  // parallel to u)
  Direction ref =
    (std::abs(u.x) < 0.9) ? Direction {1, 0, 0} : Direction {0, 1, 0};
  Direction n1 = u.cross(ref) / u.cross(ref).norm();
  Direction n2 = u.cross(n1) / u.cross(n1).norm();

  // The particle is at the entry point on the sphere surface after traveling
  // d_stoch. Compute sphere center C such that the particle is on the sphere
  // surface:
  //   d_perp^2 + d_parallel^2 = R^2
  // OpenMC's Surface::sense() handles on-surface classification via direction
  // tiebreaker.
  double d_perp = R * std::sqrt(xi);
  double d_parallel = R * std::sqrt(1.0 - xi);
  Direction perp =
    n1 * (d_perp * std::cos(phi)) + n2 * (d_perp * std::sin(phi));

  return r + d_parallel * u + perp;
}

Position StochasticMedia::compute_source_sphere_center(
  Position r, Direction u) const
{
  double R = radius();

  uint64_t seed_r =
    compute_stochastic_seed(r, u, grid_size_, StochasticRNGRole::SOURCE_RADIUS, seed_);
  uint64_t seed_theta =
    compute_stochastic_seed(r, u, grid_size_, StochasticRNGRole::SOURCE_POLAR, seed_);
  uint64_t seed_phi = compute_stochastic_seed(
    r, u, grid_size_, StochasticRNGRole::SOURCE_AZIMUTH, seed_);

  double xi_r = prn(&seed_r);
  double xi_theta = prn(&seed_theta);
  double xi_phi = prn(&seed_phi);

  // Uniform volume sampling: d = R * xi^(1/3) * scale
  double scale = 1.0 - 10.0 * FP_REL_PRECISION;
  double d = R * std::cbrt(xi_r) * scale;

  // Isotropic direction in global coordinates
  double cos_theta = 1.0 - 2.0 * xi_theta;
  double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
  double phi = 2.0 * M_PI * xi_phi;

  return r + d * Position{sin_theta * std::cos(phi),
                           sin_theta * std::sin(phi),
                           cos_theta};
}

} // namespace openmc
