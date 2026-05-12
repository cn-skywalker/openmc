#include "openmc/geometry.h"

#include <fmt/core.h>
#include <fmt/ostream.h>

#include "openmc/array.h"
#include "openmc/cell.h"
#include "openmc/constants.h"
#include "openmc/error.h"
#include "openmc/lattice.h"
#include "openmc/settings.h"
#include "openmc/simulation.h"
#include "openmc/stochastic_media.h"
#include "openmc/string_utils.h"
#include "openmc/surface.h"

namespace openmc {

//==============================================================================
// Global variables
//==============================================================================

namespace model {

int root_universe {-1};
int n_coord_levels;

vector<int64_t> overlap_check_count;

} // namespace model

//==============================================================================
// Non-member functions
//==============================================================================

bool check_cell_overlap(GeometryState& p, bool error)
{
  int n_coord = p.n_coord();

  // Loop through each coordinate level
  for (int j = 0; j < n_coord; j++) {
    Universe& univ = *model::universes[p.coord(j).universe()];

    // Loop through each cell on this level
    for (auto index_cell : univ.cells_) {
      Cell& c = *model::cells[index_cell];
      if (c.contains(p.coord(j).r(), p.coord(j).u(), p.surface())) {
        if (index_cell != p.coord(j).cell()) {
          if (error) {
            fatal_error(
              fmt::format("Overlapping cells detected: {}, {} on universe {}",
                c.id_, model::cells[p.coord(j).cell()]->id_, univ.id_));
          }
          return true;
        }
#pragma omp atomic
        ++model::overlap_check_count[index_cell];
      }
    }
  }

  return false;
}

//==============================================================================

int cell_instance_at_level(const GeometryState& p, int level)
{
  // throw error if the requested level is too deep for the geometry
  if (level > model::n_coord_levels) {
    fatal_error(fmt::format("Cell instance at level {} requested, but only {} "
                            "levels exist in the geometry.",
      level, p.n_coord()));
  }

  // determine the cell instance
  Cell& c {*model::cells[p.coord(level).cell()]};

  // quick exit if this cell doesn't have distribcell instances
  if (c.distribcell_index_ == C_NONE)
    return C_NONE;

  // compute the cell's instance
  int instance = 0;
  for (int i = 0; i < level; i++) {
    const auto& c_i {*model::cells[p.coord(i).cell()]};
    if (c_i.type_ == Fill::UNIVERSE) {
      instance += c_i.offset_[c.distribcell_index_];
    } else if (c_i.type_ == Fill::LATTICE) {
      instance += c_i.offset_[c.distribcell_index_];
      auto& lat {*model::lattices[p.coord(i + 1).lattice()]};
      const auto& i_xyz {p.coord(i + 1).lattice_index()};
      if (lat.are_valid_indices(i_xyz)) {
        instance += lat.offset(c.distribcell_index_, i_xyz);
      }
    } else if (c_i.type_ == Fill::STOCHASTIC_MEDIA) {
      instance += c_i.offset_[c.distribcell_index_];
    }
  }
  return instance;
}

//==============================================================================

bool find_cell_inner(
  GeometryState& p, const NeighborList* neighbor_list, bool verbose)
{
  // Find which cell of this universe the particle is in.  Use the neighbor list
  // to shorten the search if one was provided.
  bool found = false;
  int32_t i_cell = C_NONE;
  if (neighbor_list) {
    for (auto it = neighbor_list->cbegin(); it != neighbor_list->cend(); ++it) {
      i_cell = *it;

      // Make sure the search cell is in the same universe.
      int i_universe = p.lowest_coord().universe();
      if (model::cells[i_cell]->universe_ != i_universe)
        continue;

      // Check if this cell contains the particle.
      Position r {p.r_local()};
      Direction u {p.u_local()};
      auto surf = p.surface();
      if (model::cells[i_cell]->contains(r, u, surf)) {
        p.lowest_coord().cell() = i_cell;
        found = true;
        break;
      }
    }

    // If we're attempting a neighbor list search and fail, we
    // now know we should return false. This will trigger an
    // exhaustive search from neighbor_list_find_cell and make
    // the result from that be appended to the neighbor list.
    if (!found) {
      return found;
    }
  }

  // Check successively lower coordinate levels until finding material fill
  for (;; ++p.n_coord()) {
    // If we did not attempt to use neighbor lists, i_cell is still C_NONE.  In
    // that case, we should now do an exhaustive search to find the right value
    // of i_cell.
    //
    // Alternatively, neighbor list searches could have succeeded, but we found
    // that the fill of the neighbor cell was another universe. As such, in the
    // code below this conditional, we set i_cell back to C_NONE to indicate
    // that.
    if (i_cell == C_NONE) {
      int i_universe = p.lowest_coord().universe();
      const auto& univ {model::universes[i_universe]};
      found = univ->find_cell(p);
    }

    if (!found) {
      return found;
    }
    i_cell = p.lowest_coord().cell();

    // Announce the cell that the particle is entering.
    if (found && verbose) {
      auto msg = fmt::format("    Entering cell {}", model::cells[i_cell]->id_);
      write_message(msg, 1);
    }

    Cell& c {*model::cells[i_cell]};
    if (c.type_ == Fill::MATERIAL) {
      // Found a material cell which means this is the lowest coord level.

      p.cell_instance() = 0;
      // Find the distribcell instance number.
      if (c.distribcell_index_ >= 0) {
        p.cell_instance() = cell_instance_at_level(p, p.n_coord() - 1);
      }

      // Set the material, temperature and density multiplier.
      p.material_last() = p.material();
      p.material() = c.material(p.cell_instance());
      p.sqrtkT_last() = p.sqrtkT();
      p.sqrtkT() = c.sqrtkT(p.cell_instance());
      p.density_mult_last() = p.density_mult();
      p.density_mult() = c.density_mult(p.cell_instance());

      return true;

    } else if (c.type_ == Fill::UNIVERSE) {
      //========================================================================
      //! Found a lower universe, update this coord level then search the next.

      // Set the lower coordinate level universe.
      auto& coord {p.coord(p.n_coord())};
      coord.universe() = c.fill_;

      // Set the position and direction.
      coord.r() = p.r_local();
      coord.u() = p.u_local();

      // Apply translation.
      coord.r() -= c.translation_;

      // Apply rotation.
      if (!c.rotation_.empty()) {
        coord.rotate(c.rotation_);
      }

    } else if (c.type_ == Fill::LATTICE) {
      //========================================================================
      //! Found a lower lattice, update this coord level then search the next.

      Lattice& lat {*model::lattices[c.fill_]};

      // Set the position and direction.
      auto& coord {p.coord(p.n_coord())};
      coord.r() = p.r_local();
      coord.u() = p.u_local();

      // Apply translation.
      coord.r() -= c.translation_;

      // Apply rotation.
      if (!c.rotation_.empty()) {
        coord.rotate(c.rotation_);
      }

      // Determine lattice indices.
      auto& i_xyz {coord.lattice_index()};
      lat.get_indices(coord.r(), coord.u(), i_xyz);

      // Get local position in appropriate lattice cell
      coord.r() = lat.get_local_position(coord.r(), i_xyz);

      // Set lattice indices.
      coord.lattice() = c.fill_;

      // Set the lower coordinate level universe.
      if (lat.are_valid_indices(i_xyz)) {
        coord.universe() = lat[i_xyz];
      } else {
        if (lat.outer_ != NO_OUTER_UNIVERSE) {
          coord.universe() = lat.outer_;
        } else {
          p.mark_as_lost(fmt::format(
            "Particle {} left lattice {}, but it has no outer definition.",
            p.id(), lat.id_));
        }
      }
    } else if (c.type_ == Fill::STOCHASTIC_MEDIA) {
      //========================================================================
      //! Stochastic media cell — handle φ allocation and set phase material.
      //!
      //! find_cell_inner always returns the matrix phase by default.
      //! For source particles (cell_last == C_NONE), a position-based RNG
      //! decides whether to enter the particle phase (φ probability).
      //! If entering particle phase, set up the next coord level and let the
      //! for loop continue into the particle universe (same as Fill::UNIVERSE).
      auto& media = *model::stochastic_media[c.fill_];

      p.cell_instance() = 0;
      if (c.distribcell_index_ >= 0) {
        p.cell_instance() = cell_instance_at_level(p, p.n_coord() - 1);
      }

      // Source particle φ allocation (position-based RNG, not particle RNG)
      bool entering_particle = false;
      if (p.cell_last(p.n_coord() - 1) == C_NONE) {
        uint64_t seed = compute_stochastic_seed(p.r_local(), p.u_local(),
          media.grid_size_, StochasticRNGRole::PHASE_ALLOCATION, media.seed_);
        if (prn(&seed) < media.pf_) {
          entering_particle = true;
          // Set up next coord level (same pattern as Fill::UNIVERSE).
          // The for loop will ++p.n_coord() and exhaustive-search the particle
          // universe.
          Position C =
            media.compute_source_sphere_center(p.r_local(), p.u_local());
          auto& coord {p.coord(p.n_coord())};
          coord.universe() = media.particle_universe_;
          coord.r() = p.r_local() - C;
          coord.u() = p.u_local();
          // lattice_ is already C_NONE from caller's coord reset
        }
      }

      if (entering_particle) {
        // Don't return — fall through to i_cell=C_NONE; found=false
        // so the loop continues into the particle universe.
      } else {
        // Matrix phase (default or (1-φ) probability)
        p.material_last() = p.material();
        p.material() = media.matrix_fill_;
        p.sqrtkT_last() = p.sqrtkT();
        p.sqrtkT() = media.matrix_sqrtkT_;
        p.density_mult_last() = p.density_mult();
        p.density_mult() = 1.0;
        return true;
      }
    }
    i_cell = C_NONE; // trip non-neighbor cell search at next iteration
    found = false;
  }

  return found;
}

//==============================================================================

bool neighbor_list_find_cell(GeometryState& p, bool verbose)
{

  // Reset all the deeper coordinate levels.
  for (int i = p.n_coord(); i < model::n_coord_levels; i++) {
    p.coord(i).reset();
  }

  // Get the cell this particle was in previously.
  auto coord_lvl = p.n_coord() - 1;
  auto i_cell = p.coord(coord_lvl).cell();
  Cell& c {*model::cells[i_cell]};

  // Search for the particle in that cell's neighbor list.  Return if we
  // found the particle.
  bool found = find_cell_inner(p, &c.neighbors_, verbose);
  if (found)
    return found;

  // The particle could not be found in the neighbor list.  Try searching all
  // cells in this universe, and update the neighbor list if we find a new
  // neighboring cell.
  found = find_cell_inner(p, nullptr, verbose);
  if (found)
    c.neighbors_.push_back(p.coord(coord_lvl).cell());
  return found;
}

bool exhaustive_find_cell(GeometryState& p, bool verbose)
{
  int i_universe = p.lowest_coord().universe();
  if (i_universe == C_NONE) {
    p.coord(0).universe() = model::root_universe;
    p.n_coord() = 1;
    i_universe = model::root_universe;
  }
  // Reset all the deeper coordinate levels.
  for (int i = p.n_coord(); i < model::n_coord_levels; i++) {
    p.coord(i).reset();
  }
  return find_cell_inner(p, nullptr, verbose);
}

//==============================================================================

void cross_stochastic_boundary(
  GeometryState& p, const BoundaryInfo& info, bool verbose)
{
  // Direction detection: n_coord_last saved in event_cross_surface.
  // event_cross_surface sets n_coord = coord_level.
  // So n_coord_last <= coord_level → entering; n_coord_last > coord_level →
  // exiting
  bool entering = (p.n_coord_last() <= info.coord_level());

  // Resolve stochastic media reference (used by both entering and exiting
  // paths)
  int stoch_level = info.coord_level() - 1;
  Cell& c_stoch {*model::cells[p.coord(stoch_level).cell()]};
  auto& media = *model::stochastic_media[c_stoch.fill_];

  if (entering) {
    // ================================================================
    // Matrix→Particle: set offset coordinates, delegate to find_cell
    // ================================================================
    if (verbose) {
      write_message(
        fmt::format("    Entering stochastic media {} particle. r=({}, {}, {})",
          media.id_, p.r().x, p.r().y, p.r().z),
        1);
    }

    Position C = media.compute_sphere_center(p.r_local(), p.u_local());

    auto& coord {p.coord(p.n_coord())};
    coord.universe() = media.particle_universe_;
    coord.r() = p.r_local() - C;
    coord.u() = p.u_local();
    coord.lattice() = C_NONE;
    ++p.n_coord();

    if (!exhaustive_find_cell(p, verbose)) {
      p.mark_as_lost(
        fmt::format("Particle {} could not be located inside particle of "
                    "stochastic media {}.",
          p.id(), media.id_));
    }

  } else {
    if (verbose) {
      write_message(
        fmt::format("    Exiting stochastic media {} particle. r=({}, {}, {})",
          media.id_, p.r().x, p.r().y, p.r().z),
        1);
    }
    // ================================================================
    // Particle→Matrix: delegate to find_cell
    // ================================================================
    for (int j = p.n_coord(); j < model::n_coord_levels; ++j) {
      p.coord(j).reset();
    }

    if (!exhaustive_find_cell(p, verbose)) {
      p.mark_as_lost(fmt::format("Particle {} could not be located after "
                                 "exiting stochastic media particle.",
        p.id()));
    }
  }
}

//==============================================================================

void cross_lattice(GeometryState& p, const BoundaryInfo& boundary, bool verbose)
{
  auto& coord {p.lowest_coord()};
  auto& lat {*model::lattices[coord.lattice()]};

  if (verbose) {
    write_message(
      fmt::format("    Crossing lattice {}. Current position ({},{},{}). r={}",
        lat.id_, coord.lattice_index()[0], coord.lattice_index()[1],
        coord.lattice_index()[2], p.r()),
      1);
  }

  // Set the lattice indices.
  coord.lattice_index()[0] += boundary.lattice_translation()[0];
  coord.lattice_index()[1] += boundary.lattice_translation()[1];
  coord.lattice_index()[2] += boundary.lattice_translation()[2];

  // Set the new coordinate position.
  const auto& upper_coord {p.coord(p.n_coord() - 2)};
  const auto& cell {model::cells[upper_coord.cell()]};
  Position r = upper_coord.r();
  r -= cell->translation_;
  if (!cell->rotation_.empty()) {
    r = r.rotate(cell->rotation_);
  }
  p.r_local() = lat.get_local_position(r, coord.lattice_index());

  if (!lat.are_valid_indices(coord.lattice_index())) {
    // The particle is outside the lattice.  Search for it from the base coords.
    p.n_coord() = 1;
    bool found = exhaustive_find_cell(p);

    if (!found) {
      p.mark_as_lost(fmt::format("Particle {} could not be located after "
                                 "crossing a boundary of lattice {}",
        p.id(), lat.id_));
    }

  } else {
    // Find cell in next lattice element.
    p.lowest_coord().universe() = lat[coord.lattice_index()];
    bool found = exhaustive_find_cell(p);

    if (!found) {
      // A particle crossing the corner of a lattice tile may not be found.  In
      // this case, search for it from the base coords.
      p.n_coord() = 1;
      bool found = exhaustive_find_cell(p);
      if (!found) {
        p.mark_as_lost(fmt::format("Particle {} could not be located after "
                                   "crossing a boundary of lattice {}",
          p.id(), lat.id_));
      }
    }
  }
}

//==============================================================================

BoundaryInfo distance_to_boundary(GeometryState& p)
{
  BoundaryInfo info;
  int32_t level_surf_cross;
  array<int, 3> level_lat_trans {};

  // Loop over each coordinate level.
  for (int i = 0; i < p.n_coord(); i++) {
    const auto& coord {p.coord(i)};
    const Position& r {coord.r()};
    const Direction& u {coord.u()};
    Cell& c {*model::cells[coord.cell()]};

    // 1. Cell surface distance (existing logic)
    auto surface_distance = c.distance(r, u, p.surface(), &p);
    double d_surf = surface_distance.first;
    level_surf_cross = surface_distance.second;

    // 2. Lattice distance (existing logic unchanged)
    double d_lat {INFINITY};
    if (coord.lattice() != C_NONE) {
      auto& lat {*model::lattices[coord.lattice()]};
      std::pair<double, array<int, 3>> lattice_distance;
      switch (lat.type_) {
      case LatticeType::rect:
        lattice_distance = lat.distance(r, u, coord.lattice_index());
        break;
      case LatticeType::hex:
        auto& cell_above {model::cells[p.coord(i - 1).cell()]};
        Position r_hex {p.coord(i - 1).r()};
        r_hex -= cell_above->translation_;
        if (coord.rotated()) {
          r_hex = r_hex.rotate(cell_above->rotation_);
        }
        r_hex.z = coord.r().z;
        lattice_distance = lat.distance(r_hex, u, coord.lattice_index());
        break;
      }
      d_lat = lattice_distance.first;
      level_lat_trans = lattice_distance.second;

      if (d_lat < 0) {
        p.mark_as_lost(fmt::format("Particle {} had a negative distance "
                                   "to a lattice boundary.",
          p.id()));
      }
    }

    // 3. Stochastic media chord length distance
    // Only sample when in matrix phase (not in particle phase)
    double d_stoch {INFINITY};
    if (c.type_ == Fill::STOCHASTIC_MEDIA) {
      auto& media = *model::stochastic_media[c.fill_];
      // Phase check: not in remainder cell means in particle (or nested
      // universe inside particle)
      bool in_particle =
        (p.n_coord() > i + 1 && p.coord(i + 1).cell() != media.remainder_cell_);
      if (!in_particle) {
        auto ctx = p.chord_context();
        d_stoch = media.sample_matrix_chord(r, u, ctx);
        p.chord_context() = ChordContext::SCATTER;
        if (settings::verbosity >= 10) {
          const char* ctx_name =
            (ctx == ChordContext::BOUNDARY_ENTER) ? "BOUNDARY_ENTER" :
            (ctx == ChordContext::PARTICLE_EXIT)  ? "PARTICLE_EXIT" :
                                                   "SCATTER";
          write_message(fmt::format(
            "    distance_to_boundary: {} (CLS), "
            "d_stoch={:.6f}, r=({}, {}, {})",
            ctx_name, d_stoch, r.x, r.y, r.z),
            1);
        }
      }
    }

    // 4. Three-way competition
    double& d = info.distance();
    if (d_stoch < d_surf - FP_COINCIDENT && d_stoch < d_lat - FP_COINCIDENT) {
      // d_stoch wins -> matrix->particle crossing
      if (d == INFINITY || (d - d_stoch) / d >= FP_REL_PRECISION) {
        d = d_stoch;
        info.surface() = SURFACE_NONE;
        info.lattice_translation() = {0, 0, 0};
        info.coord_level() = i + 1;
        info.stochastic_boundary() = true;
      }
    } else if (d_surf < d_lat - FP_COINCIDENT) {
      // d_surf wins
      if (d == INFINITY || (d - d_surf) / d >= FP_REL_PRECISION) {
        d = d_surf;

        // Boundary sphere detection (particle phase exit -> matrix)
        info.stochastic_boundary() = false;
        if (i > 0) {
          Cell& c_above {*model::cells[p.coord(i - 1).cell()]};
          if (c_above.type_ == Fill::STOCHASTIC_MEDIA) {
            auto& media_above = *model::stochastic_media[c_above.fill_];
            if (std::abs(level_surf_cross) - 1 ==
                  media_above.boundary_surface_ &&
                level_surf_cross > 0) {
              info.stochastic_boundary() = true;
              info.coord_level() = i;
              info.surface() = SURFACE_NONE;
              info.lattice_translation() = {0, 0, 0};
            }
          }
        }

        if (!info.stochastic_boundary()) {
          // Normal surface crossing (same as original code)
          if (c.is_simple() || d == INFTY) {
            info.surface() = level_surf_cross;
          } else {
            Position r_hit = r + d_surf * u;
            Surface& surf {*model::surfaces[std::abs(level_surf_cross) - 1]};
            Direction norm = surf.normal(r_hit);
            info.surface() = (u.dot(norm) > 0) ? std::abs(level_surf_cross)
                                               : -std::abs(level_surf_cross);
          }
          info.lattice_translation() = {0, 0, 0};
          info.coord_level() = i + 1;
        }
      }
    } else {
      // d_lat wins (existing logic unchanged)
      if (d == INFINITY || (d - d_lat) / d >= FP_REL_PRECISION) {
        d = d_lat;
        info.surface() = SURFACE_NONE;
        info.lattice_translation() = level_lat_trans;
        info.coord_level() = i + 1;
        info.stochastic_boundary() = false;
      }
    }
  }
  return info;
}

//==============================================================================
// C API
//==============================================================================

extern "C" int openmc_find_cell(
  const double* xyz, int32_t* index, int32_t* instance)
{
  GeometryState geom_state;

  geom_state.r() = Position {xyz};
  geom_state.u() = {0.0, 0.0, 1.0};

  if (!exhaustive_find_cell(geom_state)) {
    set_errmsg(
      fmt::format("Could not find cell at position {}.", geom_state.r()));
    return OPENMC_E_GEOMETRY;
  }

  *index = geom_state.lowest_coord().cell();
  *instance = geom_state.cell_instance();
  return 0;
}

extern "C" int openmc_global_bounding_box(double* llc, double* urc)
{
  auto bbox = model::universes.at(model::root_universe)->bounding_box();

  // set lower left corner values
  llc[0] = bbox.xmin;
  llc[1] = bbox.ymin;
  llc[2] = bbox.zmin;

  // set upper right corner values
  urc[0] = bbox.xmax;
  urc[1] = bbox.ymax;
  urc[2] = bbox.zmax;

  return 0;
}

} // namespace openmc
