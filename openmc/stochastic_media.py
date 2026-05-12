from numbers import Real

import lxml.etree as ET

import openmc
import openmc.checkvalue as cv


class StochasticMedia:
    """Abstract base class for stochastic media transport methods.

    Parameters
    ----------
    stochastic_media_id : int, optional
        Unique ID for the stochastic media. If not specified, an identifier
        will automatically be assigned.
    name : str, optional
        Name of the stochastic media.

    Attributes
    ----------
    id : int or None
        Unique identifier for the stochastic media.
    name : str
        Name of the stochastic media.
    pack_fraction : float
        Volume fraction occupied by particles.
    particle_universe : openmc.Universe
        Universe containing the particle cell and remainder cell.
    particle_cell : openmc.Cell
        Cell defining the particle geometry (must belong to particle_universe).
    matrix_material : openmc.Material
        Material filling the matrix (inter-particle) region.
    """

    def __init__(self, stochastic_media_id=None, name=''):
        self.id = stochastic_media_id
        self.name = name
        self._pf = None
        self._particle_universe = None
        self._particle_cell = None
        self._matrix_material = None
        self._seed = None

    @property
    def pack_fraction(self):
        return self._pf

    @pack_fraction.setter
    def pack_fraction(self, val):
        cv.check_type('pack fraction', val, Real)
        cv.check_greater_than('pack fraction', val, 0.0, equality=True)
        if val > 1.0:
            raise ValueError('Pack fraction must be <= 1.0')
        self._pf = val

    @property
    def particle_universe(self):
        return self._particle_universe

    @particle_universe.setter
    def particle_universe(self, univ):
        cv.check_type('particle universe', univ, openmc.Universe)
        self._particle_universe = univ

    @property
    def particle_cell(self):
        return self._particle_cell

    @particle_cell.setter
    def particle_cell(self, cell):
        cv.check_type('particle cell', cell, openmc.Cell)
        self._particle_cell = cell

    @property
    def matrix_material(self):
        return self._matrix_material

    @matrix_material.setter
    def matrix_material(self, mat):
        cv.check_type('matrix material', mat, openmc.Material)
        self._matrix_material = mat

    @property
    def seed(self):
        return self._seed

    @seed.setter
    def seed(self, val):
        cv.check_type('seed', val, int)
        cv.check_greater_than('seed', val, 0, equality=True)
        self._seed = val


class CLSMedia(StochasticMedia):
    """Chord Length Sampling stochastic media.

    Parameters
    ----------
    stochastic_media_id : int, optional
        Unique ID for the stochastic media. If not specified, an identifier
        will automatically be assigned.
    name : str, optional
        Name of the stochastic media.

    Attributes
    ----------
    id : int or None
        Unique identifier for the stochastic media.
    name : str
        Name of the stochastic media.
    pack_fraction : float
        Volume fraction occupied by particles.
    particle_universe : openmc.Universe
        Universe containing the particle cell and remainder cell.
    particle_cell : openmc.Cell
        Cell defining the particle geometry (must belong to particle_universe).
    matrix_material : openmc.Material
        Material filling the matrix (inter-particle) region.
    """

    def __init__(self, stochastic_media_id=None, name=''):
        super().__init__(stochastic_media_id, name)

    def to_xml_element(self, memo=None):
        """Return XML representation of the CLS media.

        Returns
        -------
        lxml.etree._Element
            XML element containing CLS media data.
        """
        element = ET.Element('cls_media')
        if self.id is not None:
            element.set('id', str(self.id))
        if self.name:
            element.set('name', self.name)
        if self._pf is not None:
            element.set('pack_fraction', str(self._pf))
        if self._particle_universe is not None:
            element.set('particle_universe', str(self._particle_universe.id))
        if self._particle_cell is not None:
            element.set('particle_cell', str(self._particle_cell.id))
        if self._matrix_material is not None:
            element.set('matrix_material', str(self._matrix_material.id))
        if self._seed is not None:
            element.set('seed', str(self._seed))
        return element

    def export_particle_universe(self, element, memo=None):
        """Export all cells in the particle universe to geometry XML.

        The particle universe is not part of the root universe hierarchy,
        so its cells must be explicitly exported.

        Parameters
        ----------
        element : lxml.etree._Element
            The geometry XML element to append to.
        memo : set, optional
            Shared memo set for surface deduplication. If None, a new set
            is created. Must be shared with the main geometry export to
            avoid duplicate surface XML elements.
        """
        if self._particle_universe is not None:
            if memo is None:
                memo = set()
            for cell in self._particle_universe.cells.values():
                cell_elem = cell.create_xml_subelement(element, memo)
                cell_elem.set("universe", str(self._particle_universe.id))
                element.append(cell_elem)
