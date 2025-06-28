from collections.abc import Iterable 
from numbers import Real, Integral
import warnings

import h5py
import lxml.etree as ET

import openmc.checkvalue as cv
from openmc._xml import get_text

_VERSION_CHORD_LENGTH = 1


class ChordLengthStats:
    """Stochastic chord length calculation specifications and results.

    Parameters
    ----------
    matrix_domain_id : int
        ID of the matrix material domain.
    stochastic_media_domain_id : int
        ID of the stochastic media material domain.
    samples : int
        Number of samples used to generate chord length statistics.
    lower_left : Iterable of float
        Lower-left coordinates of bounding box used to sample points.
    upper_right : Iterable of float
        Upper-right coordinates of bounding box used to sample points.
    tally_bins : Iterable of float
        Bin edges for tallying chord lengths.

    Attributes
    ----------
    matrix_domain_id : int
        ID of the matrix material domain.
    stochastic_media_domain_id : int
        ID of the stochastic media material domain.
    samples : int
        Number of samples used to generate chord length statistics.
    lower_left : Iterable of float
        Lower-left coordinates of bounding box used to sample points.
    upper_right : Iterable of float
        Upper-right coordinates of bounding box used to sample points.
    tally_bins : Iterable of float
        Bin edges for tallying chord lengths.
    chord_length : dict
        Dictionary mapping bin indices to chord length frequencies.
    """

    def __init__(self, matrix_domain_id, stochastic_media_domain_id, samples,tally_bins,
                 lower_left=None, upper_right=None):
        self.domain_type = 'material'
        self.matrix_domain_id = matrix_domain_id
        self.stochastic_media_domain_id = stochastic_media_domain_id
        self.samples = samples
        self.lower_left = lower_left
        self.upper_right = upper_right
        self.tally_bins = tally_bins
        self.chord_length = {}

    @property
    def samples(self):
        return self._samples

    @samples.setter
    def samples(self, samples):
        cv.check_type('number of samples', samples, Integral)
        cv.check_greater_than('number of samples', samples, 0)
        self._samples = samples

    @property
    def lower_left(self):
        return self._lower_left

    @lower_left.setter
    def lower_left(self, lower_left):
        name = 'lower-left bounding box coordinates'
        cv.check_type(name, lower_left, Iterable, Real)
        cv.check_length(name, lower_left, 3)
        self._lower_left = lower_left

    @property
    def upper_right(self):
        return self._upper_right

    @upper_right.setter
    def upper_right(self, upper_right):
        name = 'upper-right bounding box coordinates'
        cv.check_type(name, upper_right, Iterable, Real)
        cv.check_length(name, upper_right, 3)
        self._upper_right = upper_right

    @property
    def tally_bins(self):
        return self._tally_bins

    @tally_bins.setter
    def tally_bins(self, tally_bins):
        cv.check_type('tally bins', tally_bins, Iterable, Real)
        cv.check_greater_than('number of tally bins', len(tally_bins), 1)
        self._tally_bins = tally_bins

    def to_xml_element(self):
        """Return XML representation of the chord length calculation.

        Returns
        -------
        element : lxml.etree._Element
            XML element containing chord length calculation data.
        """
        element = ET.Element("chord_length_stats")
        ET.SubElement(element, "domain_type").text = self.domain_type
        ET.SubElement(element, "matrix_domain_id").text = str(self.matrix_domain_id)
        ET.SubElement(element, "stochastic_media_domain_id").text = str(self.stochastic_media_domain_id)
        ET.SubElement(element, "samples").text = str(self.samples)
        ET.SubElement(element, "lower_left").text = ' '.join(map(str, self.lower_left))
        ET.SubElement(element, "upper_right").text = ' '.join(map(str, self.upper_right))
        ET.SubElement(element, "tally_bins").text = ' '.join(map(str, self.tally_bins))
        return element

    @classmethod
    def from_xml_element(cls, elem):
        """Generate chord length calculation object from an XML element.

        Parameters
        ----------
        elem : lxml.etree._Element
            XML element.

        Returns
        -------
        openmc.ChordLengthStats
            Chord length calculation object.
        """
        matrix_domain_id = int(get_text(elem, "matrix_domain_id"))
        stochastic_media_domain_id = int(get_text(elem, "stochastic_media_domain_id"))
        samples = int(get_text(elem, "samples"))
        lower_left = list(map(float, get_text(elem, "lower_left").split()))
        upper_right = list(map(float, get_text(elem, "upper_right").split()))
        tally_bins = list(map(float, get_text(elem, "tally_bins").split()))
        return cls(matrix_domain_id, stochastic_media_domain_id, samples, lower_left, upper_right, tally_bins)

    def to_hdf5(self, filename):
        """Write chord length statistics to an HDF5 file.

        Parameters
        ----------
        filename : str
            Path to the output HDF5 file.
        """
        with h5py.File(filename, 'w') as f:
            f.attrs['filetype'] = 'chord_length_stats'
            f.attrs['version'] = _VERSION_CHORD_LENGTH
            f.attrs['matrix_domain_id'] = self.matrix_domain_id
            f.attrs['stochastic_media_domain_id'] = self.stochastic_media_domain_id
            f.attrs['samples'] = self.samples
            f.attrs['lower_left'] = self.lower_left
            f.attrs['upper_right'] = self.upper_right
            f.create_dataset('tally_bins', data=self.tally_bins)
            f.create_dataset('chord_length', data=list(self.chord_length.items()))

    @classmethod
    def from_hdf5(cls, filename):
        """Load chord length statistics from an HDF5 file.

        Parameters
        ----------
        filename : str
            Path to the HDF5 file.

        Returns
        -------
        openmc.ChordLengthStats
            Chord length statistics object.
        """
        with h5py.File(filename, 'r') as f:
            cv.check_filetype_version(f, 'chord_length_stats', _VERSION_CHORD_LENGTH)
            matrix_domain_id = f.attrs['matrix_domain_id']
            stochastic_media_domain_id = f.attrs['stochastic_media_domain_id']
            samples = f.attrs['samples']
            lower_left = f.attrs['lower_left']
            upper_right = f.attrs['upper_right']
            tally_bins = f['tally_bins'][:]
            chord_length = dict(f['chord_length'][:])
        stats = cls(matrix_domain_id, stochastic_media_domain_id, samples, lower_left, upper_right, tally_bins)
        stats.chord_length = chord_length
        return stats