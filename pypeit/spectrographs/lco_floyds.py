"""
Module for LCO FLOYDS

.. include:: ../include/links.rst
"""
import numpy as np

from astropy.coordinates import SkyCoord
from astropy import units

from pypeit import msgs
from pypeit import telescopes
from pypeit import io
from pypeit.core import parse
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container
from pypeit import dataPaths

from IPython import embed

class LCOFLOYDSSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle LCO/FLOYDS specific code
    """
    ndet = 1
    # This is north telescope, should implement the south 2m at some point
    telescope = telescopes.LCOFTNTelescopePar() 
    pypeline = 'Echelle'
    url = 'https://lco.global/observatory/instruments/floyds/'
    ech_fixed_format = True
    #header_name = # ???

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}

        self.meta['ra'] = dict(ext=0, card="RA")
        self.meta['dec'] = dict(ext=0, card='DEC')
        self.meta['target'] = dict(ext=0, card='OBJECT')
        self.meta['mjd'] = dict(ext=0, card='MJD-OBS')
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        self.meta['binning'] = dict(ext=0, card=None, default='1,1') # Check if we need to implement binning
        self.meta['instrument'] = dict(ext=0, card='INSTRUME')
        self.meta['slitwidth'] = dict(ext=0, card='APERWID')
        
        #self.meta['decker'] = dict(ext=0, card='SLIT') # Do I need this?

    def compound_meta(self, headarr, meta_key):
        """
        Methods to generate metadata requiring interpretation of the header
        data, instead of simply reading the value of a header card.

        Args:
            headarr (:obj:`list`):
                List of `astropy.io.fits.Header`_ objects.
            meta_key (:obj:`str`):
                Metadata keyword to construct.

        Returns:
            object: Metadata value read from the header(s).
        """
        # Don't know if we need this, placeholder for now
        # This comes from not_alfosc
        #if meta_key == 'binning':
        #    # PypeIt frame
        #    binspatial = headarr[0]['DETXBIN']
        #    binspec = headarr[0]['DETYBIN']
        #    return parse.binning2string(binspec, binspatial)
        #else :
        #    msgs.error("Not ready for this compound meta")
        return {}

    def configuration_keys(self):
        """
        Return the metadata keys that define a unique instrument
        configuration.

        This list is used by :class:`~pypeit.metadata.PypeItMetaData` to
        identify the unique configurations among the list of frames read
        for a given reduction.

        Returns:
            :obj:`list`: List of keywords of data pulled from file headers
            and used to construct the :class:`~pypeit.metadata.PypeItMetaData`
            object.
        """
        return ['slitwidth']

class LCOFLOYDSNorthSpectrograph(LCOFLOYDSSpectrograph):
    """
    Child to handle VLT/XSHOOTER specific code
    """

    name = 'lco_floyds_north'
    #camera = 'XShooter_VIS'
    supported = False
    #comment = 'See :doc:`xshooter`'

    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.

        Args:
            det (:obj:`int`):
                1-indexed detector number.
            hdu (`astropy.io.fits.HDUList`_, optional):
                The open fits file with the raw image of interest.  If not
                provided, frame-dependent parameters are set to a default.

        Returns:
            :class:`~pypeit.images.detector_container.DetectorContainer`:
            Object with the detector metadata.
        """
        # Binning
        # TODO: Could this be detector dependent??
        binning = '1,1' if hdu is None else self.get_meta_value(self.get_headarr(hdu), 'binning')

        # Detector 1
        # These are still X-shooter, get started!
        detector_dict = dict(
            binning         = binning,
            det              =1,
            dataext         = 0,
            specaxis        = 0,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.16, # average from order 17 and order 30, see manual
            darkcurr        = 0.0,  # e-/pixel/hour
            saturation      = 65535.,
            nonlinear       = 0.86,
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d(0.595), # FITS format is flipped: PrimaryHDU  (2106, 4000) w/respect to Python
            ronoise         = np.atleast_1d(3.1), # raw unbinned images are (4000,2106) (spec, spat)
            datasec=np.atleast_1d('[:,11:2058]'),  # pre and oscan are in the spatial direction
            oscansec=np.atleast_1d('[:,2059:2106]'),
        )
        return detector_container.DetectorContainer(**detector_dict)





















        
