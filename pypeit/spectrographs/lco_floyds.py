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
    header_name = "en12"

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
        self.meta['idname'] = dict(ext=0, card='OBSTYPE')
        self.meta['mjd'] = dict(ext=0, card='MJD-OBS')
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        self.meta['binning'] = dict(ext=0, card=None, default='1,1') # Check if we need to implement binning
        self.meta['instrument'] = dict(ext=0, card='INSTRUME')
        self.meta['slitwid'] = dict(ext=0, card='APERWID')
        self.meta['decker'] = dict(ext=0, card='APERWID')
        self.meta['dispname'] = dict(ext=0, card='INSTRUME')
        
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
        # Single setup for grism? So slitwidth will be good enough.
        return ['dispname',"decker"]


    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.

        Args:
            ftype (:obj:`str`):
                Type of frame to check. Must be a valid frame type; see
                frame-type :ref:`frame_type_defs`.
            fitstbl (`astropy.table.Table`_):
                The table with the metadata for one or more frames to check.
            exprng (:obj:`list`, optional):
                Range in the allowed exposure time for a frame of type
                ``ftype``. See
                :func:`pypeit.core.framematch.check_frame_exptime`.

        Returns:
            `numpy.ndarray`_: Boolean array with the flags selecting the
            exposures in ``fitstbl`` that are ``ftype`` type frames.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)

        if ftype == 'science':
            return good_exp & (fitstbl['idname'] == 'SPECTRUM') 
        if ftype == 'standard':
            return good_exp & (fitstbl['target'] == 'STD,FLUX')
        # Seem to lack biases, come back here.
        if ftype == 'bias':
            return good_exp & (fitstbl['target'] == 'BIAS')
        if ftype == 'dark':
            return good_exp & (fitstbl['target'] == 'DARK')
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            return good_exp & ((fitstbl['idname'] == 'LAMPFLAT') | (fitstbl['target'] == 'LAMP,TRACE'))
        if ftype in ['arc', 'tilt']:
            return good_exp & ((fitstbl['idname'] == 'ARC') | (fitstbl['target'] == 'LAMP,WAVE'))

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

class LCOFLOYDSNorthSpectrograph(LCOFLOYDSSpectrograph):
    """
    Child to handle VLT/XSHOOTER specific code
    """

    name = 'lco_floyds_north'
    #camera = 'floyds'
    supported = False
    #comment = 'See :doc:`lco_floyds`' # update this when needed

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
        # Lot of info coming from headers, should probably re-write to read headers.
        detector_dict = dict(
            binning         = binning,
            det              =1,
            dataext         = 0, # is this used?
            specaxis        = 1,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.34, # https://lco.global/observatory/instruments/floyds/
            darkcurr        = 0.0,  # https://lco.global/observatory/instruments/floyds/
            saturation      = 38400., # This is in the file headers
            nonlinear       = 0.989, # Also from the file headers
            mincounts       = -1e10, # placeholder
            numamplifiers   = 1,
            gain            = np.atleast_1d(2.0), # https://lco.global/observatory/instruments/floyds/
            ronoise         = np.atleast_1d(3.3), # https://lco.global/observatory/instruments/floyds/
            datasec=np.atleast_1d('[1:2079,1:512]'),  # Taken from file headers
            oscansec=np.atleast_1d('[2049:2079,1:512]'), # Taken from file headers
        )
        return detector_container.DetectorContainer(**detector_dict)





















        
