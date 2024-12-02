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
            datasec= np.atleast_1d('[1:512,1:2048]'),  # Taken from file headers
            oscansec= np.atleast_1d('[1:512,2049:2079]'), # Taken from file headers
        )
        return detector_container.DetectorContainer(**detector_dict)


    # All the below is copied from not_nte as of August 2nd
    # Need to update for lco_floyds

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.

        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of ``PypeIt`` methods.
        """
        par = super().default_pypeit_par()

        # Turn off bias, illumflat, darks. Turn on overscan
        turn_off = dict(
                        #use_pixelflat=True,
                        #use_illumflat=False,
                        use_biasimage=False,
                        use_overscan=True,
                        use_darkimage=False,
                        trim=True,
                        )
        par.reset_all_processimages_par(**turn_off)

        # This works
        par['calibrations']['slitedges']['edge_thresh'] = 20.0
        par['calibrations']['slitedges']['fit_order'] = 3
        par['calibrations']['slitedges']['max_shift_adj'] = 0.5
        par['calibrations']['slitedges']['trace_thresh'] = 10
        par['calibrations']['slitedges']['fit_min_spec_length'] = 0.1
        par['calibrations']['slitedges']['det_min_spec_length'] = 0.1
        par['calibrations']['slitedges']['match_tol'] = 40
        par['calibrations']['slitedges']['det_buffer'] = 1
        par['calibrations']['slitedges']['auto_pca'] = True
        #par['calibrations']['slitedges']['left_right_pca'] = False
        #par['calibrations']['slitedges']['add_slits'] = "1:1650:365:470"
        par['calibrations']['slitedges']['sync_predict'] = "nearest"
        #par['calibrations']['slitedges']['smash_range'] = [0.3,0.7]
        #par['calibrations']['slitedges']['sobel_mode'] = "constant"
        #par['calibrations']['slitedges']['bound_detector'] = True

        # Think we need this - does not work
        #par["calibrations"]["arcframe"]["process"]["use_pixelflat"] = True

        # updating now.
        par['calibrations']['wavelengths']['lamps'] = ["HgAr_LCO"]
        par['calibrations']['wavelengths']['sigdetect'] = 5.0
        #par['calibrations']['wavelengths']['rms_thresh_frac_fwhm'] = 0.4
        #par['calibrations']['wavelengths']['fwhm'] = 5.0
        par['calibrations']['wavelengths']['n_final'] = [3,2]# [2, 4, 4, 4, 4, 4, 4, 4]
        par['calibrations']['wavelengths']['n_final'] = [5,3]# [2, 4, 4, 4, 4, 4, 4, 4]
        #par['calibrations']['wavelengths']['nreid_min'] = 1 # important
        
        par['calibrations']['wavelengths']['reference'] = 'arc'
        par['calibrations']['wavelengths']['reid_arxiv'] = 'lco_floyds_north.fits'
        par['calibrations']['wavelengths']['method'] = 'full_template'
        par['calibrations']['wavelengths']['nsnippet'] = 1 # important
        #par['calibrations']['wavelengths']['match_toler'] = 2.0
        
        #par['calibrations']['wavelengths']['reid_cont_sub'] = False

        #par['calibrations']['wavelengths']['sigrej_first'] = 1.0
        #par['calibrations']['wavelengths']['sigrej_final'] = 2.0

        # Echelle parameters
        #par['calibrations']['wavelengths']['echelle'] = True
        #par['calibrations']['wavelengths']['ech_nspec_coeff'] = 5
        #par['calibrations']['wavelengths']['ech_norder_coeff'] = 5
        #par['calibrations']['wavelengths']['ech_sigrej'] = 3.0
        #par['calibrations']['wavelengths']['ech_2dfit'] = True
        #par['calibrations']['wavelengths']['ech_sigrej'] = 3.0
        #par['calibrations']['wavelengths']['bad_orders_maxfrac'] = 0.8


        # tilts
        #par['calibrations']['tilts']['spat_order'] =  3
        
        # Flat
        par['calibrations']['flatfield']['slit_illum_finecorr'] = False # turn off for now

        # skysub
        par['reduce']['skysub']['bspline_spacing'] = 1

        # extraction
        par['reduce']['findobj']['maxnumber_sci'] = 1
        par['reduce']['findobj']['maxnumber_std'] = 1


        # Sensitivity function parameters
        par['sensfunc']['algorithm'] = 'IR'
        #par['sensfunc']['polyorder'] = [9, 11, 11, 9, 9, 8, 8, 7, 7, 7, 7, 7, 7, 7, 7]
        #par['sensfunc']['IR']['telgridfile'] = 'TelFit_Paranal_VIS_4900_11100_R25000.fits'

        return par

    @property
    def norders(self):
        """
        Number of orders observed for this spectograph.
        """
        return 2

    @property
    def order_spat_pos(self):
        """
        Return the expected spatial position of each echelle order.
        """

        return np.array([0.18,0.74])

        #np.array([91.470514, 378.12496]) were the positions used for lco floyds n
        # normalised by the detector height

    @property
    def orders(self):
        """
        Return the order number for each echelle order.
        """
        return np.arange(1, 3, 1, dtype=int) # orders 1 and 2

    @property
    def spec_min_max(self):
        """
        Return the minimum and maximum spectral pixel expected for the
        spectral range of each order.
        """
        spec_min = np.asarray([0,300])
        spec_max = np.asarray([2048,2048])
        return np.vstack((spec_min, spec_max))


    def order_platescale(self, order_vec, binning=None):
        """
        Return the platescale for each echelle order.

        This routine is only defined for echelle spectrographs, and it is
        undefined in the base class.

        Args:
            order_vec (`numpy.ndarray`_):
                The vector providing the order numbers.
            binning (:obj:`str`, optional):
                The string defining the spectral and spatial binning.

        Returns:
            `numpy.ndarray`_: An array with the platescale for each order
            provided by ``order``.
        """
        # No binning, but for an instrument with binning we would do this
        binspectral, binspatial = parse.parse_binning(binning)
        
        # Assume constant
        plate_scale = np.ones(2) * 0.34
        return plate_scale*binspatial

        # Not sure about this, commenting out
##    @property
##    def dloglam(self):
##        """
##        Return the logarithmic step in wavelength for output spectra.
##        """
##        # This number was computed by taking the mean of the dloglam for all
##        # the X-shooter orders. The specific loglam across the orders deviates
##        # from this value by +-7% from this first to final order. This is the
##        # unbinned value. It was actually measured to be 1.69207e-5 from a 2x1
##        # data and then divided by two.
##        return 8.46035e-06

    @property
    def loglam_minmax(self):
        """
        Return the base-10 logarithm of the first and last wavelength for
        ouput spectra.
        """
        return np.log10(3200), np.log10(10000)


















        
