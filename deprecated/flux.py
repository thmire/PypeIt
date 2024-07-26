""" Module for fluxing routines
"""
import glob
import numpy as np
import os
import scipy

from pkg_resources import resource_filename

from astropy import units
from astropy import constants
from astropy import coordinates
from astropy.table import Table, Column
from astropy.io import ascii
from astropy.io import fits
from astropy.stats import sigma_clipped_stats

from matplotlib import pyplot as plt

from linetools.spectra.xspectrum1d import XSpectrum1D

from pypeit.core import pydl
from pypeit import msgs
from pypeit import utils
from pypeit import debugger
from pypeit.wavemodel import conv2res

TINY = 1e-15
MAGFUNC_MAX = 25.0
MAGFUNC_MIN = -25.0
SN2_MAX = (20.0) ** 2
PYPEIT_FLUX_SCALE = 1e-17





def apply_sensfunc(spec_obj, sens_dict, airmass, exptime, extinct_correct=True, telluric_correct = False,
                   longitude=None, latitude=None):
    """ Apply the sensitivity function to the data
    We also correct for extinction.

    Parameters
    ----------
    spec_obj : dict
      SpecObj
    sens_dict : dict
      Sens Function dict
    airmass : float
      Airmass
    exptime : float
      Exposure time in seconds
    longitude : float
      longitude in degree for observatory
    latitude: float
      latitude in degree for observatory
      Used for extinction correction
    """

    # Loop on extraction modes
    for extract_type in ['boxcar', 'optimal']:
        extract = getattr(spec_obj, extract_type)
        if len(extract) == 0:
            continue
        msgs.info("Fluxing {:s} extraction for:".format(extract_type) + msgs.newline() + "{}".format(spec_obj))
        try:
            wave = np.copy(np.array(extract['WAVE_GRID']))
        except KeyError:
            wave = np.copy(np.array(extract['WAVE']))
        wave_sens = sens_dict['wave']
        sensfunc = sens_dict['sensfunc'].copy()

        # Did the user request a telluric correction from the same file?
        if telluric_correct and 'telluric' in sens_dict.keys():
            # This assumes there is a separate telluric key in this dict.
            telluric = sens_dict['telluric']
            msgs.info('Applying telluric correction')
            sensfunc = sensfunc*(telluric > 1e-10)/(telluric + (telluric < 1e-10))

        sensfunc_obs = scipy.interpolate.interp1d(wave_sens, sensfunc, bounds_error = False, fill_value='extrapolate')(wave)

        if extinct_correct:
            if longitude is None or latitude is None:
                msgs.error('You must specify longitude and latitude if we are extinction correcting')
            # Apply Extinction if optical bands
            msgs.info("Applying extinction correction")
            msgs.warn("Extinction correction applyed only if the spectra covers <10000Ang.")
            extinct = load_extinction_data(longitude,latitude)
            ext_corr = extinction_correction(wave* units.AA, airmass, extinct)
            senstot = sensfunc_obs * ext_corr
        else:
            senstot = sensfunc_obs.copy()

        flam = extract['COUNTS'] * senstot/ exptime
        flam_sig = (senstot/exptime)/ (np.sqrt(extract['COUNTS_IVAR']))
        flam_var = extract['COUNTS_IVAR'] / (senstot / exptime) **2

        # Mask bad pixels
        msgs.info(" Masking bad pixels")
        msk = np.zeros_like(senstot).astype(bool)
        msk[senstot <= 0.] = True
        msk[extract['COUNTS_IVAR'] <= 0.] = True
        flam[msk] = 0.
        flam_sig[msk] = 0.
        flam_var[msk] = 0.

        extract['FLAM'] = flam
        extract['FLAM_SIG'] = flam_sig
        extract['FLAM_IVAR'] = flam_var


def get_standard_spectrum(star_type=None, star_mag=None, ra=None, dec=None):
    '''
    Get the standard spetrum using given information of your standard/telluric star.

    Parameters:
      star_type: str
         Spectral type of your standard/telluric star
      star_mag: float
       Apparent magnitude of the telluric star
      ra: str
        Standard right-ascension in hh:mm:ss string format (e.g.,'05:06:36.6').
      dec: str
        Object declination in dd:mm:ss string format (e.g., 52:52:01.0')
    Return: dict
        Dictionary containing the information you provided and the standard/telluric spectrum.
    '''
    # Create star model
    if (ra is not None) and (dec is not None) and (star_mag is None) and (star_type is None):
        # Pull star spectral model from archive
        msgs.info("Get standard model")
        # Grab closest standard within a tolerance
        std_dict = find_standard_file(ra, dec)
        if std_dict is not None:
            # Load standard
            load_standard_file(std_dict)
        else:
            msgs.error('No spectrum found in our database for your standard star. Please use another standard star \
                       or consider add it into out database.')
    elif (star_mag is not None) and (star_type is not None):
        ## using vega spectrum
        if 'A0' in star_type:
            msgs.info('Using vega spectrum to correct telluric')
            std_dict={'stellar_type':star_type , 'Vmag': star_mag}

            ## TODO: we should get a higher resolution spectra and then convolve it to the resolution of your data!
            ##       the Tspec one works well for GNIRS and NIRES, but apparently not X-SHOOTER.
            ## vega spectrum from STSCI
            #vega_file = resource_filename('pypeit', '/data/standards/vega_04_to_06.dat')
            #vega_data = Table.read(vega_file, comment='#', format='ascii')
            #std_dict = dict(cal_file='Vega_04_to_06', name=star_type, std_source='vega',
            #                std_ra=None, std_dec=None)
            #std_dict['wave'] = vega_data['col1'] * units.AA
            ## Vega model from TSPECTOOL
            vega_file = resource_filename('pypeit', '/data/standards/vega_tspectool_vacuum.dat')
            vega_data = Table.read(vega_file, comment='#', format='ascii')
            std_dict = dict(cal_file='vega_tspectool_vacuum', name=star_type, std_source='vega',
                            std_ra=None, std_dec=None)
            std_dict['wave'] = vega_data['col1'] * units.AA

            # ToDo: we should redden the vega spectra based on E(B-V) of the telluric star
            # vega is V=0.03
            std_dict['flux'] = vega_data['col2'] * 10**(0.4*(0.03-star_mag)) / PYPEIT_FLUX_SCALE * \
                               units.erg / units.s / units.cm ** 2 / units.AA
        ## using Kurucz stellar model
        else:
            # Create star spectral model
            msgs.info("Creating standard model")
            # Create star model
            star_loglam, star_flux, std_dict = telluric_sed(star_mag, star_type)
            star_lam = 10 ** star_loglam
            # Generate a dict matching the output of find_standard_file
            std_dict = dict(cal_file='KuruczTelluricModel', name=star_type, std_source='KuruczModel',
                            std_ra=None, std_dec=None)
            std_dict['wave'] = star_lam * units.AA
            std_dict['flux'] = star_flux / PYPEIT_FLUX_SCALE * units.erg / units.s / units.cm ** 2 / units.AA
    else:
        debugger.set_trace()
        msgs.error('Insufficient information provided for fluxing. '
                   'Either the coordinates of the standard or a stellar type and magnitude are needed.')


    return std_dict

def generate_sensfunc(wave, counts, counts_ivar, airmass, exptime, longitude, latitude, telluric=True, star_type=None,
                      star_mag=None, ra=None, dec=None, std_file = None, poly_norder=4, BALM_MASK_WID=5., nresln=20.,
                      resolution=3000.,trans_thresh=0.9,polycorrect=True, debug=False):

    """ Function to generate the sensitivity function.
        This can work in different regimes:
    - If telluric=False and RA=None and Dec=None
      the code creates a synthetic standard star spectrum (or VEGA spectrum if type=A0) using the Kurucz models,
      and from this it generates a sens func using nresln=20.0 and masking out telluric regions.
    - If telluric=False and RA and Dec are assigned
      the standard star spectrum is extracted from the archive, and a sens func
      is generated using nresln=20.0 and masking out telluric regions.
    - If telluric=True
      the code creates a sintetic standard star spectrum  (or VEGA spectrum if type=A0) using the Kurucz models,
      the sens func is a pixelized sensfunc (not smooth) for correcting both throughput and telluric lines.
      if you set polycorrect=True, the sensfunc in the Hydrogen recombination line region (often seen in star spectra)
      will be replaced by a smoothed polynomial function.

    Parameters:
    ----------
    wave : array
      Wavelength of the star [no longer with units]
    counts : array
      Flux (in counts) of the star
    counts_ivar : array
      Inverse variance of the star
    airmass : float
      Airmass
    exptime : float
      Exposure time in seconds
    spectrograph : dict
      Instrument specific dict
      Used for extinction correction
    telluric : bool
      if True performs a telluric correction
    star_type : str
      Spectral type of the telluric star (used if telluric=True)
    star_mag : float
      Apparent magnitude of telluric star (used if telluric=True)
    RA : float
      deg, RA of the telluric star
      if assigned, the standard star spectrum will be extracted from
      the archive
    DEC : float
      deg, DEC of the telluric star
      if assigned, the standard star spectrum will be extracted from
      the archive
    BALM_MASK_WID : float
      Mask parameter for Balmer absorption. A region equal to
      BALM_MASK_WID*resln is masked where resln is the estimate
      for the spectral resolution.
    polycorrect: bool
      Whether you want to correct the sensfunc with polynomial in the Balmer absortion line regions
    poly_norder: int
      Order number of polynomial fit.

    Returns:
    -------
    sens_dict : dict
      sensitivity function described by a dict
    """
    # Create copy of the arrays to avoid modification and convert to
    # electrons / s
    wave_star = wave.copy()
    flux_star = counts.copy() / exptime
    ivar_star = counts_ivar.copy() * exptime ** 2

    # Units
    if not isinstance(wave_star, units.Quantity):
        wave_star = wave_star * units.AA

    # Extinction correction
    msgs.info("Applying extinction correction")
    extinct = load_extinction_data(longitude,latitude)
    ext_corr = extinction_correction(wave * units.AA, airmass, extinct)
    # Correct for extinction
    flux_star = flux_star * ext_corr
    ivar_star = ivar_star / ext_corr ** 2

    std_dict =  get_standard_spectrum(star_type=star_type, star_mag=star_mag, ra=ra, dec=dec)
    # Interpolate the standard star onto the current set of observed wavelengths
    flux_true = scipy.interpolate.interp1d(std_dict['wave'], std_dict['flux'],bounds_error=False,
                                           fill_value='extrapolate')(wave_star)
    # Do we need to extrapolate? TODO Replace with a model or a grey body?
    if np.min(flux_true) <= 0.:
        msgs.warn('Your spectrum extends beyond calibrated standard star, extrapolating the spectra with polynomial.')
        mask_model = flux_true <= 0
        msk_poly, poly_coeff = utils.robust_polyfit_djs(std_dict['wave'].value, std_dict['flux'].value,8,function='polynomial',
                                                    invvar=None, guesses=None, maxiter=50, inmask=None, \
                                                    lower=3.0, upper=3.0, maxdev=None, maxrej=3, groupdim=None,
                                                    groupsize=None,groupbadpix=False, grow=0, sticky=True, use_mad=True)
        star_poly = utils.func_val(poly_coeff, wave_star.value, 'polynomial')
        #flux_true[mask_model] = star_poly[mask_model]
        flux_true = star_poly.copy()
        if debug:
            plt.plot(std_dict['wave'], std_dict['flux'],'bo',label='Raw Star Model')
            plt.plot(std_dict['wave'],  utils.func_val(poly_coeff, std_dict['wave'].value, 'polynomial'), 'k-',label='robust_poly_fit')
            plt.plot(wave_star,flux_true,'r-',label='Your Final Star Model used for sensfunc')
            plt.show()

    # Get masks from observed star spectrum. True = Good pixels
    msk_bad, msk_star, msk_tell = get_mask(wave_star.value, flux_star, ivar_star, mask_star=True, mask_tell=True,
                                           BALM_MASK_WID=BALM_MASK_WID, trans_thresh=0.9)

    # Get sensfunc
    LBLRTM = False
    if LBLRTM:
        # sensfunc = lblrtm_sensfunc() ???
        msgs.develop('fluxing and telluric correction based on LBLRTM model is under developing.')
    else:
        sensfunc, mask_sens = standard_sensfunc(wave_star.value, flux_star, ivar_star, flux_true, msk_bad=msk_bad, msk_star=msk_star,
                                msk_tell=msk_tell, maxiter=35,upper=3.0, lower=3.0, poly_norder=poly_norder,
                                BALM_MASK_WID=BALM_MASK_WID,nresln=nresln,telluric=telluric, resolution=resolution,
                                polycorrect= polycorrect, debug=debug, show_QA=False)

    if debug:
        plt.plot(wave_star.value[mask_sens], flux_true[mask_sens], color='k',lw=2,label='Reference Star')
        plt.plot(wave_star.value[mask_sens], flux_star[mask_sens]*sensfunc[mask_sens], color='r',label='Fluxed Observed Star')
        plt.xlabel(r'Wavelength [$\AA$]')
        plt.ylabel('Flux [erg/s/cm2/Ang.]')
        plt.legend(fancybox=True, shadow=True)
        plt.show()


    # Add in wavemin,wavemax
    sens_dict = {}
    sens_dict['wave'] = wave_star.value
    sens_dict['sensfunc'] = sensfunc
    sens_dict['wave_min'] = np.min(wave_star)
    sens_dict['wave_max'] = np.max(wave_star)
    sens_dict['exptime']= exptime
    sens_dict['airmass']= airmass
    sens_dict['std_file']= std_file
    # Get other keys from standard dict
    sens_dict['std_ra'] = std_dict['std_ra']
    sens_dict['std_dec'] = std_dict['std_dec']
    sens_dict['std_name'] = std_dict['name']
    sens_dict['cal_file'] = std_dict['cal_file']
    sens_dict['flux_true'] = flux_true
    sens_dict['mask_sens'] = mask_sens
    #sens_dict['std_dict'] = std_dict
    #sens_dict['msk_star'] = msk_star
    #sens_dict['mag_set'] = mag_set

    return sens_dict

def get_mask(wave_star,flux_star, ivar_star, mask_star=True, mask_tell=True, BALM_MASK_WID=10., trans_thresh=0.9):
    '''
    Get couple of masks from your observed standard spectrum.
    Parameters
    ----------
      wave_star: numpy array
        wavelength array of your spectrum
      flux_star: numpy array
        flux array of your spectrum
      ivar_star:
        ivar array of your spectrum
      mask_star: bool
        whether you need to mask Hydrogen recombination line region. If False, the returned msk_star are all good.
      mask_tell: bool
        whether you need to mask telluric region. If False, the returned msk_tell are all good.
      trans_thresh: float
        parameter for selecting telluric regions.
    returns
    ----------
      msk_bad: bool type numpy array
        mask for bad pixels.
      msk_star: bool type numpy array
        mask for recombination lines in star spectrum.
      msk_tell: bool type numpy array
        mask for telluric regions.
    '''

    # Mask (True = good pixels)
    # mask for bad pixels
    msk_bad = np.ones_like(flux_star).astype(bool)
    # mask for recombination lines
    msk_star = np.ones_like(flux_star).astype(bool)
    # mask for telluric regions
    msk_tell = np.ones_like(flux_star).astype(bool)

    # masking bad entries
    msgs.info(" Masking bad pixels")
    msk_bad[ivar_star <= 0.] = False
    msk_bad[flux_star <= 0.] = False
    # Mask edges
    msgs.info(" Masking edges")
    msk_bad[:1] = False
    msk_bad[-1:] = False
    # Mask Atm. cutoff
    msgs.info(" Masking Below the atmospheric cutoff")
    atms_cutoff = wave_star <= 3000.0
    msk_bad[atms_cutoff] = False

    if mask_star:
        # Mask Balmer, Paschen, Brackett, and Pfund recombination lines
        msgs.info("Masking recombination lines:")
        # Mask Balmer
        msgs.info(" Masking Balmer")
        lines_balm = np.array([3836.4, 3969.6, 3890.1, 4102.8, 4102.8, 4341.6, 4862.7, 5407.0,
                               6564.6, 8224.8, 8239.2])
        for line_balm in lines_balm:
            ibalm = np.abs(wave_star - line_balm) <= BALM_MASK_WID
            msk_star[ibalm] = False
        # Mask Paschen
        msgs.info(" Masking Paschen")
        # air wavelengths from:
        # https://www.subarutelescope.org/Science/Resources/lines/hi.html
        lines_pasc = np.array([8203.6, 8440.3, 8469.6, 8504.8, 8547.7, 8600.8, 8667.4, 8752.9,
                               8865.2, 9017.4, 9229.0, 9546.0, 10049.4, 10938.1,
                               12818.1, 18751.0])
        for line_pasc in lines_pasc:
            ipasc = np.abs(wave_star - line_pasc) <= BALM_MASK_WID
            msk_star[ipasc] = False
        # Mask Brackett
        msgs.info(" Masking Brackett")
        # air wavelengths from:
        # https://www.subarutelescope.org/Science/Resources/lines/hi.html
        lines_brac = np.array([14584.0, 18174.0, 19446.0, 21655.0,26252.0, 40512.0])
        for line_brac in lines_brac:
            ibrac = np.abs(wave_star - line_brac) <= BALM_MASK_WID
            msk_star[ibrac] = False
        # Mask Pfund
        msgs.info(" Masking Pfund")
        # air wavelengths from:
        # https://www.subarutelescope.org/Science/Resources/lines/hi.html
        lines_pfund = np.array([22788.0, 32961.0, 37395.0, 46525.0,74578.0])
        for line_pfund in lines_pfund:
            ipfund = np.abs(wave_star - line_pfund) <= BALM_MASK_WID
            msk_star[ipfund] = False

    if mask_tell:
        ## Mask telluric region in the optical
        tell_opt = np.any([((wave_star >= 6270.00) & (wave_star <= 6290.00)), # H2O
                       ((wave_star >= 6850.00) & (wave_star <= 6960.00)), #O2 telluric band
                       ((wave_star >= 7580.00) & (wave_star <= 7750.00)), #O2 telluric band
                       ((wave_star >= 7160.00) & (wave_star <= 7340.00)), #H2O
                       ((wave_star >= 8150.00) & (wave_star <= 8250.00))],axis=0) #H2O
        msk_tell[tell_opt] = False
        ## Mask near-infrared telluric region
        if np.max(wave_star)>9100.0:
            # ToDo: should use the specific atmosphere transmission after FBD get the grid.
            ## Read atmosphere transmission
            '''
            if watervp <1.5:
                skytrans_file = resource_filename('pypeit', '/data/skisim/'+'mktrans_zm_10_10.dat')
            elif (watervp>=1.5 and watervp<2.3):
                skytrans_file = resource_filename('pypeit', '/data/skisim/'+'mktrans_zm_16_10.dat')
            elif (watervp>=2.3 and watervp<4.0):
                skytrans_file = resource_filename('pypeit', '/data/skisim/' + 'mktrans_zm_30_10.dat')
            else:
                skytrans_file = resource_filename('pypeit', '/data/skisim/' + 'mktrans_zm_50_10.dat')
            '''
            skytrans_file = resource_filename('pypeit', '/data/skisim/' + 'mktrans_zm_10_10.dat')
            skytrans = ascii.read(skytrans_file)
            wave_trans, trans = skytrans['wave']*10000.0, skytrans['trans']
            trans_use = (wave_trans>=np.min(wave_star)-100.0) & (wave_trans<=np.max(wave_star)+100.0)
            # Estimate the resolution of your spectra.
            # I assumed 3 pixels per resolution. This gives an approximate right resolution at the middle point.
            resolution = np.median(wave_star) / np.median(wave_star - np.roll(wave_star, 1)) / 3
            trans_convolved, px_sigma, px_bin = conv2res(wave_trans[trans_use], trans[trans_use], resolution,
                                                         central_wl='midpt', debug=False)
            trans_final = scipy.interpolate.interp1d(wave_trans[trans_use], trans_convolved,bounds_error=False,
                                                     fill_value='extrapolate')(wave_star)
            tell_nir = (trans_final<trans_thresh) & (wave_star>9100.0)
            msk_tell[tell_nir] = False
        else:
            msgs.info('Your spectrum is bluer than 9100A, only optical telluric regions are masked.')

    return msk_bad, msk_star, msk_tell


def standard_sensfunc(wave, flux, ivar, flux_std, msk_bad=None, msk_star=None, msk_tell=None,
                 maxiter=35, upper=2, lower=2, poly_norder=5, BALM_MASK_WID=50., nresln=20., telluric=True,
                 resolution=2700., polycorrect=True, debug=False, show_QA=False):
    """
    Generate a sensitivity function based on observed flux and standard spectrum.

    Parameters
    ----------
    wave : ndarray
      wavelength as observed
    flux : ndarray
      counts/s as observed
    ivar : ndarray
      inverse variance
    flux_std : Quantity array
      standard star true flux (erg/s/cm^2/A)
    msk_bad : ndarray
      mask for bad pixels. True is good.
    msk_star: ndarray
      mask for hydrogen recombination lines. True is good.
    msk_tell:ndarray
      mask for telluric regions. True is good.
    maxiter : integer
      maximum number of iterations for polynomial fit
    upper : integer
      number of sigma for rejection in polynomial
    lower : integer
      number of sigma for rejection in polynomial
    poly_norder : integer
      order of polynomial fit
    BALM_MASK_WID : float
      in units of angstrom
      Mask parameter for Balmer absorption. A region equal to
      BALM_MASK_WID is masked.
    resolution: integer/float.
      spectra resolution
      This paramters should be removed in the future. The resolution should be estimated from spectra directly.
    debug : bool
      if True shows some dubugging plots

    Returns
    -------
    sensfunc
    """
    # Create copy of the arrays to avoid modification
    wave_obs = wave.copy()
    flux_obs = flux.copy()
    ivar_obs = ivar.copy()
    # preparing arrays
    if np.all(~np.isfinite(ivar_obs)):
        msgs.warn("NaN are present in the inverse variance")

    # check masks
    if msk_bad is None:
        msk_bad = np.ones_like(wave_obs,dtype=bool)
    if msk_tell is None:
        msk_tell = np.ones_like(wave_obs,dtype=bool)
    if msk_star is None:
        msk_star = np.ones_like(wave_obs, dtype=bool)

    # Removing outliers
    # Calculate log of flux_obs setting a floor at TINY
    logflux_obs = 2.5 * np.log10(np.maximum(flux_obs, TINY))
    # Set a fix value for the variance of logflux
    logivar_obs = np.ones_like(logflux_obs) * (10.0 ** 2)
    # Calculate log of flux_std model setting a floor at TINY
    logflux_std = 2.5 * np.log10(np.maximum(flux_std, TINY))
    # Calculate ratio setting a floor at MAGFUNC_MIN and a ceiling at
    # MAGFUNC_MAX
    magfunc = logflux_std - logflux_obs
    magfunc = np.maximum(np.minimum(magfunc, MAGFUNC_MAX), MAGFUNC_MIN)
    msk_magfunc = (magfunc < 0.99 * MAGFUNC_MAX) & (magfunc > 0.99 * MAGFUNC_MIN)

    # Define two new masks, True is good and False is masked pixel
    # mask for all bad pixels on sensfunc
    masktot = msk_bad & msk_magfunc & np.isfinite(ivar_obs) & np.isfinite(logflux_obs) & np.isfinite(logflux_std)
    logivar_obs[np.invert(masktot)] = 0.0
    # mask used for polynomial fit
    msk_fit_sens = masktot & msk_tell & msk_star

    # Polynomial fitting to derive a smooth sensfunc (i.e. without telluric)
    _, poly_coeff = utils.robust_polyfit_djs(wave_obs[msk_fit_sens], magfunc[msk_fit_sens], poly_norder, \
                                                    function='polynomial', invvar=None, guesses=None, maxiter=maxiter, \
                                                    inmask=None, lower=lower, upper=upper, maxdev=None, \
                                                    maxrej=None, groupdim=None, groupsize=None, groupbadpix=False, \
                                                    grow=0, sticky=True, use_mad=True)
    magfunc_poly = utils.func_val(poly_coeff, wave_obs, 'polynomial')

    # Polynomial corrections on Hydrogen Recombination lines
    if ((sum(msk_fit_sens) > 0.5 * len(msk_fit_sens)) & polycorrect):
        ## Only correct Hydrogen Recombination lines with polyfit in the telluric free region
        balmer_clean = np.zeros_like(wave_obs, dtype=bool)
        # Commented out the bluest recombination lines since they are weak for spectroscopic standard stars.
        #836.4, 3969.6, 3890.1, 4102.8, 4102.8, 4341.6, 4862.7,   \
        lines_hydrogen = np.array([5407.0, 6564.6, 8224.8, 8239.2, 8203.6, 8440.3, 8469.6, 8504.8, 8547.7, 8600.8, \
                                   8667.4, 8752.9, 8865.2, 9017.4, 9229.0, 10049.4, 10938.1, 12818.1, 21655.0])
        for line_hydrogen in lines_hydrogen:
            ihydrogen = np.abs(wave_obs - line_hydrogen) <= BALM_MASK_WID
            balmer_clean[ihydrogen] = True
        msk_clean = ((balmer_clean) | (magfunc == MAGFUNC_MAX) | (magfunc == MAGFUNC_MIN)) & \
                    (magfunc_poly > MAGFUNC_MIN) & (magfunc_poly < MAGFUNC_MAX)
        magfunc[msk_clean] = magfunc_poly[msk_clean]
        msk_badpix = np.isfinite(ivar_obs) & (ivar_obs > 0)
        magfunc[~msk_badpix] = magfunc_poly[~msk_badpix]
    else:
        ## if half more than half of your spectrum is masked (or polycorrect=False) then do not correct it with polyfit
        msgs.warn('No polynomial corrections performed on Hydrogen Recombination line regions')

    if not telluric:
        # Apply mask to ivar
        #logivar_obs[~msk_fit_sens] = 0.

        # ToDo
        # Compute an effective resolution for the standard. This could be improved
        # to setup an array of breakpoints based on the resolution. At the
        # moment we are using only one number
        msgs.work("Should pull resolution from arc line analysis")
        msgs.work("At the moment the resolution is taken as the PixelScale")
        msgs.work("This needs to be changed!")
        std_pix = np.median(np.abs(wave_obs - np.roll(wave_obs, 1)))
        std_res = np.median(wave_obs/resolution) # median resolution in units of Angstrom.
        #std_res = std_pix
        #resln = std_res
        if (nresln * std_res) < std_pix:
            msgs.warn("Bspline breakpoints spacing shoud be larger than 1pixel")
            msgs.warn("Changing input nresln to fix this")
            nresln = std_res / std_pix

        # Fit magfunc with bspline
        kwargs_bspline = {'bkspace': std_res * nresln}
        kwargs_reject = {'maxrej': 5}
        msgs.info("Initialize bspline for flux calibration")
        init_bspline = pydl.bspline(wave_obs, bkspace=kwargs_bspline['bkspace'])
        fullbkpt = init_bspline.breakpoints

        # TESTING turning off masking for now
        # remove masked regions from breakpoints
        msk_obs = np.ones_like(wave_obs).astype(bool)
        msk_obs[~masktot] = False
        msk_bkpt = scipy.interpolate.interp1d(wave_obs, msk_obs, kind='nearest', fill_value='extrapolate')(fullbkpt)
        init_breakpoints = fullbkpt[msk_bkpt > 0.999]

        # init_breakpoints = fullbkpt
        msgs.info("Bspline fit on magfunc. ")
        bset1, bmask = pydl.iterfit(wave_obs, magfunc, invvar=logivar_obs, inmask=msk_fit_sens, upper=upper, lower=lower,
                                    fullbkpt=init_breakpoints, maxiter=maxiter, kwargs_bspline=kwargs_bspline,
                                    kwargs_reject=kwargs_reject)
        logfit1, _ = bset1.value(wave_obs)
        logfit_bkpt, _ = bset1.value(init_breakpoints)

        if debug:
            # Check for calibration
            plt.figure(1)
            plt.plot(wave_obs, magfunc, drawstyle='steps-mid', color='black', label='magfunc')
            plt.plot(wave_obs, logfit1, color='cornflowerblue', label='logfit1')
            plt.plot(wave_obs[~msk_fit_sens], magfunc[~msk_fit_sens], '+', color='red', markersize=5.0,
                     label='masked magfunc')
            plt.plot(wave_obs[~msk_fit_sens], logfit1[~msk_fit_sens], '+', color='red', markersize=5.0,
                     label='masked logfit1')
            plt.plot(init_breakpoints, logfit_bkpt, '.', color='green', markersize=4.0, label='breakpoints')
            plt.plot(init_breakpoints, np.interp(init_breakpoints, wave_obs, magfunc), '.', color='green',
                     markersize=4.0,
                     label='breakpoints')
            plt.plot(wave_obs, 1.0 / np.sqrt(logivar_obs), color='orange', label='sigma')
            plt.legend()
            plt.xlabel('Wavelength [ang]')
            plt.ylim(0.0, 1.2 * MAGFUNC_MAX)
            plt.title('1st Bspline fit')
            plt.show()
        # Create sensitivity function
        magfunc = np.maximum(np.minimum(logfit1, MAGFUNC_MAX), MAGFUNC_MIN)
        if ((sum(msk_fit_sens) > 0.5 * len(msk_fit_sens)) & polycorrect):
            msk_clean = ((magfunc==MAGFUNC_MAX) | (magfunc==MAGFUNC_MIN)) & \
                        (magfunc_poly>MAGFUNC_MIN) & (magfunc_poly<MAGFUNC_MAX)
            magfunc[msk_clean] = magfunc_poly[msk_clean]
            msk_badpix = np.isfinite(ivar_obs)& (ivar_obs>0)
            magfunc[~msk_badpix] = magfunc_poly[~msk_badpix]
        else:
            ## if half more than half of your spectrum is masked (or polycorrect=False) then do not correct it with polyfit
            msgs.warn('No polynomial corrections performed on Hydrogen Recombination line regions')

    # Calculate sensfunc
    sensfunc = 10.0 ** (0.4 * magfunc)

    if debug:
        plt.figure()
        magfunc_raw = logflux_std - logflux_obs
        plt.plot(wave_obs[masktot],magfunc_raw[masktot] , 'k-',lw=3,label='Raw Magfunc')
        plt.plot(wave_obs[masktot],magfunc_poly[masktot] , 'c-',lw=3,label='Polynomial Fit')
        plt.plot(wave_obs[np.invert(msk_tell)], magfunc_raw[np.invert(msk_tell)], 's',
                 color='0.7',label='Telluric Region')
        plt.plot(wave_obs[np.invert(msk_star)], magfunc_raw[np.invert(msk_star)], 'r+',label='Recombination Line region')
        plt.plot(wave_obs[masktot], magfunc[masktot],'b-',label='Final Magfunc')
        plt.legend(fancybox=True, shadow=True)
        plt.xlim([0.995*np.min(wave_obs[masktot]),1.005*np.max(wave_obs[masktot])])
        plt.ylim([0.,1.2*np.max(magfunc[masktot])])
        plt.show()
        plt.close()

    return sensfunc, masktot

def extinction_correction(wave, airmass, extinct):
    """
    Derive extinction correction
    Based on algorithm in LowRedux (long_extinct)

    Parameters
    ----------
    wave : ndarray
      Wavelengths for interpolation. Should be sorted
      Assumes Angstroms
    airmass : float
      Airmass
    extinct : Table
      Table of extinction values

    Returns
    -------
    flux_corr : ndarray
      Flux corrections at the input wavelengths
    """
    # Checks
    if airmass < 1.:
        msgs.error("Bad airmass value in extinction_correction")
    # Interpolate
    f_mag_ext = scipy.interpolate.interp1d(extinct['wave'],extinct['mag_ext'], bounds_error=False, fill_value=0.)
    mag_ext = f_mag_ext(wave)#.to('AA').value)

    # Deal with outside wavelengths
    gdv = np.where(mag_ext > 0.)[0]

    if len(gdv) == 0:
        msgs.warn("No valid extinction data available at this wavelength range. Extinction correction not applied")
    elif gdv[0] != 0:  # Low wavelengths
        mag_ext[0:gdv[0]] = mag_ext[gdv[0]]
        msgs.warn("Extrapolating at low wavelengths using last valid value")
    elif gdv[-1] != (mag_ext.size - 1):  # High wavelengths
        mag_ext[gdv[-1] + 1:] = mag_ext[gdv[-1]]
        msgs.warn("Extrapolating at high wavelengths using last valid value")
    else:
        msgs.info("Extinction data covered the whole spectra. Correct it!")
    # Evaluate
    flux_corr = 10.0 ** (0.4 * mag_ext * airmass)
    # Return
    return flux_corr


def find_standard_file(ra, dec, toler=20.*units.arcmin, check=False):
    """
    Find a match for the input file to one of the archived
    standard star files (hopefully).  Priority is by order of search.

    Args:
        ra (str or float):
            Object right-ascension in hh:mm:ss string format or in degrees (e.g.,
            '05:06:36.6' or 76.6525). astropy.coordinates.SkyCoord is used to parse the coordinates
        dec (str or float):
            Object declination in dd:mm:ss string format or in degrees (e.g.,
            52:52:01.0' or 52.86694) astropy.coordinates.SkyCoord is used to parse the coordinates
        toler (:class:`astropy.units.quantity.Quantity`, optional):
            Tolerance on matching archived standards to input.  Expected
            to be in arcmin.
        check (:obj:`bool`, optional):
            If True, the routine will only check to see if a standard
            star exists within the input ra, dec, and toler range.

    Returns:
        dict, bool: If check is True, return True or False depending on
        if the object is matched to a library standard star.  If check
        is False and no match is found, return None.  Otherwise, return
        a dictionary with the matching standard star with the following
        meta data::
            - 'file': str -- Filename
              table
            - 'name': str -- Star name
            - 'ra': str -- RA(J2000)
            - 'dec': str -- DEC(J2000)
    """
    # Priority
    std_sets = [load_xshooter, load_calspec, load_esofil]
    std_file_source = ['xshooter', 'calspec', 'eso']  # XSHOOTER ASCII format; Calspec style FITS binary table; ESO ASCII format.

    # SkyCoord
    try:
        ra, dec = float(ra), float(dec)
        obj_coord = coordinates.SkyCoord(ra, dec, unit=(units.deg, units.deg))
    except:
        obj_coord = coordinates.SkyCoord(ra, dec, unit=(units.hourangle, units.deg))

    # Loop on standard sets
    closest = dict(sep=999 * units.deg)
    for qq, sset in enumerate(std_sets):
        # Stars
        path, star_tbl = sset()
        star_coords = coordinates.SkyCoord(star_tbl['RA_2000'], star_tbl['DEC_2000'],
                                           unit=(units.hourangle, units.deg))
        # Match
        idx, d2d, d3d = coordinates.match_coordinates_sky(obj_coord, star_coords, nthneighbor=1)
        if d2d < toler:
            if check:
                return True
            else:
                # Generate a dict
                _idx = int(idx)
                std_dict = dict(cal_file=os.path.join(path,star_tbl[_idx]['File']),
                                name=star_tbl[_idx]['Name'], std_source=std_file_source[qq],
                                std_ra=star_tbl[_idx]['RA_2000'],
                                std_dec=star_tbl[_idx]['DEC_2000'])
                # Return
                msgs.info("Using standard star {:s}".format(std_dict['name']))
                return std_dict
        else:
            # Save closest found so far
            imind2d = np.argmin(d2d)
            mind2d = d2d[imind2d]
            if mind2d < closest['sep']:
                closest['sep'] = mind2d
                closest.update(dict(name=star_tbl[int(idx)]['Name'],
                                    ra=star_tbl[int(idx)]['RA_2000'],
                                    dec=star_tbl[int(idx)]['DEC_2000']))

    # Standard star not found
    if check:
        return False

    msgs.warn("No standard star was found within a tolerance of {:g}".format(toler))
    msgs.info("Closest standard was {:s} at separation {:g}".format(closest['name'],
                                                                    closest['sep'].to('arcmin')))
    msgs.warn("Flux calibration will not be performed")
    return None


def load_calspec():
    """
    Load the list of calspec standards

    Parameters
    ----------

    Returns
    -------
    calspec_path : str
      Path from pypeitdir to calspec standard star files
    calspec_stds : Table
      astropy Table of the calspec standard stars (file, Name, RA, DEC)
    """
    # Read
    calspec_path = 'data/standards/calspec/'
    calspec_file = resource_filename('pypeit', calspec_path + 'calspec_info.txt')
    calspec_stds = Table.read(calspec_file, comment='#', format='ascii')
    # Return
    return calspec_path, calspec_stds

def load_esofil():
    """
    Load the list of ESO standards

    Parameters
    ----------

    Returns
    -------
    esofil_path : str
      Path from pypeitdir to calspec standard star files
    esofil_stds : Table
      astropy Table of the calspec standard stars (file, Name, RA, DEC)
    """
    # Read
    esofil_path = 'data/standards/ESOFIL/'
    esofil_file = resource_filename('pypeit', esofil_path + 'esofil_info.txt')
    esofil_stds = Table.read(esofil_file, comment='#', format='ascii')
    # Return
    return esofil_path, esofil_stds

def load_xshooter():
    """
    Load the list of ESO standards

    Parameters
    ----------

    Returns
    -------
    esofil_path : str
      Path from pypeitdir to calspec standard star files
    esofil_stds : Table
      astropy Table of the calspec standard stars (file, Name, RA, DEC)
    """
    # Read
    xshooter_path = 'data/standards/xshooter/'
    xshooter_file = resource_filename('pypeit', xshooter_path + 'xshooter_info.txt')
    xshooter_stds = Table.read(xshooter_file, comment='#', format='ascii')
    # Return
    return xshooter_path, xshooter_stds

def load_extinction_data(longitude, latitude, toler=5. * units.deg):
    """
    Find the best extinction file to use, based on longitude and latitude
    Loads it and returns a Table

    Parameters
    ----------
    toler : Angle, optional
      Tolerance for matching detector to site (5 deg)

    Returns
    -------
    ext_file : Table
      astropy Table containing the 'wavelength', 'extinct' data for AM=1.
    """
    # Mosaic coord
    mosaic_coord = coordinates.SkyCoord(longitude, latitude, frame='gcrs', unit=units.deg)
    # Read list
    extinct_path = resource_filename('pypeit', '/data/extinction/')
    extinct_summ = extinct_path + 'README'
    extinct_files = Table.read(extinct_summ, comment='#', format='ascii')
    # Coords
    ext_coord = coordinates.SkyCoord(extinct_files['Lon'], extinct_files['Lat'], frame='gcrs',
                                     unit=units.deg)
    # Match
    idx, d2d, d3d = coordinates.match_coordinates_sky(mosaic_coord, ext_coord, nthneighbor=1)
    if d2d < toler:
        extinct_file = extinct_files[int(idx)]['File']
        msgs.info("Using {:s} for extinction corrections.".format(extinct_file))
    else:
        msgs.warn("No file found for extinction corrections.  Applying none")
        msgs.warn("You should generate a site-specific file")
        return None
    # Read
    extinct = Table.read(extinct_path + extinct_file, comment='#', format='ascii',
                         names=('iwave', 'mag_ext'))
    wave = Column(np.array(extinct['iwave']) * units.AA, name='wave')
    extinct.add_column(wave)
    # Return
    return extinct[['wave', 'mag_ext']]

def load_filter_file(filter):
    """
    Load a system response curve for a given filter

    Args:
        filter (str): Name of filter

    Returns:
        ndarray, ndarray: wavelength, instrument throughput

    """
    '''
    # Optical filters
    BASS_MZLS_filters = ['BASS-MZLS-{}'.format(i) for i in ['G', 'R','Z']]
    CFHT_filters = ['CFHT-{}'.format(i) for i in ['U', 'G', 'R', 'I', 'Z']]
    DECAM_filters = ['DECAM-{}'.format(i) for i in ['U', 'G', 'R', 'I', 'Z', 'Y']]
    HSC_filters = ['HSC-{}'.format(i) for i in ['G', 'R', 'I', 'Z', 'Y']]
    LSST_filters = ['LSST-{}'.format(i) for i in ['U', 'G', 'R', 'I', 'Z', 'Y']]
    PS1_filters = ['PS1-{}'.format(i) for i in ['G', 'R', 'I', 'Z', 'Y']]
    SDSS_filters = ['SDSS-{}'.format(i) for i in ['U', 'G', 'R', 'I', 'Z']]

    # NIR filters
    UKIDSS_filters = ['UKIRT-{}'.format(i) for i in ['Y', 'J', 'H', 'K']]
    VISTA_filters = ['VISTA-{}'.format(i) for i in ['Z', 'Y', 'J', 'H', 'K']]
    TMASS_filters = ['TMASS-{}'.format(i) for i in ['J', 'H', 'K']]

    # Other filters
    GAIA_filters = ['GAIA-{}'.format(i) for i in ['G', 'B', 'R']]
    GALEX_filters = ['GALEX-{}'.format(i) for i in ['F', 'N']]
    WISE_filters = ['WISE-{}'.format(i) for i in ['W1', 'W2', 'W3', 'W4']]

    allowed_options = BASS_MZLS_filters + CFHT_filters + DECAM_filters + HSC_filters \
                      + LSST_filters + PS1_filters + SDSS_filters + UKIDSS_filters\
                      + VISTA_filters + TMASS_filters + GAIA_filters + GALEX_filters + WISE_filters
    '''
    filter_file = resource_filename('pypeit', os.path.join('data', 'filters', 'filter_list.ascii'))
    tbl = Table.read(filter_file, format='ascii')

    allowed_options = tbl['filter'].data

    # Check
    if filter not in allowed_options:
        msgs.error("PypeIt is not ready for filter = {}".format(filter))

    trans_file = resource_filename('pypeit', os.path.join('data', 'filters', 'filtercurves.fits'))
    trans = fits.open(trans_file)
    wave = trans[filter].data['lam']  # Angstroms
    instr = trans[filter].data['Rlam']  # Am keeping in atmospheric terms
    keep = instr > 0.
    # Parse
    wave = wave[keep]
    instr = instr[keep]

    # Return
    return wave, instr

def load_standard_file(std_dict):
    """Load standard star data

    Parameters
    ----------
    std_dict : dict
      Info on standard star indcluding filename in 'file'
      May be compressed

    Returns
    -------
    wave, flux: Quantity, Quantity filled in place in std_dict
      Wavelengths of standard star array
      Flux of standard star in flambda, cgs with scaling of 1e-17
    """
    root = resource_filename('pypeit', std_dict['cal_file'] + '*')
    fil = glob.glob(root)
    if len(fil) == 0:
        msgs.error("No standard star file: {:s}".format(fil))
    else:
        fil = fil[0]
        msgs.info("Loading standard star file: {:s}".format(fil))
        msgs.info("Fluxes are flambda, normalized to 1e-17")

    if std_dict['std_source'] == 'xshooter': # XSHOOTER files
        std_spec = Table.read(fil, format='ascii')
        # Load
        std_dict['wave'] = std_spec['col1'] * units.AA
        std_dict['flux'] = std_spec['col2'] / PYPEIT_FLUX_SCALE  * units.erg / units.s / units.cm ** 2 / units.AA
    elif std_dict['std_source'] == 'calspec': # Calspec
        std_spec = fits.open(fil)[1].data
        # Load
        std_dict['wave'] = std_spec['WAVELENGTH'] * units.AA
        std_dict['flux'] = std_spec['FLUX'] / PYPEIT_FLUX_SCALE * units.erg / units.s / units.cm ** 2 / units.AA
    elif std_dict['std_source'] == 'eso': # ESO files
        std_spec = Table.read(fil, format='ascii')
        # Load
        std_dict['wave'] = std_spec['col1'] * units.AA
        std_dict['flux'] = std_spec['col2'] * units.erg / units.s / units.cm ** 2 / units.AA
    else:
        msgs.error("Bad Standard Star Format")
    return


def find_standard(specobj_list):
    """
    Take the median boxcar and then the max object as the standard

    Parameters
    ----------
    specobj_list : list

    Returns
    -------
    mxix : int
      Index of the standard star

    """
    # Repackage as necessary (some backwards compatability)
    # Do it
    medfx = []
    for indx, spobj in enumerate(specobj_list):
        if spobj is None:
            medfx.append(0.)
        else:
            medfx.append(np.median(spobj.boxcar['COUNTS']))
    try:
        mxix = np.argmax(np.array(medfx))
    except:
        debugger.set_trace()
    msgs.info("Putative standard star {} has a median boxcar count of {}".format(specobj_list[mxix],
                                                                                 np.max(medfx)))
    # Return
    return mxix

def scale_in_filter(xspec, scale_dict):
    """
    Scale spectra to input magnitude in given filter

    scale_dict has data model:
      'filter' (str): name of filter
      'mag' (float): magnitude
      'mag_type' (str, optional): type of magnitude.  Assumed 'AB'
      'masks' (list, optional): Wavelength ranges to mask in calculation

    Args:
        xspec (linetools.spectra.xspectrum1d.XSpectrum1D):
        scale_dict (dict):

    Returns:
        linetools.spectra.xspectrum1d.XSpectrum1D, float:  Scaled spectrum


    """
    # Parse the spectrum
    sig = xspec.sig
    gdx = sig > 0.
    wave = xspec.wavelength.value[gdx]
    flux = xspec.flux.value[gdx]

    # Mask further?
    if 'masks' in scale_dict:
        if scale_dict['masks'] is not None:
            gdp = np.ones_like(wave, dtype=bool)
            for mask in scale_dict['masks']:
                bad = (wave > mask[0]) & (wave < mask[1])
                gdp[bad] = False
            # Cut again
            wave = wave[gdp]
            flux = flux[gdp]

    if ('mag_type' in scale_dict) | (scale_dict['mag_type'] is not None):
        mag_type = scale_dict['mag_type']
    else:
        mag_type = 'AB'

    # Grab the instrument response function
    fwave, trans = load_filter_file(scale_dict['filter'])
    tfunc = scipy.interpolate.interp1d(fwave, trans, bounds_error=False, fill_value=0.)

    # Convolve
    allt = tfunc(wave)
    wflam = np.sum(flux*allt) / np.sum(allt) * PYPEIT_FLUX_SCALE * units.erg/units.s/units.cm**2/units.AA

    mean_wv = np.sum(fwave*trans)/np.sum(trans) * units.AA

    #
    if mag_type == 'AB':
        # Convert flam to AB magnitude
        fnu = wflam * mean_wv**2 / constants.c
        # Apparent AB
        AB = -2.5 * np.log10(fnu.to('erg/s/cm**2/Hz').value) - 48.6
        # Scale factor
        Dm = AB - scale_dict['mag']
        scale = 10**(Dm/2.5)
        msgs.info("Scaling spectrum by {}".format(scale))
        # Generate
        new_spec = XSpectrum1D.from_tuple((xspec.wavelength, xspec.flux*scale, xspec.sig*scale))
    else:
        msgs.error("Need a magnitude for scaling")

    return new_spec,scale

def telluric_params(sptype):
    """Compute physical parameters for a given stellar type.
    This is used by telluric_sed(V, sptype) to create the model spectrum.

    Parameters:
    ----------
    sptype: str
      Spectral type of telluric star

    Returns:
    ----------
    tell_param: dict
      Star parameters
    """

    # log(g) of the Sun
    logg_sol = np.log10(6.67259e-8) + np.log10(1.989e33) - 2.0 * np.log10(6.96e10)

    # Load Schmidt-Kaler (1982) table
    sk82_file = resource_filename('pypeit', 'data/standards/kurucz93/schmidt-kaler_table.txt')
    sk82_tab = ascii.read(sk82_file, names=('Sp', 'logTeff', 'Teff', '(B-V)_0', 'M_V', 'B.C.', 'M_bol', 'L/L_sol'))

    # Match input type
    mti = np.where(sptype == sk82_tab['Sp'])[0]
    if len(mti) != 1:
        raise ValueError('Not ready to interpolate yet.')

    # Calculate final quantities
    # Relation between radius, temp, and bolometric luminosity
    logR = 0.2 * (42.26 - sk82_tab['M_bol'][mti[0]] - 10.0 * sk82_tab['logTeff'][mti[0]])

    # Mass-bolometric luminosity relation from schimdt-kaler p28 valid for M_bol < 7.5
    logM = 0.46 - 0.10 * sk82_tab['M_bol'][mti[0]]
    logg = logM - 2.0 * logR + logg_sol
    M_V = sk82_tab['M_V'][mti[0]]
    tell_param = dict(logR=logR, logM=logM, logg=logg, M_V=M_V,
                      T=sk82_tab['Teff'][mti[0]])

    # Return
    return tell_param


def telluric_sed(V, sptype):
    """Parse Kurucz SED given T and g
    Also convert absolute/apparent magnitudes

    Parameters:
    ----------
    V: float
      Apparent magnitude of the telluric star
    sptype: str
      Spectral type of the telluric star

    Returns:
    ----------
    loglam: ndarray
      log wavelengths
    flux: ndarray
      SED f_lambda (cgs units, I think, probably per Ang)
    """

    # Grab telluric star parameters
    tell_param = telluric_params(sptype)

    # Flux factor (absolute/apparent V mag)

    # Define constants
    parsec = constants.pc.cgs  # 3.086e18
    R_sol = constants.R_sun.cgs  # 6.96e10

    # Distance modulus
    logd = 0.2 * (V - tell_param['M_V']) + 1.0
    D = parsec * 10. ** logd
    R = R_sol * 10. ** tell_param['logR']

    # Factor converts the kurucz surface flux densities to flux observed on Earth
    flux_factor = (R / D.value) ** 2

    # Grab closest T in Kurucz SEDs
    T1 = 3000. + np.arange(28) * 250
    T2 = 10000. + np.arange(6) * 500
    T3 = 13000. + np.arange(22) * 1000
    T4 = 35000. + np.arange(7) * 2500
    Tk = np.concatenate([T1, T2, T3, T4])
    indT = np.argmin(np.abs(Tk - tell_param['T']))

    # Grab closest g in Kurucz SEDs
    loggk = np.arange(11) * 0.5
    indg = np.argmin(np.abs(loggk - tell_param['logg']))

    # Grab Kurucz filename
    std_file = resource_filename('pypeit', '/data/standards/kurucz93/kp00/kp00_{:d}.fits.gz'.format(int(Tk[indT])))
    std = Table.read(std_file)

    # Grab specific spectrum
    loglam = np.array(np.log10(std['WAVELENGTH']))
    gdict = {0: 'g00', 1: 'g05', 2: 'g10', 3: 'g15', 4: 'g20',
             5: 'g25', 6: 'g30', 7: 'g35', 8: 'g40', 9: 'g45',
             10: 'g50'}
    flux = std[gdict[indg]]

    # Generate the standard star dict
    std_dict = dict(stellar_type=sptype, Vmag=V)

    # Return
    return loglam, flux.data * flux_factor, std_dict



# TODO I believe this function is now deprecated.
def apply_sensfunc_spec(wave, counts, ivar, sensfunc, airmass, exptime, mask=None, extinct_correct=True, telluric=None,
                        longitude=None, latitude=None, debug=False):

    if mask is None:
        mask = ivar > 0.0

    # Did the user request a telluric correction from the same file?
    if telluric is not None:
        # This assumes there is a separate telluric key in this dict.
        msgs.info('Applying telluric correction')
        sensfunc = sensfunc*(telluric > 1e-10)/(telluric + (telluric < 1e-10))

    if extinct_correct:
        if longitude is None or latitude is None:
            msgs.error('You must specify longitude and latitude if we are extinction correcting')
        # Apply Extinction if optical bands
        msgs.info("Applying extinction correction")
        msgs.warn("Extinction correction applyed only if the spectra covers <10000Ang.")
        extinct = load_extinction_data(longitude,latitude)
        ext_corr = extinction_correction(wave* units.AA, airmass, extinct)
        senstot = sensfunc * ext_corr
    else:
        senstot = sensfunc.copy()

    flam = counts * senstot/ exptime
    flam_ivar = ivar / (senstot / exptime) **2

    # Mask bad pixels
    msgs.info(" Masking bad pixels")
    outmask =  mask & (senstot>0.)

    # debug
    if debug:
        wave_mask = wave > 1.0
        fig = plt.figure(figsize=(12, 8))
        ymin, ymax = coadd.get_ylim(flam, flam_ivar, outmask)
        plt.plot(wave[wave_mask], flam[wave_mask], color='black', drawstyle='steps-mid', zorder=1, alpha=0.8)
        plt.plot(wave[wave_mask], np.sqrt(utils.calc_ivar(flam_ivar[wave_mask])), zorder=2, color='red', alpha=0.7,
                       drawstyle='steps-mid', linestyle=':')
        plt.ylim([ymin,ymax])
        plt.xlim([wave[wave_mask].min(),wave[wave_mask].max()])
        plt.xlabel('Wavelength (Angstrom)')
        plt.ylabel('Flux')
        plt.show()

    return flam, flam_ivar, outmask

def apply_sensfunc_specobjs(specobjs, sens_meta, sens_table, airmass, exptime, extinct_correct=True, tell_correct=False,
                            longitude=None, latitude=None, debug=False, show=False):

    # TODO This function should operate on a single object
    func = sens_meta['FUNC'][0]
    polyorder_vec = sens_meta['POLYORDER_VEC'][0]
    nimgs = len(specobjs)

    if show:
        fig = plt.figure(figsize=(12, 8))
        xmin, xmax = [], []
        ymin, ymax = [], []

    for ispec in range(nimgs):
        # get the ECH_ORDER, ECH_ORDERINDX, WAVELENGTH from your science
        sobj_ispec = specobjs[ispec]
        ## TODO Comment on the logich here. Hard to follow
        try:
            ech_order, ech_orderindx, idx = sobj_ispec.ech_order, sobj_ispec.ech_orderindx, sobj_ispec.idx
            msgs.info('Applying sensfunc to Echelle data')
        except:
            ech_orderindx = 0
            idx = sobj_ispec.idx
            msgs.info('Applying sensfunc to Longslit/Multislit data')

        # Hotfix for multi-slit data where ech_orderindx is populated with 999
        if ech_orderindx == 999:
            msgs.info('Switch to applying sensfunc to Longslit/Multislit data')
            ech_orderindx = 0
            polyorder_vec = [polyorder_vec]

        for extract_type in ['boxcar', 'optimal']:
            extract = getattr(sobj_ispec, extract_type)

            if len(extract) == 0:
                continue
            msgs.info("Fluxing {:s} extraction for:".format(extract_type) + msgs.newline() + "{}".format(idx))
            wave = extract['WAVE'].value.copy()
            wave_mask = wave > 1.0
            counts = extract['COUNTS'].copy()
            counts_ivar = extract['COUNTS_IVAR'].copy()
            mask = extract['MASK'].copy()

            # get sensfunc from the sens_table
            coeff = sens_table[ech_orderindx]['OBJ_THETA'][0:polyorder_vec[ech_orderindx] + 2]
            wave_min = sens_table[ech_orderindx]['WAVE_MIN']
            wave_max = sens_table[ech_orderindx]['WAVE_MAX']
            sensfunc = np.zeros_like(wave)
            sensfunc[wave_mask] = np.exp(utils.func_val(coeff, wave[wave_mask], func,
                                             minx=wave_min, maxx=wave_max))

            # get telluric from the sens_table
            if tell_correct:
                msgs.work('Evaluate telluric!')
                telluric = None
            else:
                telluric = None

            flam, flam_ivar, outmask = apply_sensfunc_spec(wave, counts, counts_ivar, sensfunc, airmass, exptime,
                                                           mask=mask, extinct_correct=extinct_correct, telluric=telluric,
                                                           longitude=longitude, latitude=latitude, debug=debug)
            flam_sig = np.sqrt(utils.inverse(flam_ivar))
            # The following will be changed directly in the specobjs, so do not need to return anything.
            extract['MASK'] = outmask
            extract['FLAM'] = flam
            extract['FLAM_SIG'] = flam_sig
            extract['FLAM_IVAR'] = flam_ivar

            if show:
                xmin_ispec = wave[wave_mask].min()
                xmax_ispec = wave[wave_mask].max()
                xmin.append(xmin_ispec)
                xmax.append(xmax_ispec)
                ymin_ispec, ymax_ispec = coadd.get_ylim(flam, flam_ivar, outmask)
                ymin.append(ymin_ispec)
                ymax.append(ymax_ispec)

                med_width = (2.0 * np.ceil(0.1 / 10.0 * np.size(wave[outmask])) + 1).astype(int)
                flam_med, flam_ivar_med = coadd.median_filt_spec(flam, flam_ivar, outmask, med_width)
                if extract_type == 'boxcar':
                    plt.plot(wave[wave_mask], flam_med[wave_mask], color='black', drawstyle='steps-mid', zorder=1, alpha=0.8)
                    #plt.plot(wave[wave_mask], np.sqrt(utils.calc_ivar(flam_ivar_med[wave_mask])), zorder=2, color='m',
                    #         alpha=0.7, drawstyle='steps-mid', linestyle=':')
                else:
                    plt.plot(wave[wave_mask], flam_med[wave_mask], color='dodgerblue', drawstyle='steps-mid', zorder=1, alpha=0.8)
                    #plt.plot(wave[wave_mask], np.sqrt(utils.calc_ivar(flam_ivar_med[wave_mask])), zorder=2, color='red',
                    #         alpha=0.7, drawstyle='steps-mid', linestyle=':')
    if show:
        xmin_final, xmax_final = np.min(xmin), np.max(xmax)
        ymax_final = 1.3*np.median(ymax)
        ymin_final = -0.15*ymax_final
        plt.xlim([xmin_final, xmax_final])
        plt.ylim([ymin_final, ymax_final])
        plt.title('Blue is Optimal extraction and Black is Boxcar extraction',fontsize=16)
        plt.xlabel('Wavelength (Angstrom)')
        plt.ylabel('Flux')
        plt.show()



def generate_sensfunc_old(wave, counts, counts_ivar, airmass, exptime, longitude, latitude, telluric=True, star_type=None,
                      star_mag=None, ra=None, dec=None, std_file = None, poly_norder=4, BALM_MASK_WID=5., nresln=20.,
                      resolution=3000.,trans_thresh=0.9,polycorrect=True, debug=False):
    """
    Function to generate the sensitivity function.

    This can work in different regimes:

        - If telluric=False and RA=None and Dec=None the code creates
          a synthetic standard star spectrum (or VEGA spectrum if
          type=A0) using the Kurucz models, and from this it
          generates a sens func using nresln=20.0 and masking out
          telluric regions.

        - If telluric=False and RA and Dec are assigned the standard
          star spectrum is extracted from the archive, and a sens
          func is generated using nresln=20.0 and masking out
          telluric regions.

        - If telluric=True the code creates a sintetic standard star
          spectrum (or VEGA spectrum if type=A0) using the Kurucz
          models, the sens func is a pixelized sensfunc (not smooth)
          for correcting both throughput and telluric lines. if you
          set polycorrect=True, the sensfunc in the Hydrogen
          recombination line region (often seen in star spectra) will
          be replaced by a smoothed polynomial function.

    Parameters
    ----------
    wave : array
        Wavelength of the star [no longer with units]
    counts : array
        Flux (in counts) of the star
    counts_ivar : array
        Inverse variance of the star
    airmass : float
        Airmass
    exptime : float
        Exposure time in seconds
    spectrograph : dict
        Instrument specific dict
        Used for extinction correction
    telluric : bool
        if True performs a telluric correction
    star_type : str
        Spectral type of the telluric star (used if telluric=True)
    star_mag : float
        Apparent magnitude of telluric star (used if telluric=True)
    RA : float
        deg, RA of the telluric star
        if assigned, the standard star spectrum will be extracted from
        the archive
    DEC : float
        deg, DEC of the telluric star
        if assigned, the standard star spectrum will be extracted from
        the archive
    BALM_MASK_WID : float
        Mask parameter for Balmer absorption. A region equal to
        BALM_MASK_WID*resln is masked where resln is the estimate
        for the spectral resolution.
    polycorrect: bool
        Whether you want to correct the sensfunc with polynomial in the Balmer absortion line regions
    poly_norder: int
        Order number of polynomial fit.

    Returns
    -------
    sens_dict : dict
        sensitivity function described by a dict
    """
    # Create copy of the arrays to avoid modification and convert to
    # electrons / s
    wave_star = wave.copy()
    flux_star = counts.copy() / exptime
    ivar_star = counts_ivar.copy() * exptime ** 2

    # Units
    if not isinstance(wave_star, units.Quantity):
        wave_star = wave_star * units.AA

    # Extinction correction
    msgs.info("Applying extinction correction")
    extinct = load_extinction_data(longitude,latitude)
    ext_corr = extinction_correction(wave * units.AA, airmass, extinct)
    # Correct for extinction
    flux_star = flux_star * ext_corr
    ivar_star = ivar_star / ext_corr ** 2

    std_dict = get_standard_spectrum(star_type=star_type, star_mag=star_mag, ra=ra, dec=dec)

    # Interpolate the standard star onto the current set of observed wavelengths
    flux_true = scipy.interpolate.interp1d(std_dict['wave'], std_dict['flux'],bounds_error=False,
                                           fill_value='extrapolate')(wave_star)
    # Do we need to extrapolate? TODO Replace with a model or a grey body?
    if np.min(flux_true) <= 0.:
        msgs.warn('Your spectrum extends beyond calibrated standard star, extrapolating the spectra with polynomial.')
        mask_model = flux_true <= 0
        msk_poly, poly_coeff = utils.robust_polyfit_djs(std_dict['wave'].value, std_dict['flux'].value,8,function='polynomial',
                                                    invvar=None, guesses=None, maxiter=50, inmask=None, \
                                                    lower=3.0, upper=3.0, maxdev=None, maxrej=3, groupdim=None,
                                                    groupsize=None,groupbadpix=False, grow=0, sticky=True, use_mad=True)
        star_poly = utils.func_val(poly_coeff, wave_star.value, 'polynomial')
        #flux_true[mask_model] = star_poly[mask_model]
        flux_true = star_poly.copy()
        if debug:
            plt.plot(std_dict['wave'], std_dict['flux'],'bo',label='Raw Star Model')
            plt.plot(std_dict['wave'],  utils.func_val(poly_coeff, std_dict['wave'].value, 'polynomial'), 'k-',label='robust_poly_fit')
            plt.plot(wave_star,flux_true,'r-',label='Your Final Star Model used for sensfunc')
            plt.show()

    # Get masks from observed star spectrum. True = Good pixels
    msk_bad, msk_star, msk_tell = get_mask(wave_star.value, flux_star, ivar_star, mask_balmer=True, mask_tell=True,
                                           BALM_MASK_WID=BALM_MASK_WID, trans_thresh=0.9)

    # Get sensfunc
    LBLRTM = False
    if LBLRTM:
        # sensfunc = lblrtm_sensfunc() ???
        msgs.develop('fluxing and telluric correction based on LBLRTM model is under developing.')
    else:
        sensfunc, mask_sens = standard_sensfunc(wave_star.value, flux_star, ivar_star, flux_true, msk_bad=msk_bad, msk_star=msk_star,
                                msk_tell=msk_tell, maxiter=35,upper=3.0, lower=3.0, poly_norder=poly_norder,
                                BALM_MASK_WID=BALM_MASK_WID,nresln=nresln,telluric=telluric, resolution=resolution,
                                polycorrect= polycorrect, debug=debug, show_QA=False)

    if debug:
        plt.plot(wave_star.value[mask_sens], flux_true[mask_sens], color='k',lw=2,label='Reference Star')
        plt.plot(wave_star.value[mask_sens], flux_star[mask_sens]*sensfunc[mask_sens], color='r',label='Fluxed Observed Star')
        plt.xlabel(r'Wavelength [$\AA$]')
        plt.ylabel('Flux [erg/s/cm2/Ang.]')
        plt.legend(fancybox=True, shadow=True)
        plt.show()


    # Add in wavemin,wavemax
    sens_dict = {}
    sens_dict['wave'] = wave_star.value
    sens_dict['sensfunc'] = sensfunc
    sens_dict['wave_min'] = np.min(wave_star)
    sens_dict['wave_max'] = np.max(wave_star)
    sens_dict['exptime']= exptime
    sens_dict['airmass']= airmass
    sens_dict['std_file']= std_file
    # Get other keys from standard dict
    sens_dict['std_ra'] = std_dict['std_ra']
    sens_dict['std_dec'] = std_dict['std_dec']
    sens_dict['std_name'] = std_dict['name']
    sens_dict['cal_file'] = std_dict['cal_file']
    sens_dict['flux_true'] = flux_true
    sens_dict['mask_sens'] = mask_sens
    #sens_dict['std_dict'] = std_dict
    #sens_dict['msk_star'] = msk_star
    #sens_dict['mag_set'] = mag_set

    return sens_dict