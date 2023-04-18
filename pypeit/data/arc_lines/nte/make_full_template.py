from pypeit.core.wavecal import templates
import glob

wavelength_soln_files = glob.glob("./wvcalib*.fits")

wavelength_soln_files = [
    "wvcalib_order8.fits","wvcalib_order9.fits"
    "wvcalib_order10.fits","wvcalib_order11.fits"
    "wvcalib_order12.fits","wvcalib_order13.fits"
    "wvcalib_order14.fits","wvcalib_order15.fits"

                         ]

print(wavelength_soln_files)

templates.build_template(wavelength_soln_files,
                         [],
                         wv_cuts,
                         binspec,
                         outroot,
                         ifiles=ifiles,
                         det_cut=det_cut, chk=True, normalize=False, lowredux=False,
                         subtract_conti=True, overwrite=overwrite, shift_wave=True)
