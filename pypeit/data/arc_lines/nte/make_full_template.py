from pypeit.core.wavecal import templates
import glob

wavelength_soln_files = glob.glob("./*.fits")

templates.build_template("./wvarxiv_not_nte_vis_order8_20230418T1401.fits",
                         slits,
                         wv_cuts,
                         binspec,
                         outroot,
                         ifiles=ifiles,
                         det_cut=det_cut, chk=True, normalize=False, lowredux=False,
                         subtract_conti=True, overwrite=overwrite, shift_wave=True)
