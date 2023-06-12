from pypeit.core.wavecal import templates
import glob

wavelength_soln_files = glob.glob("./wvcalib*.fits")

wavelength_soln_files = [
    "wvcalib_order15.fits", # [4800.0,     4100.0]
    "wvcalib_order14.fits", # [5100.0,     4300.0],
    "wvcalib_order13.fits", # [5500.0,     4500.0],
    "wvcalib_order12.fits", # [6100.0,     4900.0],
    "wvcalib_order11.fits", # [6500.0,     5300.0],
    "wvcalib_order10.fits", # [7300.0,     5900.0],
    "wvcalib_order9.fits", #   [8100.0,     6500.0],
    "wvcalib_order8.fits" #     [7300.0,     9100.0],

                         ]

slits = [
         #49,188,313,429,538,641,739,836 # THIS FOR UNMASKED
         #82,198,323,439,548,651,749,845 # THIS FOR MASKED
         182,198,313,429,538,641,739,836  # This is a blend
         
         ]

wv_cuts = [
    4400, 4950, 5300, 5800, 6300, 7000, 7800
    ]

binspec = 1
outroot = "not_nte_vis.fits"
#print(wavelength_soln_files)

templates.build_template(wavelength_soln_files,
                         slits,
                         wv_cuts,
                         binspec,
                         outroot,
                         #ifiles=ifiles,
                         #det_cut=det_cut,
                         chk=True,
                         normalize=False,
                         lowredux=False,
                         subtract_conti=True,
                         #overwrite=overwrite,
                         shift_wave=True,
                         in_vac=False)
