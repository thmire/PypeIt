"""
Module to run tests on SlitTraceSet
"""
from pathlib import Path

from IPython import embed

import numpy as np

from pypeit.slittrace import SlitTraceSet, SlitTraceBitMask
from pypeit.tests.tstutils import data_output_path


def test_bits():
    # Make sure bits are correct
    bm = SlitTraceBitMask()
    assert bm.bits['SHORTSLIT'] == 0, 'Bits changed'
    assert bm.bits['BOXSLIT'] == 1, 'Bits changed'
    assert bm.bits['USERIGNORE'] == 2, 'Bits changed'
    assert bm.bits['BADWVCALIB'] == 3, 'Bits changed'
    assert bm.bits['BADTILTCALIB'] == 4, 'Bits changed'
    assert bm.bits['BADALIGNCALIB'] == 5, 'Bits changed'
    assert bm.bits['SKIPFLATCALIB'] == 6, 'Bits changed'
    assert bm.bits['BADFLATCALIB'] == 7, 'Bits changed'
    assert bm.bits['BADREDUCE'] == 8, 'Bits changed'
    assert bm.bits['BADSKYSUB'] == 9, 'Bits changed'
    assert bm.bits['BADEXTRACT'] == 10, 'Bits changed'


def test_init():
    slits = SlitTraceSet(left_init=np.full((1000,3), 2, dtype=float),
                         right_init=np.full((1000,3), 8, dtype=float),
                         pypeline='MultiSlit',
                         nspat=10, PYP_SPEC='dummy')

    left, right, _ = slits.select_edges()
    center = (left+right)/2
    assert np.all(center == 5), 'Bad center'


def test_io():

    slits = SlitTraceSet(np.full((1000,3), 2, dtype=float), np.full((1000,3), 8, dtype=float),
                         'MultiSlit', nspat=10, PYP_SPEC='dummy')
    slits.set_paths(data_output_path(''), 'A', '1', 'DET01')
    ofile = Path(slits.get_path()).absolute()

    # Try to save it
    slits.to_file(overwrite=True)
    assert ofile.exists(), 'File not written'

    # Instantiate an empty SlitTraceSet with the same file, and
    # indicate it should be reused
    _slits = SlitTraceSet.from_file(ofile)
    assert np.array_equal(_slits.left_init, np.full((1000,3), 2, dtype=float)), 'Bad left read'
    # And that it's the same as the existing one
    assert np.array_equal(_slits.left_init, slits.left_init), 'Bad left read'

    ofile.unlink()


def test_io_single():
    # NOTE this is just a test string. The file itself is not actually read, so it is not required to run the test.
    file = '/home/xavier/Projects/PypeIt-development-suite/RAW_DATA/keck_deimos/830G_M_8500/DE.20100913.57006.fits.gz'
    slits = SlitTraceSet(np.full((1000, 1), 2, dtype=float), np.full((1000, 1), 8, dtype=float),
                         'MultiSlit',
                         nspat=10, PYP_SPEC='dummy',
                         maskfile=file)
    slits.set_paths(data_output_path(''), 'A', '1', 'DET01')
    ofile = Path(slits.get_path()).absolute()

    # Try to save it
    slits.to_file()

    # And read it back in
    _slits = SlitTraceSet.from_file(ofile)

    assert np.array_equal(_slits.left_init, np.full((1000, 1), 2, dtype=float)), 'Bad left read'
    # And that it's the same as the existing one
    assert np.array_equal(_slits.left_init, slits.left_init), 'Bad left read'

    ofile.unlink()


