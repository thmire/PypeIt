import pytest

from IPython import embed

import numpy

from astropy.io import fits

from pypeit.bitmask import BitMask
from pypeit.slittrace import SlitTraceBitMask

#-----------------------------------------------------------------------------

class ImageBitMask(BitMask):
    def __init__(self):
        bits = {'BPM': 'Pixel is part of a bad-pixel mask',
                'COSMIC': 'Pixel is contaminated by a cosmic ray',
                'SATURATED': 'Pixel is saturated.'}
        super(ImageBitMask, self).__init__(list(bits.keys()), descr=list(bits.values()))


def test_new():
    image_bm = ImageBitMask()
    assert list(image_bm.bits.keys()) == ['BPM', 'COSMIC', 'SATURATED']
    assert image_bm.keys() == ['BPM', 'COSMIC', 'SATURATED']
    assert list(image_bm.bits.values()) == [0, 1, 2]


def test_flagging():

    n = 1024
    shape = (n,n)

    image_bm = ImageBitMask()
    mask = numpy.zeros(shape, dtype=image_bm.minimum_dtype())

    cosmics_indx = numpy.zeros(shape, dtype=bool)
    cosmics_indx[numpy.random.randint(0,high=n,size=100),
                 numpy.random.randint(0,high=n,size=100)] = True
    mask[cosmics_indx] = image_bm.turn_on(mask[cosmics_indx], 'COSMIC')

    saturated_indx = numpy.zeros(shape, dtype=bool)
    saturated_indx[numpy.random.randint(0,high=n,size=100),
                   numpy.random.randint(0,high=n,size=100)] = True
    mask[saturated_indx] = image_bm.turn_on(mask[saturated_indx], 'SATURATED')

    assert numpy.sum(image_bm.flagged(mask, flag='BPM')) == 0
    assert numpy.sum(image_bm.flagged(mask, flag='COSMIC')) == numpy.sum(cosmics_indx)
    assert numpy.sum(image_bm.flagged(mask, flag='SATURATED')) == numpy.sum(saturated_indx)

    assert image_bm.flagged_bits(1) == ['BPM']
    assert image_bm.flagged_bits(2) == ['COSMIC']
    assert image_bm.flagged_bits(4) == ['SATURATED']

    unique_flags = numpy.sort(numpy.unique(numpy.concatenate(
                        [image_bm.flagged_bits(b) for b in numpy.unique(mask)]))).tolist()
    assert unique_flags == ['COSMIC', 'SATURATED']

    mask[saturated_indx] = image_bm.turn_off(mask[saturated_indx], 'SATURATED')
    assert numpy.sum(image_bm.flagged(mask, flag='COSMIC')) == numpy.sum(cosmics_indx)
    assert numpy.sum(image_bm.flagged(mask, flag='SATURATED')) == 0

    unique_flags = numpy.sort(numpy.unique(numpy.concatenate(
                        [image_bm.flagged_bits(b) for b in numpy.unique(mask)]))).tolist()
    assert unique_flags == ['COSMIC']

    b_indx, c_indx, s_indx = image_bm.unpack(mask)
    assert numpy.sum(b_indx) == 0
    assert numpy.sum(c_indx) == numpy.sum(cosmics_indx)
    assert numpy.sum(s_indx) == 0


def test_hdr_io():
    
    image_bm = ImageBitMask()
    hdr = fits.Header()
    image_bm.to_header(hdr)

    assert list(hdr.keys()) == ['BIT0', 'BIT1', 'BIT2']
    assert list(hdr.values()) == ['BPM', 'COSMIC', 'SATURATED']

    bm = BitMask.from_header(hdr)
    assert list(bm.bits.keys()) == ['BPM', 'COSMIC', 'SATURATED']
    assert list(bm.bits.values()) == [0, 1, 2]


def test_wrong_bits():
    image_bm = ImageBitMask()

    n = 1024
    shape = (n,n)

    image_bm = ImageBitMask()
    mask = numpy.zeros(shape, dtype=image_bm.minimum_dtype())

    cosmics_indx = numpy.zeros(shape, dtype=bool)
    cosmics_indx[numpy.random.randint(0,high=n,size=100),
                 numpy.random.randint(0,high=n,size=100)] = True
    mask[cosmics_indx] = image_bm.turn_on(mask[cosmics_indx], 'COSMIC')
    
    # Fails with all bad flags
    with pytest.raises(ValueError):
        out = image_bm.flagged(mask, flag='JUNK')

    # Fails with mix of good and bad flags
    with pytest.raises(ValueError):
        out = image_bm.flagged(mask, flag=['COSMIC', 'JUNK'])

    assert numpy.sum(image_bm.flagged(mask, flag='COSMIC')) == numpy.sum(cosmics_indx)

def test_flag_order():

    bm = ImageBitMask()

    flags = bm.keys()
    assert bm.correct_flag_order(flags), 'Flags should not be mismatched'

    flags += ['NEWBIT']
    assert bm.correct_flag_order(flags), 'Appending flags should be fine'

    flags = bm.keys()[:-1]
    assert bm.correct_flag_order(flags), 'Checking a subset of the flags should be fine'

    flags = bm.keys()[::-1]
    assert not bm.correct_flag_order(flags), 'Reordering the flags is not okay'


def test_exclude_and_not():
    
    n = 1024
    shape = (n,n)

    rng = numpy.random.default_rng(99)

    image_bm = ImageBitMask()
    mask = numpy.zeros(shape, dtype=image_bm.minimum_dtype())

    cosmics_indx = numpy.zeros(shape, dtype=bool)
    cosmics_indx[rng.integers(0,high=n,size=9000), rng.integers(0,high=n,size=9000)] = True
    mask[cosmics_indx] = image_bm.turn_on(mask[cosmics_indx], 'COSMIC')

    saturated_indx = numpy.zeros(shape, dtype=bool)
    saturated_indx[rng.integers(0,high=n,size=9000), rng.integers(0,high=n,size=9000)] = True
    mask[saturated_indx] = image_bm.turn_on(mask[saturated_indx], 'SATURATED')

    # NOTE: Want to make sure there are pixels flagged as both COSMIC and
    # SATURATED.  Otherwise the `and_not` test is not useful.
    assert numpy.sum(cosmics_indx & saturated_indx) > 0, 'Bad test setup'

    assert numpy.array_equal(image_bm.flagged(mask), cosmics_indx | saturated_indx), \
            'Mask incorrect'
    assert numpy.array_equal(image_bm.flagged(mask, exclude='SATURATED'), cosmics_indx), \
            'Exclude incorrect'

    assert numpy.array_equal(image_bm.flagged(mask, and_not='SATURATED'),
                             cosmics_indx & numpy.logical_not(saturated_indx)), 'Expunge incorrect'

def test_boxslit():
    """
    Tests old vs. new bpm after adding `and_not` functionality.
    """
    bm = SlitTraceBitMask()
    v = numpy.array([10,0,2])

    desired_bpm = (v > 2) & (numpy.invert(bm.flagged(v, flag=bm.exclude_for_reducing)))

    new_bpm = bm.flagged(v, exclude='BOXSLIT', and_not=bm.exclude_for_reducing)

    assert numpy.all(new_bpm == desired_bpm)


