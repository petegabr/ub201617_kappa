from unittest import TestCase

from Bio.SeqRecord import SeqRecord

from main import get_binding_sites
from parse import TDPBindingSite


class TestBindingSites(TestCase):
    @staticmethod
    def _build_site(start, end, direction):
        if direction == TDPBindingSite.DIRECTION_POSITIVE:
            sign = '+'
        else:
            sign = '-'
        return TDPBindingSite.from_string(
            '1	%s	%s	na	0	%s' % (start, end, sign))

    @staticmethod
    def _build_sequence(sequence, direction, start=0):
        seq = SeqRecord(sequence)
        seq.start = start
        seq.end = len(sequence) - 1
        seq.direction = direction
        return seq

    def test_binding_site_contained_within(self):
        sequence = self._build_sequence(
            'MARYHADALITTLELAMBITSFLEECEASWHITEASSNOW',
            TDPBindingSite.DIRECTION_POSITIVE)

        binding_site = self._build_site(4, 6, TDPBindingSite.DIRECTION_POSITIVE)
        self.assertTrue(
            binding_site.is_contained_within(sequence),
            'Valid binding site detected as invalid')

        binding_site = self._build_site(100, 102, TDPBindingSite.DIRECTION_POSITIVE)
        self.assertFalse(
            binding_site.is_contained_within(sequence),
            'Binding site not in gene range not detected as invalid')

        binding_site = self._build_site(4, 6, TDPBindingSite.DIRECTION_NEGATIVE)
        self.assertFalse(
            binding_site.is_contained_within(sequence),
            'Binding site with wrong direction not detected as invalid')

    def test_get_binding_sites(self):
        direction = TDPBindingSite.DIRECTION_POSITIVE
        sequence = self._build_sequence(
            'MARYHADALITTLELAMBITSFLEECEASWHITEASSNOW', direction, start=20)
        binding_site = self._build_site(24, 27, direction)

        self.assertEquals(get_binding_sites([binding_site], [sequence]), ['HAD'])
