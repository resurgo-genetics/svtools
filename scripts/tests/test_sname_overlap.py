from unittest import TestCase
import sname_overlap as so

class OverlapTests(TestCase):
    def test_set_from_string(self):
        self.assertEqual(set(), so.set_from_string(''))
        self.assertEqual(set(['HI']), so.set_from_string('HI'))
        self.assertEqual(set(['25', '36']), so.set_from_string('25,36'))

    def test_overlapping_ids(self):
        filter_list = [ (0, set([1, 2])), (1, set([1])), (2, set([3])) ]
        self.assertEqual([], so.overlapping_ids(set([0]), filter_list))
        self.assertEqual([0, 1], so.overlapping_ids(set([1]), filter_list))
        self.assertEqual([0], so.overlapping_ids(set([2]), filter_list))
        self.assertEqual([2], so.overlapping_ids(set([3]), filter_list))


