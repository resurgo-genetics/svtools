from unittest import TestCase, main
import tempfile
import os
import sys
import difflib
import svtools.lmerge as lmerge

class LmergeUnitTest(TestCase):
    def test_null_format_string(self):
        self.assertEqual(lmerge.null_format_string('GT:GQ:AD'), './.:.:.')
        self.assertEqual(lmerge.null_format_string('GQ:AD'), '.:.')

    def test_clip_out_tag(self):
        self.assertEqual(lmerge.clip_out_tag('SOMETAG=T;OTHERTAG=G', 'SOMETAG='), 'OTHERTAG=G')
        self.assertEqual(lmerge.clip_out_tag('SOMETAG=T;OTHERTAG=G', 'OTHERTAG='), 'SOMETAG=T')
        self.assertEqual(lmerge.clip_out_tag('SOMETAG=T;OTHERTAG=G', 'MISSINGTAG='), 'SOMETAG=T;OTHERTAG=G')

    def test_update_sname(self):
        self.assertEqual(lmerge.update_sname('SNAME=NAME;OTHERTAG', '12'), ('SNAME=NAME:12;OTHERTAG', 'NAME'))
        self.assertEqual(lmerge.update_sname('OTHERTAG;SNAME=NAME', '12'), ('OTHERTAG;SNAME=NAME:12', 'NAME'))

class LmergeIntegrationTest(TestCase):
    def run_integration_test(self):
        test_directory = os.path.dirname(os.path.abspath(__file__))
        self.test_data_dir = os.path.join(test_directory, 'test_data', 'lmerge')
        input_file = os.path.join(self.test_data_dir, 'input')
        expected_result = os.path.join(self.test_data_dir, 'expected.vcf')
        temp_descriptor, temp_output_path = tempfile.mkstemp(suffix='.vcf')
        with os.fdopen(temp_descriptor, 'w') as output_handle:
            # FIXME this is pretty hacky
            temp_handle, sys.stdout = sys.stdout, output_handle

            lmerge.l_cluster_by_line(input_file, 0.0, 20, True)
            output_handle.flush()
            expected_lines = open(expected_result).readlines()
            produced_lines = open(temp_output_path).readlines()
            diff = difflib.unified_diff(produced_lines, expected_lines, fromfile=temp_output_path, tofile=expected_result)
            result = ''.join(diff)
            if result != '':
                for line in result:
                    sys.stdout.write(line)
                self.assertFalse(result)
        os.remove(temp_output_path)

if __name__ == "__main__":
    main()
