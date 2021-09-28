import unittest

from projman_filler.models.db_models import SampleResult
from projman_filler.sample_level_statistics import calculate_sample_statistics
from projman_filler.bcl2fastq_run_stats_parser import Bcl2fastqRunStatsParser

from tests.test_utils import *
from unittest.mock import MagicMock


class TestSampleLevelStatistics(unittest.TestCase):
    flowcell_id = "foo"
    reads_and_cycles = {1: 151, 2: 151}
    error_rates = {1: {1: 2.0, 2: 2.0}, 2: {1: 2.0, 2: 2.0}}
    densities = {
        1:
            {1: {"raw_density": 100000, "pass_filter_density": 100000},
             2: {"raw_density": 100000, "pass_filter_density": 100000}},
        2:
            {1: {"raw_density": 100000, "pass_filter_density": 100000},
             2: {"raw_density": 100000, "pass_filter_density": 100000}},
    }
    q30s = {1: {1: 95.0, 2: 92.0}, 2: {1: 96.1, 2: 94.2}}

    class SampleSheetMock(object):
        def __init__(self):
            self.project_dict = {"Sample_A": "Project1",
                                 "Sample_B": "Project2",
                                 "Sample_C": "Project1",
                                 "Sample_D": "Project2"}

        def project_for_sample(self, sample_id, lane):
            return self.project_dict[sample_id]

        def library_name_for_sample(self, sample_id, lane):
            return "{}.library".format(sample_id)

    samplesheet_mock = SampleSheetMock()

    def preprocess_conversion_results(self, cr):
        Bcl2fastqRunStatsParser.__init__ = MagicMock(return_value=None)
        stats_parser = Bcl2fastqRunStatsParser("foo")
        stats_parser._stats = {'ConversionResults': cr}
        return stats_parser.get_conversion_results()

    def test_calculate_sample_level_statistics(self):
        cr = self.preprocess_conversion_results(conversion_results)
        actual = list(calculate_sample_statistics(flowcell_name=self.flowcell_id,
                                                  conversion_results=cr,
                                                  reads_and_cycles=self.reads_and_cycles,
                                                  samplesheet=self.samplesheet_mock))

        # One row for per sample, index (in the test data there are 3 per sample), and 2 reads.
        self.assertEqual(len(actual), 4*3*2)

        actual_sample_a = list(filter(lambda x: x.sample_name == 'A', actual))
        list_of_values_for_a = [
            {'flowcell_id': 'foo', 'project_id': 'Project1', 'sample_name': 'A',
             'tag_seq': 'GTAGAGGA-CTCTCTAT', 'lane_num': 1, 'read_num': 1, 'cycles': 151,
             'pct_lane': 49.91040361971908, 'pf_clusters': 81217423.0,
             'pct_q30': 98.02429332249935, 'pct_tag_err': 0.671112157794024,
             'library_name': 'Sample_A.library', 'mean_q': 38.84148990743496},
            {'flowcell_id': 'foo', 'project_id': 'Project1', 'sample_name': 'A',
             'tag_seq': 'GTAGAGGA-CTCTCTAT', 'lane_num': 1, 'read_num': 2, 'cycles': 151,
             'pct_lane': 49.91040361971908, 'pf_clusters': 81217423.0,
             'pct_q30': 96.45192508767363, 'pct_tag_err': 0.671112157794024,
             'library_name': 'Sample_A.library', 'mean_q': 38.373262536376345},
            {'flowcell_id': 'foo', 'project_id': 'Project1', 'sample_name': 'A',
             'tag_seq': 'TAGGCATG-CTCTCTAT', 'lane_num': 1, 'read_num': 1, 'cycles': 151,
             'pct_lane': 49.91040361971908, 'pf_clusters': 81217423.0,
             'pct_q30': 98.02429332249935, 'pct_tag_err': 0.7880181078880083,
             'library_name': 'Sample_A.library', 'mean_q': 38.84148990743496},
            {'flowcell_id': 'foo', 'project_id': 'Project1', 'sample_name': 'A',
             'tag_seq': 'TAGGCATG-CTCTCTAT', 'lane_num': 1, 'read_num': 2, 'cycles': 151,
             'pct_lane': 49.91040361971908, 'pf_clusters': 81217423.0,
             'pct_q30': 96.45192508767363, 'pct_tag_err': 0.7880181078880083,
             'library_name': 'Sample_A.library', 'mean_q': 38.373262536376345},
            {'flowcell_id': 'foo', 'project_id': 'Project1', 'sample_name': 'A',
             'tag_seq': 'TCCTGAGC-CTCTCTAT', 'lane_num': 1, 'read_num': 1, 'cycles': 151,
             'pct_lane': 49.91040361971908, 'pf_clusters': 81217423.0,
             'pct_q30': 98.02429332249935, 'pct_tag_err': 0.7687463809335591,
             'library_name': 'Sample_A.library', 'mean_q': 38.84148990743496},
            {'flowcell_id': 'foo', 'project_id': 'Project1', 'sample_name': 'A',
             'tag_seq': 'TCCTGAGC-CTCTCTAT', 'lane_num': 1, 'read_num': 2, 'cycles': 151,
             'pct_lane': 49.91040361971908, 'pf_clusters': 81217423.0,
             'pct_q30': 96.45192508767363, 'pct_tag_err': 0.7687463809335591,
             'library_name': 'Sample_A.library', 'mean_q': 38.373262536376345}]

        expected_sample_a = list(map(lambda x: SampleResult(**x), list_of_values_for_a))
        self.assertListEqual(expected_sample_a, actual_sample_a)

    def test_calculate_sample_level_statistics_without_index_metrics(self):
        cr = self.preprocess_conversion_results(conversion_results_without_index_metrics)
        actual = list(calculate_sample_statistics(flowcell_name=self.flowcell_id,
                                                  conversion_results=cr,
                                                  reads_and_cycles=self.reads_and_cycles,
                                                  samplesheet=self.samplesheet_mock))

        # One row for per sample, lane, index (in this case only one 'unknown'), and 2 reads.
        self.assertEqual(len(actual), 4)

        actual_sample_a = list(filter(lambda x: x.sample_name == 'A', actual))
        list_of_values_for_a = [
            {'flowcell_id': 'foo', 'project_id': 'Project1', 'sample_name': 'A',
             'tag_seq': 'unknown', 'lane_num': 1, 'read_num': 1, 'cycles': 151,
             'pct_lane': 49.91040361971908, 'pf_clusters': 81217423.0,
             'pct_q30': 98.02429332249935, 'pct_tag_err': None,
             'library_name': 'Sample_A.library', 'mean_q': 38.84148990743496},
            {'flowcell_id': 'foo', 'project_id': 'Project1', 'sample_name': 'A',
             'tag_seq': 'unknown', 'lane_num': 1, 'read_num': 2, 'cycles': 151,
             'pct_lane': 49.91040361971908, 'pf_clusters': 81217423.0,
             'pct_q30': 96.45192508767363, 'pct_tag_err': None,
             'library_name': 'Sample_A.library', 'mean_q': 38.373262536376345}]

        expected_sample_a = list(map(lambda x: SampleResult(**x), list_of_values_for_a))
        self.assertListEqual(expected_sample_a, actual_sample_a)

    def test_calculate_sample_level_statistics_sample_with_no_reads(self):
        cr = self.preprocess_conversion_results(conversion_results_sample_with_no_reads)
        actual = list(calculate_sample_statistics(flowcell_name=self.flowcell_id,
                                                  conversion_results=cr,
                                                  reads_and_cycles=self.reads_and_cycles,
                                                  samplesheet=self.samplesheet_mock))

        # One row per sample, lane , index and read
        # In this case: 2 lanes, 2 reads, 3 samples with 3 indices and 1 with 1 index
        self.assertEqual(len(actual), 3*3*2+1*1*2)

        actual_sample_a = list(filter(lambda x: x.sample_name == 'A', actual))
        list_of_values_for_a = [
            {'flowcell_id': 'foo', 'project_id': 'Project1', 'sample_name': 'A',
             'tag_seq': 'GTAGAGGA-CTCTCTAT', 'lane_num': 1, 'read_num': 1, 'cycles': 151,
             'pct_lane': 0, 'pf_clusters': 0,
             'pct_q30': None, 'pct_tag_err': None,
             'library_name': 'Sample_A.library', 'mean_q': None},
            {'flowcell_id': 'foo', 'project_id': 'Project1', 'sample_name': 'A',
             'tag_seq': 'GTAGAGGA-CTCTCTAT', 'lane_num': 1, 'read_num': 2, 'cycles': 151,
             'pct_lane': 0, 'pf_clusters': 0,
             'pct_q30': None, 'pct_tag_err': None,
             'library_name': 'Sample_A.library', 'mean_q': None}]

        expected_sample_a = list(map(lambda x: SampleResult(**x), list_of_values_for_a))
        self.assertListEqual(expected_sample_a, actual_sample_a)

    def test_calculate_sample_level_statistics_samples_with_multiple_sample_ids(self):
        #For this test we need to update the list of Sample IDs
        self.samplesheet_mock_multiple_sampleIDs = self.SampleSheetMock()
        self.samplesheet_mock_multiple_sampleIDs.project_dict = {"SI-GA-D1_1": "Project1",
                                                                 "SI-GA-D1_2": "Project1",
                                                                 "SI-GA-D1_3": "Project1",
                                                                 "SI-GA-D1_4": "Project1",
                                                                 "SI-GA-F2_1": "Project2",
                                                                 "SI-GA-F2_2": "Project2",
                                                                 "SI-GA-F2_3": "Project2",
                                                                 "SI-GA-F2_4": "Project2",
                                                                 "SI-GA-E1_1": "Project3",
                                                                 "SI-GA-E1_2": "Project3",
                                                                 "SI-GA-E1_3": "Project3",
                                                                 "SI-GA-E1_4": "Project3",
                                                                 "SI-GA-F1_1": "Project4",
                                                                 "SI-GA-F1_2": "Project4",
                                                                 "SI-GA-F1_3": "Project4",
                                                                 "SI-GA-F1_4": "Project4"}

        cr = self.preprocess_conversion_results(conversion_results_multiple_sampleIDs_per_sampleName)
        actual = list(calculate_sample_statistics(flowcell_name=self.flowcell_id,
                                                  conversion_results=cr,
                                                  reads_and_cycles=self.reads_and_cycles,
                                                  samplesheet=self.samplesheet_mock_multiple_sampleIDs))

        # One row per sample, lane, index and read
        self.assertEqual(len(actual), 4*1*4*2)

        actual_sample_a = list(filter(lambda x: x.sample_name == 'A', actual))

        list_of_values_for_a = [
            {'flowcell_id': 'foo', 'project_id': 'Project1', 'sample_name': 'A',
             'lane_num': 1, 'pct_lane': 6.115905919800121, 'library_name': 'SI-GA-D1_1.library',
             'tag_seq': 'CACTCGGA', 'pct_tag_err': 3.621424494761074, 'read_num': 1, 'cycles': 151,
             'mean_q': 39.16345840568938, 'pct_q30': 92.45344394187265, 'pf_clusters': 2360121.0},
            {'flowcell_id': 'foo', 'project_id': 'Project1', 'sample_name': 'A',
             'lane_num': 1, 'pct_lane': 6.115905919800121, 'library_name': 'SI-GA-D1_1.library',
             'tag_seq': 'CACTCGGA', 'pct_tag_err': 3.621424494761074, 'read_num': 2, 'cycles': 151,
             'mean_q': 36.559150347300495, 'pct_q30': 84.04189940076341, 'pf_clusters': 2360121.0},
            {'flowcell_id': 'foo', 'project_id': 'Project1', 'sample_name': 'A',
             'lane_num': 1, 'pct_lane': 5.757983102514638, 'library_name': 'SI-GA-D1_2.library',
             'tag_seq': 'GCTGAATT', 'pct_tag_err': 3.720388713046226, 'read_num': 1, 'cycles': 151,
             'mean_q': 39.23153569054157, 'pct_q30': 92.65040530937226, 'pf_clusters': 2221999.0},
            {'flowcell_id': 'foo', 'project_id': 'Project1', 'sample_name': 'A',
             'lane_num': 1, 'pct_lane': 5.757983102514638, 'library_name': 'SI-GA-D1_2.library',
             'tag_seq': 'GCTGAATT', 'pct_tag_err': 3.720388713046226, 'read_num': 2, 'cycles': 151,
             'mean_q': 36.748773746773196, 'pct_q30': 84.63309285113054, 'pf_clusters': 2221999.0},
            {'flowcell_id': 'foo', 'project_id': 'Project1', 'sample_name': 'A',
             'lane_num': 1, 'pct_lane': 5.187123651369359, 'library_name': 'SI-GA-D1_3.library',
             'tag_seq': 'TGAAGTAC', 'pct_tag_err': 3.5747525234737383, 'read_num': 1, 'cycles': 151,
             'mean_q': 39.18503968413285, 'pct_q30': 92.50719688617771, 'pf_clusters': 2001705.0},
            {'flowcell_id': 'foo', 'project_id': 'Project1', 'sample_name': 'A',
             'lane_num': 1, 'pct_lane': 5.187123651369359, 'library_name': 'SI-GA-D1_3.library',
             'tag_seq': 'TGAAGTAC', 'pct_tag_err': 3.5747525234737383, 'read_num': 2, 'cycles': 151,
             'mean_q': 36.757323706043906, 'pct_q30': 84.6543801541636, 'pf_clusters': 2001705.0},
            {'flowcell_id': 'foo', 'project_id': 'Project1', 'sample_name': 'A',
             'lane_num': 1, 'pct_lane': 6.201171788958993, 'library_name': 'SI-GA-D1_4.library',
             'tag_seq': 'ATGCTCCG', 'pct_tag_err': 3.622402607578274, 'read_num': 1, 'cycles': 151,
             'mean_q': 39.17663714585525, 'pct_q30': 92.50440798869728, 'pf_clusters': 2393025.0},
            {'flowcell_id': 'foo', 'project_id': 'Project1', 'sample_name': 'A',
             'lane_num': 1, 'pct_lane': 6.201171788958993, 'library_name': 'SI-GA-D1_4.library',
             'tag_seq': 'ATGCTCCG', 'pct_tag_err': 3.622402607578274, 'read_num': 2, 'cycles': 151,
              'mean_q': 36.580552589683414, 'pct_q30': 84.11983530225224, 'pf_clusters': 2393025.0}]

        expected_sample_a = list(map(lambda x: SampleResult(**x), list_of_values_for_a))
        self.assertListEqual(expected_sample_a, actual_sample_a)


    def test_calculate_sample_level_statistics_empty_lane(self):
        cr = self.preprocess_conversion_results(conversion_results_empty_lane)
        actual = list(calculate_sample_statistics(flowcell_name=self.flowcell_id,
                                                  conversion_results=cr,
                                                  reads_and_cycles=self.reads_and_cycles,
                                                  samplesheet=self.samplesheet_mock))

        # One row per sample, lane , index and read
        # In this case: 2 lanes, 2 reads, 3 samples with 3 indices and 1 with 1 index
        self.assertEqual(len(actual), 3*3*2+1*1*2)

        actual_sample_a = list(filter(lambda x: x.sample_name == 'A', actual))
        list_of_values_for_a = [
            {'flowcell_id': 'foo', 'project_id': 'Project1', 'sample_name': 'A',
             'tag_seq': 'GTAGAGGA-CTCTCTAT', 'lane_num': 1, 'read_num': 1, 'cycles': 151,
             'pct_lane': None, 'pf_clusters': 0,
             'pct_q30': None, 'pct_tag_err': None,
             'library_name': 'Sample_A.library', 'mean_q': None},
            {'flowcell_id': 'foo', 'project_id': 'Project1', 'sample_name': 'A',
             'tag_seq': 'GTAGAGGA-CTCTCTAT', 'lane_num': 1, 'read_num': 2, 'cycles': 151,
             'pct_lane': None, 'pf_clusters': 0,
             'pct_q30': None, 'pct_tag_err': None,
             'library_name': 'Sample_A.library', 'mean_q': None}]

        expected_sample_a = list(map(lambda x: SampleResult(**x), list_of_values_for_a))
        self.assertListEqual(expected_sample_a, actual_sample_a)

if __name__ == '__main__':
    unittest.main(),
