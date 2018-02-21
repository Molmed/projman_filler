import unittest

from projman_filler.models.db_models import SampleResult
from projman_filler.sample_level_statistics import calculate_sample_statistics

from tests.test_utils import conversion_results, conversion_results_without_index_metrics, conversion_results_sample_with_no_reads


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
            self.project_dict = {"A": "Project1",
                                 "B": "Project2",
                                 "C": "Project1",
                                 "D": "Project2"}

        def project_for_sample(self, sample_name, lane):
            return self.project_dict[sample_name]

        def library_name_for_sample(self, sample_name, lane):
            return "{}.library".format(sample_name)

    samplesheet_mock = SampleSheetMock()

    def test_calculate_sample_level_statistics(self):
        actual = list(calculate_sample_statistics(flowcell_name=self.flowcell_id,
                                                  conversion_results=conversion_results,
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
             'library_name': 'A.library', 'mean_q': 38.84148990743496},
            {'flowcell_id': 'foo', 'project_id': 'Project1', 'sample_name': 'A',
             'tag_seq': 'GTAGAGGA-CTCTCTAT', 'lane_num': 1, 'read_num': 2, 'cycles': 151,
             'pct_lane': 49.91040361971908, 'pf_clusters': 81217423.0,
             'pct_q30': 96.45192508767363, 'pct_tag_err': 0.671112157794024,
             'library_name': 'A.library', 'mean_q': 38.373262536376345},
            {'flowcell_id': 'foo', 'project_id': 'Project1', 'sample_name': 'A',
             'tag_seq': 'TAGGCATG-CTCTCTAT', 'lane_num': 1, 'read_num': 1, 'cycles': 151,
             'pct_lane': 49.91040361971908, 'pf_clusters': 81217423.0,
             'pct_q30': 98.02429332249935, 'pct_tag_err': 0.7880181078880083,
             'library_name': 'A.library', 'mean_q': 38.84148990743496},
            {'flowcell_id': 'foo', 'project_id': 'Project1', 'sample_name': 'A',
             'tag_seq': 'TAGGCATG-CTCTCTAT', 'lane_num': 1, 'read_num': 2, 'cycles': 151,
             'pct_lane': 49.91040361971908, 'pf_clusters': 81217423.0,
             'pct_q30': 96.45192508767363, 'pct_tag_err': 0.7880181078880083,
             'library_name': 'A.library', 'mean_q': 38.373262536376345},
            {'flowcell_id': 'foo', 'project_id': 'Project1', 'sample_name': 'A',
             'tag_seq': 'TCCTGAGC-CTCTCTAT', 'lane_num': 1, 'read_num': 1, 'cycles': 151,
             'pct_lane': 49.91040361971908, 'pf_clusters': 81217423.0,
             'pct_q30': 98.02429332249935, 'pct_tag_err': 0.7687463809335591,
             'library_name': 'A.library', 'mean_q': 38.84148990743496},
            {'flowcell_id': 'foo', 'project_id': 'Project1', 'sample_name': 'A',
             'tag_seq': 'TCCTGAGC-CTCTCTAT', 'lane_num': 1, 'read_num': 2, 'cycles': 151,
             'pct_lane': 49.91040361971908, 'pf_clusters': 81217423.0,
             'pct_q30': 96.45192508767363, 'pct_tag_err': 0.7687463809335591,
             'library_name': 'A.library', 'mean_q': 38.373262536376345}]

        expected_sample_a = list(map(lambda x: SampleResult(**x), list_of_values_for_a))
        self.assertListEqual(expected_sample_a, actual_sample_a)

    def test_calculate_sample_level_statistics_without_index_metrics(self):
        actual = list(calculate_sample_statistics(flowcell_name=self.flowcell_id,
                                                  conversion_results=conversion_results_without_index_metrics,
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
             'library_name': 'A.library', 'mean_q': 38.84148990743496},
            {'flowcell_id': 'foo', 'project_id': 'Project1', 'sample_name': 'A',
             'tag_seq': 'unknown', 'lane_num': 1, 'read_num': 2, 'cycles': 151,
             'pct_lane': 49.91040361971908, 'pf_clusters': 81217423.0,
             'pct_q30': 96.45192508767363, 'pct_tag_err': None,
             'library_name': 'A.library', 'mean_q': 38.373262536376345}]

        expected_sample_a = list(map(lambda x: SampleResult(**x), list_of_values_for_a))
        self.assertListEqual(expected_sample_a, actual_sample_a)

    def test_calculate_sample_level_statistics_sample_with_no_reads(self):
        actual = list(calculate_sample_statistics(flowcell_name=self.flowcell_id,
                                                  conversion_results=conversion_results_sample_with_no_reads,
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
             'library_name': 'A.library', 'mean_q': None},
            {'flowcell_id': 'foo', 'project_id': 'Project1', 'sample_name': 'A',
             'tag_seq': 'GTAGAGGA-CTCTCTAT', 'lane_num': 1, 'read_num': 2, 'cycles': 151,
             'pct_lane': 0, 'pf_clusters': 0,
             'pct_q30': None, 'pct_tag_err': None,
             'library_name': 'A.library', 'mean_q': None}]

        expected_sample_a = list(map(lambda x: SampleResult(**x), list_of_values_for_a))
        self.assertListEqual(expected_sample_a, actual_sample_a)

if __name__ == '__main__':
    unittest.main(),
