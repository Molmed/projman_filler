import unittest

from tests.test_utils import conversion_results

from projman_filler.sample_level_statistics import calculate_sample_statistics,SampleLevelStats


class TestSampleLevelStatistics(unittest.TestCase):
    flowcell_name = "foo"
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
        actual = list(calculate_sample_statistics(flowcell_name=self.flowcell_name,
                                                  conversion_results=conversion_results,
                                                  reads_and_cycles=self.reads_and_cycles,
                                                  samplesheet=self.samplesheet_mock))

        # One row for per sample, index (in the test data there are 3 per sample), and 2 reads.
        self.assertEqual(len(actual), 4*3*2)

        actual_sample_a = list(filter(lambda x: x.sample_name == 'A', actual))
        list_of_values_for_a = [
            {'flowcell_name': 'foo', 'sample_project': 'Project1', 'sample_name': 'A',
                           'sample_tag': 'GTAGAGGA+CTCTCTAT', 'lane_nbr': 1, 'read_nbr': 1, 'cycles': 151,
                           'percent_of_lane': 49.91040361971908, 'pass_filter_clusters': 81217423.0,
                           'percent_q_30': 98.02429332249935, 'percent_tag_error': 0.671112157794024,
                           'sample_library_name': 'A.library', 'mean_q': 38.84148990743496},
                          {'flowcell_name': 'foo', 'sample_project': 'Project1', 'sample_name': 'A',
                           'sample_tag': 'GTAGAGGA+CTCTCTAT', 'lane_nbr': 1, 'read_nbr': 2, 'cycles': 151,
                           'percent_of_lane': 49.91040361971908, 'pass_filter_clusters': 81217423.0,
                           'percent_q_30': 96.45192508767363, 'percent_tag_error': 0.671112157794024,
                           'sample_library_name': 'A.library', 'mean_q': 38.373262536376345},
                          {'flowcell_name': 'foo', 'sample_project': 'Project1', 'sample_name': 'A',
                           'sample_tag': 'TAGGCATG+CTCTCTAT', 'lane_nbr': 1, 'read_nbr': 1, 'cycles': 151,
                           'percent_of_lane': 49.91040361971908, 'pass_filter_clusters': 81217423.0,
                           'percent_q_30': 98.02429332249935, 'percent_tag_error': 0.7880181078880083,
                           'sample_library_name': 'A.library', 'mean_q': 38.84148990743496},
                          {'flowcell_name': 'foo', 'sample_project': 'Project1', 'sample_name': 'A',
                           'sample_tag': 'TAGGCATG+CTCTCTAT', 'lane_nbr': 1, 'read_nbr': 2, 'cycles': 151,
                           'percent_of_lane': 49.91040361971908, 'pass_filter_clusters': 81217423.0,
                           'percent_q_30': 96.45192508767363, 'percent_tag_error': 0.7880181078880083,
                           'sample_library_name': 'A.library', 'mean_q': 38.373262536376345},
                          {'flowcell_name': 'foo', 'sample_project': 'Project1', 'sample_name': 'A',
                           'sample_tag': 'TCCTGAGC+CTCTCTAT', 'lane_nbr': 1, 'read_nbr': 1, 'cycles': 151,
                           'percent_of_lane': 49.91040361971908, 'pass_filter_clusters': 81217423.0,
                           'percent_q_30': 98.02429332249935, 'percent_tag_error': 0.7687463809335591,
                           'sample_library_name': 'A.library', 'mean_q': 38.84148990743496},
                          {'flowcell_name': 'foo', 'sample_project': 'Project1', 'sample_name': 'A',
                           'sample_tag': 'TCCTGAGC+CTCTCTAT', 'lane_nbr': 1, 'read_nbr': 2, 'cycles': 151,
                           'percent_of_lane': 49.91040361971908, 'pass_filter_clusters': 81217423.0,
                           'percent_q_30': 96.45192508767363, 'percent_tag_error': 0.7687463809335591,
                           'sample_library_name': 'A.library', 'mean_q': 38.373262536376345}]

        expected_sample_a = list(map(lambda x: SampleLevelStats(**x), list_of_values_for_a))
        self.assertListEqual(expected_sample_a, actual_sample_a)

if __name__ == '__main__':
    unittest.main(),
