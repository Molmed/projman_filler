import unittest

from projman_filler.lane_level_statistics import calculate_lane_statistics
from projman_filler.lane_level_statistics import FlowcellLaneResult

from tests.test_utils import conversion_results, conversion_results_without_undetermined, conversion_results_empty_lane


class TestLaneLevelStatistics(unittest.TestCase):
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
    q30s = {1: {1: 95.0, 2: 92.0}, 2: {1: 96.1, 2: 94.0}}

    def test_calculate_lane_level_stats(self):
        actual = list(calculate_lane_statistics(flowcell_name=self.flowcell_name,
                                                conversion_results=conversion_results,
                                                reads_and_cycles=self.reads_and_cycles,
                                                error_rates=self.error_rates,
                                                densities=self.densities,
                                                q30s=self.q30s))
        # The test data has
        self.assertEqual(len(actual), 4)
        expected = [
            FlowcellLaneResult(flowcell_id='foo', lane_num=1, read_num=1, raw_density=100000, pf_density=100000,
                               error_rate=2.0, raw_clusters=168865204, pf_clusters=162726440, cycles=151,
                               pct_q30=0.95, mean_q=38.778728523122744),
            FlowcellLaneResult(flowcell_id='foo', lane_num=1, read_num=2, raw_density=100000, pf_density=100000,
                               error_rate=2.0, raw_clusters=168865204, pf_clusters=162726440, cycles=151,
                               pct_q30=0.92, mean_q=38.23724053948388),
            FlowcellLaneResult(flowcell_id='foo', lane_num=2, read_num=1, raw_density=100000, pf_density=100000,
                               error_rate=2.0, raw_clusters=170966905, pf_clusters=164470667, cycles=151,
                               pct_q30=0.961, mean_q=38.74961137133973),
            FlowcellLaneResult(flowcell_id='foo', lane_num=2, read_num=2, raw_density=100000, pf_density=100000,
                               error_rate=2.0, raw_clusters=170966905, pf_clusters=164470667, cycles=151,
                               pct_q30=0.94, mean_q=38.18423664140531)
        ]
        self.assertListEqual(actual, expected)

    def test_calculate_lane_level_stats_no_undetermined(self):
        actual = list(calculate_lane_statistics(flowcell_name=self.flowcell_name,
                                                conversion_results=conversion_results_without_undetermined,
                                                reads_and_cycles=self.reads_and_cycles,
                                                error_rates=self.error_rates,
                                                densities=self.densities,
                                                q30s=self.q30s))
        # The test data has
        self.assertEqual(len(actual), 4)
        expected = [
            FlowcellLaneResult(flowcell_id='foo', lane_num=1, read_num=1, raw_density=100000, pf_density=100000,
                               error_rate=2.0, raw_clusters=168865204, pf_clusters=162726440, cycles=151,
                               pct_q30=0.95, mean_q=38.838484242499376),
            FlowcellLaneResult(flowcell_id='foo', lane_num=1, read_num=2, raw_density=100000, pf_density=100000,
                               error_rate=2.0, raw_clusters=168865204, pf_clusters=162726440, cycles=151,
                               pct_q30=0.92, mean_q=38.347822285262204),
            FlowcellLaneResult(flowcell_id='foo', lane_num=2, read_num=1, raw_density=100000, pf_density=100000,
                               error_rate=2.0, raw_clusters=170966905, pf_clusters=164470667, cycles=151,
                               pct_q30=0.961, mean_q=38.81856298771521),
            FlowcellLaneResult(flowcell_id='foo', lane_num=2, read_num=2, raw_density=100000, pf_density=100000,
                               error_rate=2.0, raw_clusters=170966905, pf_clusters=164470667, cycles=151,
                               pct_q30=0.94, mean_q=38.309461018362995)
        ]
        self.assertListEqual(actual, expected)

    def test_calculate_lane_level_stats_empty_lane(self):
        actual = list(calculate_lane_statistics(flowcell_name=self.flowcell_name,
                                                conversion_results=conversion_results_empty_lane,
                                                reads_and_cycles=self.reads_and_cycles,
                                                error_rates=self.error_rates,
                                                densities=self.densities,
                                                q30s={1: {1: 0, 2: 0}, 2: {1: 96.1, 2: 94.0}}))
        # The test data has
        self.assertEqual(len(actual), 4)
        expected = [
            FlowcellLaneResult(flowcell_id='foo', lane_num=1, read_num=1, raw_density=100000, pf_density=100000,
                               error_rate=2.0, raw_clusters=168865204, pf_clusters=0, cycles=151,
                               pct_q30=0.0, mean_q=None),
            FlowcellLaneResult(flowcell_id='foo', lane_num=1, read_num=2, raw_density=100000, pf_density=100000,
                               error_rate=2.0, raw_clusters=168865204, pf_clusters=0, cycles=151,
                               pct_q30=0.0, mean_q=None),
            FlowcellLaneResult(flowcell_id='foo', lane_num=2, read_num=1, raw_density=100000, pf_density=100000,
                               error_rate=2.0, raw_clusters=170966905, pf_clusters=164470667, cycles=151,
                               pct_q30=0.961, mean_q=38.74961137133973),
            FlowcellLaneResult(flowcell_id='foo', lane_num=2, read_num=2, raw_density=100000, pf_density=100000,
                               error_rate=2.0, raw_clusters=170966905, pf_clusters=164470667, cycles=151,
                               pct_q30=0.94, mean_q=38.18423664140531)
        ]
        self.assertListEqual(actual, expected)

if __name__ == '__main__':
    unittest.main()
