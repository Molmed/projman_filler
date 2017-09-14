import unittest

from projman_filler.lane_level_statistics import calculate_lane_statistics
from projman_filler.lane_level_statistics import LaneLevelStats

from tests.test_utils import conversion_results


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
    q30s = {1: {1: 95.0, 2: 92.0}, 2: {1: 96.1, 2: 94.2}}

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
            LaneLevelStats(flowcell_id='foo', lane_nbr=1, read_nbr=1, raw_density=100000, pf_density=100000,
                           error_rate=2.0, total_clusters_raw=168865204, total_clusters_pf=162726440, cycles=151,
                           percent_q30=95.0, mean_q=38.778728523122744),
            LaneLevelStats(flowcell_id='foo', lane_nbr=1, read_nbr=2, raw_density=100000, pf_density=100000,
                           error_rate=2.0, total_clusters_raw=168865204, total_clusters_pf=162726440, cycles=151,
                           percent_q30=92.0, mean_q=38.23724053948388),
            LaneLevelStats(flowcell_id='foo', lane_nbr=2, read_nbr=1, raw_density=100000, pf_density=100000,
                           error_rate=2.0, total_clusters_raw=170966905, total_clusters_pf=164470667, cycles=151,
                           percent_q30=96.1, mean_q=38.74961137133973),
            LaneLevelStats(flowcell_id='foo', lane_nbr=2, read_nbr=2, raw_density=100000, pf_density=100000,
                           error_rate=2.0, total_clusters_raw=170966905, total_clusters_pf=164470667, cycles=151,
                           percent_q30=94.2, mean_q=38.18423664140531)
        ]
        self.assertListEqual(actual, expected)


if __name__ == '__main__':
    unittest.main()
