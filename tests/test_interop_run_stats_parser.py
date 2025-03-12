import unittest
from unittest.mock import MagicMock
from projman_filler.interop_run_stats_parser import InteropRunStatsParser
from interop import py_interop_run, py_interop_run_metrics, py_interop_summary
import interop

class TestRunStatsParsers(unittest.TestCase):
    def _setup_mocks(self):
        py_interop_run.info = MagicMock()
        py_interop_run_metrics.run_metrics = MagicMock()
        py_interop_summary.summarize_run_metrics = MagicMock()
        interop.summary = MagicMock()

    def test_interop_standardize_read_numbers(self):
        self._setup_mocks()

        non_index_reads = [0, 2, 3]
        iop = InteropRunStatsParser("foo", non_index_reads)
        remapped = iop._standardize_read_numbers()
        expected = {1: 0, 2: 2, 3: 3}
        self.assertEqual(remapped, expected)

        reads = iop.get_reads()
        expected = [1, 2, 3]
        self.assertEqual(reads, expected)

    def test_lanes_total_clusters(self):
        non_index_reads = [0, 2, 3]
        runfolder = "tests/resources/200624_A00834_0183_BHMTFYTINY"
        iop = InteropRunStatsParser(runfolder, non_index_reads)
        for lane in iop._conversion_results:
            assert lane._total_clusters_pf is not None and lane._total_clusters_pf != 0
            assert lane._total_clusters_raw is not None and lane._total_clusters_raw != 0

if __name__ == '__main__':
    unittest.main()
