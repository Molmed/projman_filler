import unittest
from unittest.mock import MagicMock
from projman_filler.interop_run_stats_parser import InteropRunStatsParser
from interop import py_interop_run, py_interop_run_metrics, py_interop_summary
import pandas as pd
import interop
from projman_filler.lane import Lane


class TestRunStatsParsers(unittest.TestCase):
    def test_interop_standardize_read_numbers(self):
        runfolder = "tests/resources/200624_A00834_0183_BHMTFYTINY"

        non_index_reads = [0, 2, 3]
        iop = InteropRunStatsParser(runfolder, non_index_reads)
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

    def test_clusters_same_across_lanes(self):
        """
        Verify that 'Reads' and 'Reads Pf' are consistently the same across all reads within the same lane
        """
        non_index_reads = [0, 2, 3]
        runfolder = "tests/resources/200624_A00834_0183_BHMTFYTINY"
        iop = InteropRunStatsParser(runfolder, non_index_reads)
        

        data = {
                'ReadNumber': [1, 1, 2, 2, 3, 3],
                'IsIndex': [78, 78, 89, 89, 89, 89],
                'Lane': [1, 2, 1, 2, 1, 2],
                'Reads': [638337024.0] * 6,
                'Reads Pf': [532464320.0, 530917568.0] * 3,
            }

        df = pd.DataFrame(data)

        for lane_index, (lane, rows) in enumerate(df.groupby('Lane')):
            rows = rows.reset_index()
            # Assert all 'Reads' and 'Reads Pf' values are consistent within the lane
            assert rows['Reads'].nunique() == 1, (
                    f"Inconsistent 'Reads' in lane {lane}"
                    )
            assert rows['Reads Pf'].nunique() == 1, (
                    f"Inconsistent 'Reads Pf' in lane {lane}"
                    )

            # These are the same for the lane across all reads
            total_clusters_pf = rows.at[0, 'Reads Pf']
            total_clusters_raw = rows.at[0, 'Reads']
            
            
            # Compare with InteropRunStatsParser results
            lane_results = iop._conversion_results[lane_index]
            assert lane_results._total_clusters_pf == total_clusters_pf, (
                    f"Mismatch in total_clusters_pf for lane {lane}"
                    )
            assert lane_results._total_clusters_raw == total_clusters_raw, (
                    f"Mismatch in total_clusters_raw for lane {lane}"
                    )
                

if __name__ == '__main__':
    unittest.main()

