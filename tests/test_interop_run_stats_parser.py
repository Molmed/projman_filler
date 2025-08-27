import unittest
import pytest
import os
from unittest.mock import MagicMock

from projman_filler.interop_run_stats_parser import InteropRunStatsParser
from projman_filler.lane_level_statistics import calculate_lane_statistics
from projman_filler.sample_level_statistics import calculate_sample_statistics
from projman_filler.samplesheet import Samplesheet
from projman_filler.bcl2fastq_run_stats_parser import Bcl2fastqRunStatsParser
import pandas as pd
import interop
from projman_filler.lane import Lane
from tests.test_utils import qc_data


class Mock_qc_data():
    def __init__(self, samplesheet, sequencing_metrics):
        self.samplesheet = samplesheet
        self.sequencing_metrics = sequencing_metrics

class TestRunStatsParsers(unittest.TestCase):
    def test_interop_standardize_read_numbers(self):
        runfolder = "tests/resources/200624_A00834_0183_BHMTFYTINY"

        non_index_reads = [0]
        iop = InteropRunStatsParser(runfolder, non_index_reads)
        remapped = iop._standardize_read_numbers()
        expected = {1: 0}
        self.assertEqual(remapped, expected)

        reads = iop.get_reads()
        expected = [1]
        self.assertEqual(reads, expected)

    def test_lanes_total_clusters(self):
        non_index_reads = [0]
        runfolder = "tests/resources/200624_A00834_0183_BHMTFYTINY"
        iop = InteropRunStatsParser(runfolder, non_index_reads)
        for lane in iop._conversion_results:
            assert lane._total_clusters_pf is not None and lane._total_clusters_pf != 0
            assert lane._total_clusters_raw is not None and lane._total_clusters_raw != 0

    def test_clusters_same_across_lanes(self):
        """
        Verify that 'Reads' and 'Reads Pf' are consistently the same across all reads within the same lane
        """
        non_index_reads = [0]
        runfolder = "tests/resources/200624_A00834_0183_BHMTFYTINY"
        iop = InteropRunStatsParser(runfolder, non_index_reads)
        
        interop_lane_summary = interop.summary(iop._run_metrics, 'Lane')
        df = pd.DataFrame(interop_lane_summary)

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
                
    def test_get_checkqc_interop_stats(self):
        runfolder = "tests/resources/200624_A00834_0184_BHMTFYTINY"
        iop = InteropRunStatsParser(runfolder)

        test_qc_data = Mock_qc_data(qc_data["samplesheet"], qc_data["sequencing_metrics"])
        flowcell_lane_results, sample_results = iop.get_checkqc_interop_stats(test_qc_data)

        assert len(flowcell_lane_results) == 2
        assert len(sample_results) == 4

    def test_bclconvert_results_similar_to_older_results(self):
        runfolder_v1 = "tests/resources/200624_A00834_0183_BHMTFYTINY"
        runfolder_v2 = "tests/resources/200624_A00834_0184_BHMTFYTINY"

        iop_v1 = InteropRunStatsParser(runfolder_v1)
        iop_v2 = InteropRunStatsParser(runfolder_v2)
        
        samplesheet_file = os.path.join(runfolder_v1, "SampleSheet.csv")
        samplesheet = Samplesheet(samplesheet_file)
        a =  f"{runfolder_v1}/Unaligned/Stats"
        bcl2fastq_stats = Bcl2fastqRunStatsParser(f"{runfolder_v1}/Unaligned/Stats")
        conversion_results = bcl2fastq_stats.get_conversion_results()
        flowcell_name = "HMTFYDRXX"
        reads_and_cycles = {1: 36}
        
        v1_lane_stats = calculate_lane_statistics(
            iop_v1, flowcell_name, conversion_results
        )
        v1_sample_stats = calculate_sample_statistics(
            flowcell_name, conversion_results, reads_and_cycles, samplesheet
        )
        v1_lane_stats = list(v1_lane_stats)
        v1_sample_stats = list(v1_sample_stats)

        test_qc_data = Mock_qc_data(qc_data["samplesheet"], qc_data["sequencing_metrics"])
        v2_lane_stats, v2_sample_stats = iop_v2.get_checkqc_interop_stats(test_qc_data)

        assert len(v1_lane_stats) == len(v2_lane_stats)
        for v1_results, v2_results in zip(v1_lane_stats, v2_lane_stats):
            assert v1_results.pct_q30 == pytest.approx(v2_results.pct_q30)
            assert v1_results.error_rate == pytest.approx(v2_results.error_rate)
            # rel(=(expected-obtained) / expected) gives percent relative difference to allow
            assert v1_results.pf_clusters == pytest.approx(v2_results.pf_clusters, rel=0.017)
            assert v1_results.pf_density == pytest.approx(v2_results.pf_density)
            assert v1_results.raw_clusters == pytest.approx(v2_results.raw_clusters)
            assert v1_results.raw_density == pytest.approx(v2_results.raw_density)

        for v1_results, v2_results in zip(v1_sample_stats, [v2_sample_stats[0], v2_sample_stats[2]]):
            assert round(v1_results.pct_q30, 2) == pytest.approx(v2_results.pct_q30, rel=0.005)
            assert round(v1_results.mean_q, 2) == pytest.approx(v2_results.mean_q)
            assert round(v1_results.pct_lane, 2) == pytest.approx(v2_results.pct_lane)
            assert round(v1_results.pct_tag_err, 2) == pytest.approx(round(v2_results.pct_tag_err, 2))
            assert v1_results.pf_clusters == pytest.approx(v2_results.pf_clusters)
            assert v1_results.project_id == pytest.approx(v2_results.project_id)
            assert v1_results.tag_seq == pytest.approx(v2_results.tag_seq)
            assert v1_results.lane_num == pytest.approx(v2_results.lane_num)
            assert v1_results.sample_name == pytest.approx(v2_results.sample_name)
            assert v1_results.library_name == pytest.approx(v2_results.library_name)

if __name__ == '__main__':
    unittest.main()

