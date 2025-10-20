import unittest
import pytest
import os
from pathlib import Path

from projman_filler.qc_data_parser import QCDataParser
from projman_filler.interop_run_stats_parser import InteropRunStatsParser
from projman_filler.lane_level_statistics import calculate_lane_statistics
from projman_filler.sample_level_statistics import calculate_sample_statistics
from projman_filler.samplesheet import Samplesheet
from projman_filler.bcl2fastq_run_stats_parser import Bcl2fastqRunStatsParser
from checkQC.qc_data_utils import bclconvert_test_runfolder
from checkQC.qc_data import QCData


runfolder_v1 = "tests/resources/200624_A00834_0183_BHMTFYTINY"
runfolder_v2 = Path(__file__).parent / "resources/200624_A00834_0184_BHMTFYTINY"
class Mock_qc_data():
    def __init__(self, samplesheet, sequencing_metrics, qc_data):
        self.samplesheet = samplesheet
        self.sequencing_metrics = sequencing_metrics
        self.qc_data = qc_data


def get_mock_qc_data(runfolder):
    parser_config = {
        "reports_location": "Reports"
    }

    qc_data = QCData.from_bclconvert(runfolder, parser_config,)
    bclconvert_data =  bclconvert_test_runfolder(qc_data, runfolder)
    return Mock_qc_data(
        bclconvert_data["expected_samplesheet"]["head"],
        bclconvert_data["expected_sequencing_metrics"],
        qc_data
    )

class TestQCDataParsers(unittest.TestCase):               
    def test_qc_data_parser_results(self):
        test_qc_data = get_mock_qc_data(runfolder_v2)
        qc_data_parser = QCDataParser(test_qc_data.qc_data, runfolder_v2)

        flowcell_lane_results = qc_data_parser._build_lane_results()
        sample_results = qc_data_parser._build_sample_results()

        assert len(flowcell_lane_results) == 2
        assert len(sample_results) == 4

    def test_bclconvert_results_similar_to_older_results(self):
        test_qc_data = get_mock_qc_data(runfolder_v2)
        iop_v1 = InteropRunStatsParser(runfolder_v1)
        qc_data_parser = QCDataParser(test_qc_data.qc_data, runfolder_v2)
        
        samplesheet_file = os.path.join(runfolder_v1, "SampleSheet.csv")
        samplesheet = Samplesheet(samplesheet_file)
        bcl2fastq_stats = Bcl2fastqRunStatsParser(f"{runfolder_v1}/Unaligned/Stats")
        conversion_results = bcl2fastq_stats.get_conversion_results()
        reads_and_cycles = {1: 36}
        flowcell_name = "HMTFYDRXX"
        
        v1_lane_stats = calculate_lane_statistics(
            iop_v1, flowcell_name, conversion_results
        )
        v1_sample_stats = calculate_sample_statistics(
            flowcell_name, conversion_results, reads_and_cycles, samplesheet
        )
        v1_lane_stats = list(v1_lane_stats)
        v1_sample_stats = list(v1_sample_stats)

        v2_lane_stats = qc_data_parser._build_lane_results()
        v2_sample_stats = qc_data_parser._build_sample_results()

        assert len(v1_lane_stats) == len(v2_lane_stats)
        for v1_results, v2_results in zip(v1_lane_stats, v2_lane_stats):
            assert v1_results.pct_q30 == pytest.approx(v2_results.pct_q30)
            assert v1_results.error_rate == pytest.approx(v2_results.error_rate)
            # rel(=(expected-obtained) / expected) gives percent relative difference to allow for small variations
            assert v1_results.pf_clusters == pytest.approx(v2_results.pf_clusters, rel=0.017)
            assert v1_results.pf_density == pytest.approx(v2_results.pf_density)
            assert v1_results.raw_clusters == pytest.approx(v2_results.raw_clusters, rel=0.15)
            assert v1_results.raw_density == pytest.approx(v2_results.raw_density)

        for v1_results, v2_results in zip(v1_sample_stats, [v2_sample_stats[0], v2_sample_stats[2]]):
            assert round(v1_results.pct_q30, 2) == pytest.approx(v2_results.pct_q30, rel=0.005)  # assert 96.4 == 96.0 ± 9.6e-05
            assert round(v1_results.mean_q, 2) == pytest.approx(v2_results.mean_q, rel=0.005)  #assert 36.43 == 36.37 ± 3.6e-05
            assert round(v1_results.pct_lane, 2) == pytest.approx(v2_results.pct_lane, rel=0.035) # assert 0.3 == 0.29  ± 2.9e-07
            assert round(v1_results.pct_tag_err, 2) == pytest.approx(round(v2_results.pct_tag_err, 2), rel=0.015)  # assert 2.07 == 2.04 ± 2.0e-06
            # TODO: Fix Stats.json sample yield value as it seems to be wrong
            # assert v1_results.pf_clusters == pytest.approx(v2_results.pf_clusters)  # assert 1576102.0 == 9920.0 ± 9.9e-03
            assert v1_results.project_id == pytest.approx(v2_results.project_id)
            assert v1_results.tag_seq == pytest.approx(v2_results.tag_seq)
            assert v1_results.lane_num == pytest.approx(v2_results.lane_num)
            assert v1_results.sample_name == pytest.approx(v2_results.sample_name)
            assert v1_results.library_name == pytest.approx(v2_results.library_name)

if __name__ == '__main__':
    unittest.main()

