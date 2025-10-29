import pandas as pd
from math import isnan

from checkQC.parsers.illumina import _read_interop_summary
from projman_filler.models.db_models import FlowcellLaneResult, SampleResult


class QCDataParser:
    def __init__(self, qc_data, runfolder):
        """
        Initializes the QCDataParser with QC data.

        :param qc_data: QC data object from checkQC parser.
        :param runfolder: Path to the sequencing run folder.
        """
        self.run_summary, index_summary, run_info = _read_interop_summary(runfolder)
        self.flowcell_id = run_info.flowcell_id()
        self.samplesheet_df = pd.DataFrame(qc_data.samplesheet)
        self.qc_data = qc_data


    def _build_lane_results(self):
        """
        Constructs lane-level statistics from checkQC sequencing metrics results 
            from illumina parser.

        :return results (List[FlowcellLaneResult]): Lane-level statistics for non-index reads.
        """
        results = []
        for lane_no, lane_data in self.qc_data.sequencing_metrics.items():
            read_no = 0
            for number, read_data in lane_data["reads"].items():
                if read_data["is_index"]:
                    continue
                read_no += 1 # counts only non-index reads
                cycles = self.run_summary.at(number - 1).read().total_cycles()
                error_rate = None if read_data["mean_error_rate"] and \
                        isnan(read_data["mean_error_rate"]) else \
                        read_data["mean_error_rate"]

                results.append(FlowcellLaneResult(
                    flowcell_id=self.flowcell_id,
                    lane_num=lane_no,
                    read_num=read_no,
                    raw_density=lane_data["raw_density"],
                    pf_density=lane_data["pf_density"],
                    error_rate=error_rate,
                    pf_clusters=lane_data["total_reads_pf"],
                    raw_clusters=lane_data["total_reads"],
                    cycles=cycles,
                    pct_q30=read_data["percent_q30"] / 100,
                ))
        return results
    
    def _build_sample_results(self):
        """
        Constructs sample-level statistics by combining sequencing metrics, 
            sample sheet data results from the checkQC interop parser

        :return results (List[SampleResult]): Sample-level statistics including 
            quality scores, indexing accuracy, and library metadata.
        """
        results = []
        for lane_no, lane_data in self.qc_data.sequencing_metrics.items():
            for sample_data in lane_data["reads_per_sample"]:
                sample_id = sample_data["sample_id"]
                sample_row = self.samplesheet_df[
                    (self.samplesheet_df['lane'] == lane_no) &
                    (self.samplesheet_df['sample_id'] == sample_id)
                ]

                library_name = \
                    sample_row['custom_description'].to_string(index=False).split(
                        "LIBRARY_NAME:"
                    )[-1].strip()
                index1 = sample_row['index'].to_string(index=False)
                index2 = sample_row['index2'].to_string(index=False)
                sample_project = sample_row['sample_project'].to_string(index=False)

                read_no = 0
                for number, read_data in lane_data["reads"].items():
                    if read_data["is_index"]:
                        continue
                    read_no += 1 # counts only non-index reads
                    cycles = self.run_summary.at(number - 1).read().total_cycles()
                    no_of_reads = len(
                        [
                            no 
                            for no, read_data in lane_data["reads"].items()
                            if not read_data["is_index"]
                        ]
                    )
                    results.append(SampleResult(
                        flowcell_id=self.flowcell_id,
                        project_id=sample_project,
                        sample_name="_".join(sample_id.split("_")[1:]),
                        tag_seq=f"{index1}-{index2}" if index2 else index1,
                        lane_num=lane_no,
                        read_num=read_no,
                        cycles=cycles,
                        pct_lane=sample_data["percent_of_lane"],
                        pf_clusters=float(sample_data["cluster_count"]/no_of_reads),
                        pct_q30=sample_data["percent_q30"],
                        pct_tag_err=100 - sample_data["percent_perfect_index_reads"],
                        library_name=library_name,
                        mean_q=sample_data["mean_q30"],
                    ))
        return results
