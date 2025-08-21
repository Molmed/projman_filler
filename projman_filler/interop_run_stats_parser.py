from interop import py_interop_run, py_interop_run_metrics, py_interop_summary
import interop as iop
import pandas as pd
from math import isnan
from projman_filler.run_stats_parser_interface import RunStatsParserInterface
from projman_filler.lane import Lane
from projman_filler.models.db_models import FlowcellLaneResult, SampleResult
from checkQC.parsers.illumina import _read_interop_summary


class InteropRunStatsParser(RunStatsParserInterface):
    # We allow overriding the Interop list of non-index-reads with
    # a custom list obtained from bcl2fastq stats.
    def __init__(self, runfolder: str, non_index_reads: list = None):
        run_info = py_interop_run.info()
        run_info.read(runfolder)

        run_metrics = py_interop_run_metrics.run_metrics()
        run_metrics.run_info()

        valid_to_load = py_interop_run.uchar_vector(py_interop_run.MetricCount, 0)
        py_interop_run_metrics.list_summary_metrics_to_load(valid_to_load)
        run_metrics.read(runfolder, valid_to_load)

        summary = py_interop_summary.run_summary()
        py_interop_summary.summarize_run_metrics(run_metrics, summary)

        self._runfolder = runfolder
        self._run_info = run_info
        self._run_metrics = run_metrics
        self._run_summary = summary
        self._num_reads = self._run_summary.size()

        # non-index reads
        if non_index_reads:
            self._non_index_reads = non_index_reads
        else:
            self._non_index_reads = self._get_non_index_reads()

        self._standard_read_map = self._standardize_read_numbers()

        # conversion results
        self._conversion_results = self._get_conversion_results()

    def get_flowcell_name(self) -> str:
        return self._run_info.flowcell_id()

    def get_non_index_reads(self) -> list:
        return self._non_index_reads

    # Returns list of reads in standardized format
    def get_reads(self) -> list:
        return list(self._standard_read_map.keys())

    def get_cycles(self, read_number) -> int:
        read = self._get_interop_read_by_index(self._interop_read_number(read_number))
        return self._clean(read.total_cycles())

    def get_error_rate(self, read_number, lane_number) -> float:
        lane = self._get_interop_lane(read_number, lane_number)
        return self._clean(lane.error_rate().mean())

    def get_raw_density(self, read_number, lane_number) -> float:
        lane = self._get_interop_lane(read_number, lane_number)
        return self._clean(lane.density().mean())

    def get_pf_density(self, read_number, lane_number) -> float:
        lane = self._get_interop_lane(read_number, lane_number)
        return self._clean(lane.density_pf().mean())

    def get_pct_gt_q30(self, read_number, lane_number) -> float:
        lane = self._get_interop_lane(read_number, lane_number)
        # This will be given as a fraction to be compatible with legacy data /JD 2017-10-04
        return self._clean(float(lane.percent_gt_q30()) / 100)

    def get_conversion_results(self) -> list:
        return self._conversion_results

    # Convert nan's to Nones
    def _clean(self, val):
        if isnan(val):
            return None
        return val

    # "Renumber" the reads to have their intuitive read numbers, i.e. Read 1 -> 1, Read 2 -> 2
    def _standardize_read_numbers(self) -> list:
        mapping = {}
        for new_nbr, old_read in enumerate(self._non_index_reads):
            mapping[new_nbr + 1] = old_read
        return mapping

    # Takes the "intuitive" read number and returns the interop read number
    def _interop_read_number(self, read_number) -> int:
        return self._standard_read_map[read_number]

    def _get_interop_read(self, read_number) -> object:
        return self._run_summary.at(self._interop_read_number(read_number))

    # read_index is the index of a read in an array,
    # not to be confused with an index read, the Illumina sequencing term.
    def _get_interop_read_by_index(self, read_index) -> object:
        return self._run_summary.at(read_index).read()

    def _get_interop_lane(self, read_number, lane_number) -> object:
        return self._get_interop_read(read_number).at(lane_number - 1)

    def _get_non_index_reads(self) -> list:
        # Get non-index reads from interop
        non_index_reads = []
        for r in range(0, self._num_reads):
            is_index_read = self._get_interop_read_by_index(r).is_index()
            if not is_index_read:
                non_index_reads.append(r)
        return non_index_reads

    def _get_conversion_results(self) -> list:
        ar = iop.summary(self._run_metrics, 'Lane')
        df = pd.DataFrame(ar)
        # Get statistics per-lane
        n_lanes = self._run_summary.lane_count()
        lanes = []

        for l in range(1, n_lanes+1):
            rows = df.loc[df['Lane'] == l][['ReadNumber','Reads', 'Reads Pf', 'IsIndex']]
            rows = rows.reset_index()
            
            # These are the same for the lane across all reads
            total_clusters_pf = rows.iloc[0].get('Reads Pf', 0)
            total_clusters_raw = rows.iloc[0].get('Reads', 0)
            
            lanes.append(Lane(l, total_clusters_raw, total_clusters_pf))
        return lanes
    
    def get_checkqc_interop_stats(self, qc_data):
        """
        Gets run stats from checkqc illumina parser and the rest from iterop

        Params:
        :qc_data (acheckQC QCData object): results from checkQC illumina parser 
            for specified demultiplexer

        Returns:
        : tuple of flowcell_lane_result, sample_results to be added to respective DB
        """
        flowcell_lane_results = []
        sample_results = []
        
        run_summary, index_summary = _read_interop_summary(self._runfolder)
        flowcell_id = self.get_flowcell_name()  #self._run_info.flowcell_id()
        samplesheet_pd = pd.DataFrame(qc_data.samplesheet)

        for lane_no, lane_dict in qc_data.sequencing_metrics.items():
            for read_no, read_dict in lane_dict["reads"].items():
                if read_dict["is_index"]:
                    continue
                cycles = run_summary.at(read_no-1).read().total_cycles()

                flowcell_lane_results.append(
                    FlowcellLaneResult(
                        flowcell_id=flowcell_id,
                        lane_num=lane_no,
                        read_num=read_no,
                        raw_density=lane_dict["raw_density"],
                        pf_density=lane_dict["pf_density"],
                        error_rate=self._clean(read_dict["mean_error_rate"]),
                        pf_clusters=lane_dict["total_cluster_pf"],
                        cycles=cycles,
                        pct_q30=read_dict["percent_q30"],
                    )
                )

            for sample_dict in lane_dict["reads_per_sample"]:
                for read_no, read_dict in lane_dict["reads"].items():
                    if read_dict["is_index"]:
                        continue
                    sample_id = sample_dict["sample_id"]
                    sample_row = \
                        samplesheet_pd[
                            (samplesheet_pd['lane'] == lane_no) &
                            (samplesheet_pd['sample_id'] == sample_id) 
                        ]
                    library_name = sample_row['description'].to_string(index=False).split("LIBRARY_NAME:")[-1].strip()
                    index1 = sample_row['index'].to_string(index=False)
                    index2 = sample_row['index'].to_string(index=False)
                    sample_project = sample_row['sample_project'].to_string(index=False)

                    sample_results.append(
                        SampleResult(
                            flowcell_id=flowcell_id,
                            project_id=sample_project,
                            sample_name="_".join(sample_id.split("_")[1:]),
                            tag_seq=f"{index1}-{index2}" \
                                if index2 else index1,
                            lane_num=lane_no,
                            read_num=read_no,
                            cycles=cycles,
                            pct_lane=sample_dict["percent_of_lane"],
                            pf_clusters=float(sample_dict["cluster_count"]),
                            pct_q30=sample_dict["percent_q30"] * 100,
                            pct_tag_err=100-sample_dict["percent_perfect_index_reads"],
                            library_name=library_name,
                            mean_q=sample_dict["mean_q30"],
                        )
                    )
    
        return flowcell_lane_results, sample_results
