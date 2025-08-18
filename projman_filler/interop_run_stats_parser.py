from interop import py_interop_run, py_interop_run_metrics, py_interop_summary
import interop as iop
import pandas as pd
from math import isnan
from projman_filler.run_stats_parser_interface import RunStatsParserInterface
from projman_filler.lane import Lane


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
    
    def get_reads_and_cycles(self) -> dict:
        """Returns number of cycles per read, e.g. {1: 151, 2: 151}."""
        reads_and_cycles = {}
        for read in self.get_reads():
            reads_and_cycles[read] = self.get_cycles(read)
        return reads_and_cycles

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

    def get_mismatch_counts(self, target_lane, target_tile) -> dict:
        """
        Get mismatch counts for each lane and tile.

        :param target_lane: Lane number to filter by.
        :param target_tile: Tile number to filter by.
        
        Returns a dictionary with keys as mismatch counts.
            Example: {"0": 100, "1": 50} means 100 perfect matches and 50 with one mismatch.
            If no mismatch counts are found, returns an empty dictionary
        """
        error_metrics = self._run_metrics.error_metric_set()
        mismatch_counts = {}

        # Iterate through error metrics to find the target lane and tile
        metrics_for_lane = error_metrics.metrics_for_lane(target_lane)

        for i in range(metrics_for_lane.size()):
            metric = error_metrics.at(i)
            tile = metric.tile()
            # Get mismatch cluster counts: index 0 = perfect, 1 = 1 mismatch, etc.
            counts = [
                metric.mismatch_cluster_count(n)
                for n in range(metric.mismatch_count())
                # for n in range(2)
                ]
            
            if counts:
                if "0" not in mismatch_counts:
                    mismatch_counts = {
                        "0": counts[0],
                        "1": counts[1],
                    }
                else:
                    mismatch_counts["0"] += counts[0]
                    mismatch_counts["1"] += counts[1]

        return mismatch_counts if mismatch_counts else None
    
    def _get_conversion_results(self) -> list:
        """
        Get conversion results from the run summary and metrics.
        Returns a list of Lane objects with detailed metrics.
        """
        # Prepare to load index, and read metrics
        run_metrics = self._run_metrics
        valid_to_load = py_interop_run.uchar_vector(py_interop_run.MetricCount, 0)
        py_interop_run_metrics.list_index_metrics_to_load(valid_to_load)
        # Reload metrics with updated valid_to_load vector
        run_metrics.read(self._runfolder, valid_to_load)

        # Summarize index metrics
        index_summary = py_interop_summary.index_flowcell_summary()
        py_interop_summary.summarize_index_metrics(run_metrics, index_summary)

        num_lanes = self._run_summary.at(0).size()
        lanes = []
        index_metric_set = run_metrics.index_metric_set()

        for lane_index in range(num_lanes):
            sample_demux_results = []
            num_reads = self._run_summary.at(0).at(lane_index).size()
            # Retrieve total reads and PF reads for the lane. index 0 is used 
            # because total reads and PF reads are the same for all reads in the lane.
            total_reads = self._run_summary.at(0).at(lane_index).reads()
            total_reads_pf = self._run_summary.at(0).at(lane_index).reads_pf()
            total_yield= total_reads_pf * self.get_cycles(1) * 2
            sample_demux = index_summary.at(lane_index)

            for sample_no in range(sample_demux_size:=sample_demux.size()):
                sample = sample_demux.at(sample_no)
                sample_index1 = f"{sample.index1()}"
                sample_index = sample_index1+f"{sample.index2()}" if sample.index2() else sample_index1
                sample_num_reads = sample.cluster_count() / len(self._non_index_reads)
                samples_yield = sample_num_reads * self.get_cycles(1) * 2
    
                mismatch_counts = self.get_mismatch_counts(
                    lane_index+1, index_metric_set.at(sample_no).tile()
                )

                # Construct index metrics
                index_metrics = [{
                    "IndexSequence": sample_index,
                    "MismatchCounts": mismatch_counts or {},
                }]
                read_metrics = []
                sample_demux_results.append({
                        "SampleName": "_".join(sample.sample_id().split("_")[1:]),
                        "SampleId": sample.sample_id(),
                        "IndexMetrics": index_metrics,
                        "ReadMetrics": read_metrics,
                        "Yield": samples_yield,
                        "NumberReads": sample_num_reads,
                    })

                for read_nbr in range(num_reads):
                    reads_per_sample = self._run_summary.at(0).at(lane_index).at(read_nbr)

                    # Get sample-level demux summary
                    try:
                        sample = sample_demux.at(read_nbr)
                        # Calculate yield in bases (from Gb)
                        read_yield = samples_yield / 2  #reads_per_sample.cluster_count().mean() * self.get_cycles(1) * 2
                        sample_demux_results[sample_no]["ReadMetrics"].append({
                            "ReadNumber": read_nbr + 1,
                            "Yield":read_yield,
                            "YieldQ30":  (
                                    reads_per_sample.yield_g() * 
                                    (reads_per_sample.percent_gt_q30() / 100.0)
                                )* 1e9,
                            "PercentQ30": reads_per_sample.percent_gt_q30(),
                            "PercentPF": reads_per_sample.percent_pf()
                        })

                    except iop.py_interop_metrics.index_out_of_bounds_exception:
                        # If the sample read is not found, continue to the next read
                        continue
            # Construct Lane object
            lanes.append(
                Lane(
                    lane_index+1, total_reads, total_reads_pf, yld=total_yield, 
                    sample_demux_results=sample_demux_results
                )
            )      

        return lanes
