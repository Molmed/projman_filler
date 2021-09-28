from os import path
from collections import defaultdict
import json
from projman_filler.run_stats_parser_interface import RunStatsParserInterface
from projman_filler.lane import Lane


class Bcl2fastqRunStatsParser(RunStatsParserInterface):
    def __init__(self, stats_dir_path: str):
        json_path = path.join(stats_dir_path, 'Stats.json')
        with open(json_path, 'r') as f:
            self._stats = json.load(f)

    def get_flowcell_name(self) -> str:
        return self._stats["Flowcell"]

    def get_reads_and_cycles(self) -> dict:
        """Returns number of cycles per read, e.g. {1: 151, 2: 151}."""
        reads_and_cycles = {}
        for read_info in self._stats["ReadInfosForLanes"][0]["ReadInfos"]:
            if not read_info["IsIndexedRead"]:
                reads_and_cycles[read_info["Number"]] = read_info["NumCycles"]
        return reads_and_cycles

    def get_non_index_reads(self) -> list:
        non_index_reads = []
        stats_json_read_info = self._stats["ReadInfosForLanes"][0]["ReadInfos"]
        for read_nbr, read_info in enumerate(stats_json_read_info):
            if not read_info["IsIndexedRead"]:
                non_index_reads.append(read_nbr)
        return non_index_reads

    def get_conversion_results(self) -> list:
        lanes = []
        for lane_dict in self._stats["ConversionResults"]:
            lane_nbr = lane_dict["LaneNumber"]
            total_clusters_raw = lane_dict["TotalClustersRaw"]
            total_clusters_pf = lane_dict["TotalClustersPF"]
            mean_q_scores = self._get_mean_q_scores(lane_dict)
            lane_yield = lane_dict["Yield"]
            sample_demux_results = lane_dict["DemuxResults"]
            lanes.append(Lane(lane_nbr, total_clusters_raw, total_clusters_pf, mean_q_scores, lane_yield, sample_demux_results))
        return lanes

    # The mean_q_scores metric has been ported here as a remnant from the old
    # Sisyphus system. It is unclear if it is used, so this logic has not been
    # implemented in the interop parser.
    def _get_mean_q_scores(self, lane_dict) -> list:
        results = defaultdict(list)

        demultiplexing_results = lane_dict["DemuxResults"]
        if lane_dict.get("Undetermined"):
            undetermined_results = lane_dict["Undetermined"]
            demultiplex_info_dicts = demultiplexing_results + [undetermined_results]
        else:
            demultiplex_info_dicts = demultiplexing_results

        for sample_dict in demultiplex_info_dicts:
            read_metrics = sample_dict["ReadMetrics"]
            for read_metric in read_metrics:
                read_nbr = read_metric["ReadNumber"]
                # If yield is 0 for a sample we do not add it to the mean score
                if read_metric["Yield"] != 0:
                    results[read_nbr].append(float(read_metric["QualityScoreSum"]) / read_metric["Yield"])

        # Calculate the mean within each read
        for k, v in results.items():
            results[k] = sum(v) / float(len(v))
        return results
