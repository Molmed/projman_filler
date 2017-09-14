
from collections import defaultdict


class LaneLevelStats(object):
    def __init__(self, flowcell_id, lane_nbr, read_nbr, raw_density, pf_density, error_rate, total_clusters_raw,
                 total_clusters_pf, cycles, percent_q30, mean_q):
        self.flowcell_id = flowcell_id
        self.lane_nbr = lane_nbr
        self.read_nbr = read_nbr
        self.raw_density = raw_density
        self.pf_density = pf_density
        self.error_rate = error_rate
        self.total_clusters_raw = total_clusters_raw
        self.total_clusters_pf = total_clusters_pf
        self.cycles = cycles
        self.percent_q30 = percent_q30
        self.mean_q = mean_q

    def __str__(self):
        return str(self.__dict__)

    def __eq__(self, other):
        if isinstance(other, LaneLevelStats):
            if self.__dict__ == other.__dict__:
                return True
        else:
            return False


def _get_mean_q_scores(lane_dict):

    results = defaultdict(list)

    demultiplexing_results = lane_dict["DemuxResults"]
    undetermined_results = lane_dict["Undetermined"]

    demultiplex_info_dicts = demultiplexing_results + [undetermined_results]
    for sample_dict in demultiplex_info_dicts:
        read_metrics = sample_dict["ReadMetrics"]
        for read_metric in read_metrics:
            read_nbr = read_metric["ReadNumber"]
            results[read_nbr].append(float(read_metric["QualityScoreSum"]) / read_metric["Yield"])

    # Calculate the mean within each read
    for k, v in results.items():
        results[k] = sum(v) / float(len(v))
    return results


def calculate_lane_statistics(flowcell_name, conversion_results, reads_and_cycles, error_rates, densities, q30s):
    for lane_dict in conversion_results:
        lane_nbr = lane_dict["LaneNumber"]
        total_clusters_raw = lane_dict["TotalClustersRaw"]
        total_clusters_pf = lane_dict["TotalClustersPF"]
        mean_q = _get_mean_q_scores(lane_dict)

        for read_nbr in reads_and_cycles.keys():
            cycles = reads_and_cycles[read_nbr]

            error_rate = error_rates[lane_nbr][read_nbr]
            raw_density = densities[lane_nbr][read_nbr]["raw_density"]
            pf_density = densities[lane_nbr][read_nbr]["pass_filter_density"]
            percent_q30 = q30s[lane_nbr][read_nbr]
            mean_q_for_read = mean_q[read_nbr]
            yield LaneLevelStats(flowcell_name, lane_nbr, read_nbr, raw_density, pf_density, error_rate,
                                 total_clusters_raw, total_clusters_pf, cycles, percent_q30, mean_q_for_read)

