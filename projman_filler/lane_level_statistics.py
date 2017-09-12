

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


def calculate_lane_statistics(flowcell_name, conversion_results, reads_and_cycles, error_rates, densities, q30s):
    for lane_dict in conversion_results:
        lane_nbr = lane_dict["LaneNumber"]
        total_clusters_raw = lane_dict["TotalClustersRaw"]
        total_clusters_pf = lane_dict["TotalClustersPF"]
        sample_demux_results = lane_dict["DemuxResults"]
        mean_q = "TODO" # TODO
        for sample_demux_result in sample_demux_results:
            sample_id = sample_demux_result["SampleId"]
            read_metrics = sample_demux_result["ReadMetrics"]
            for read_metric in read_metrics:
                read_nbr = read_metric["ReadNumber"]
                cycles = reads_and_cycles[read_nbr]

                error_rate = error_rates[lane_nbr][read_nbr]
                raw_density = densities[lane_nbr][read_nbr]["raw_density"]
                pf_density = densities[lane_nbr][read_nbr]["pass_filter_density"]
                percent_q30 = q30s[lane_nbr][read_nbr]
                yield LaneLevelStats(flowcell_name, lane_nbr, read_nbr, raw_density, pf_density, error_rate,
                                     total_clusters_raw, total_clusters_pf, cycles, percent_q30, mean_q)

