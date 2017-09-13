

class SampleLevelStats(object):

    def __init__(self, flowcell_name, sample_project, sample_name, sample_tag, lane_nbr, read_nbr, cycles,
                 percent_of_lane, pass_filter_clusters, percent_q_30, percent_tag_error, sample_library_name,
                 mean_q):
        self.flowcell_name = flowcell_name
        self.sample_project = sample_project
        self.sample_name = sample_name
        self.sample_tag = sample_tag
        self.lane_nbr = lane_nbr
        self.read_nbr = read_nbr
        self.cycles = cycles
        self.percent_of_lane = percent_of_lane
        self.pass_filter_clusters = pass_filter_clusters
        self.percent_q_30 = percent_q_30
        self.percent_tag_error = percent_tag_error
        self.sample_library_name = sample_library_name
        self.mean_q = mean_q

    def __str__(self):
        return str(self.__dict__)


def fetch_index_sequence(sample_demux_result):
    indexes = []
    for index in sample_demux_result["IndexMetrics"]:
        indexes.append(index["IndexSequence"])
    return "-".join(indexes)


def calculate_sample_statistics(flowcell_name, conversion_results, reads_and_cycles, samplesheet):
    for lane_dict in conversion_results:

        lane_nbr = lane_dict["LaneNumber"]
        total_clusters_raw = lane_dict["TotalClustersRaw"]
        total_clusters_pf = lane_dict["TotalClustersPF"]
        lane_yield = lane_dict["Yield"]
        sample_demux_results = lane_dict["DemuxResults"]

        for sample_demux_result in sample_demux_results:
            # In samplesheets from the lims this is prefixed with the project name
            sample_name = sample_demux_result["SampleName"]
            sample_tag = fetch_index_sequence(sample_demux_result)
            sample_project = samplesheet.project_for_sample(sample_name, lane_nbr)
            sample_library_name = samplesheet.library_name_for_sample(sample_name, lane_nbr)

            read_metrics = sample_demux_result["ReadMetrics"]
            for read_metric in read_metrics:
                read_nbr = read_metric["ReadNumber"]
                cycles = reads_and_cycles[read_nbr]
                read_yield = read_metric["Yield"]
                mean_q = read_metric["QualityScoreSum"] / read_metric["Yield"]
                percent_q_30 = read_metric["YieldQ30"] / read_metric["Yield"]
                percent_of_lane = read_yield / lane_yield

                # TODO Is this good enough? The cluster info is only easily available on a per lane basis
                pass_filter_clusters = total_clusters_pf * percent_of_lane
                percent_tag_error = "TODO" # TODO How to get this?
                yield SampleLevelStats(flowcell_name, sample_project, sample_name, sample_tag, lane_nbr, read_nbr, cycles,
                                       percent_of_lane, pass_filter_clusters, percent_q_30, percent_tag_error,
                                       sample_library_name, mean_q)

                # TODO How to get % tag error?
                #   -> Not sure what the definition of this actually is?
                # TODO How to get library name?
                #   -> This probably needs to be picked up from the sample sheet
