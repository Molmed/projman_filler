
from projman_filler.models.db_models import SampleResult

def sum_of_mismatch_counts(mismatch_counts):
    result = 0
    for k, v in mismatch_counts.items():
        if int(k) != 0:
            result = result + v
    return result


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

            sample_project = samplesheet.project_for_sample(sample_name, lane_nbr)
            sample_library_name = samplesheet.library_name_for_sample(sample_name, lane_nbr)

            number_of_reads = sample_demux_result["NumberReads"]
            sample_yield = float(sample_demux_result["Yield"])
            fraction_of_lane = sample_yield / lane_yield
            percent_of_lane = fraction_of_lane * 100

            read_metrics = sample_demux_result["ReadMetrics"]

            index_metrics = sample_demux_result["IndexMetrics"]
            for index in index_metrics:
                sample_tag = index["IndexSequence"]
                mismatch_counts = sum_of_mismatch_counts(index["MismatchCounts"])
                percent_tag_error = (float(mismatch_counts) / float(number_of_reads))*100
                for read_metric in read_metrics:
                    read_nbr = read_metric["ReadNumber"]
                    cycles = reads_and_cycles[read_nbr]
                    mean_q = read_metric["QualityScoreSum"] / read_metric["Yield"]
                    percent_q_30 = (float(read_metric["YieldQ30"]) / read_metric["Yield"])*100
                    pass_filter_clusters = total_clusters_pf * fraction_of_lane
                    yield SampleResult(flowcell_id=flowcell_name, project_id=sample_project, sample_name=sample_name,
                                       tag_seq=sample_tag, lane_num=lane_nbr, read_num=read_nbr, cycles=cycles,
                                       pct_lane=percent_of_lane, pf_clusters=pass_filter_clusters, pct_q30=percent_q_30,
                                       pct_tag_err=percent_tag_error, library_name=sample_library_name,
                                       mean_q=mean_q)

