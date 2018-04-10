
import copy

from projman_filler.models.db_models import SampleResult


def sum_of_mismatch_counts(mismatch_counts):
    result = 0
    for k, v in mismatch_counts.items():
        if int(k) != 0:
            result = result + v
    return result


def _add_tag_and_tag_error(sample_result, sample_demux_result):
    sample_result_copy = copy.deepcopy(sample_result)
    if sample_demux_result.get("IndexMetrics"):
        index_metrics = sample_demux_result["IndexMetrics"]
        number_of_reads = sample_demux_result["NumberReads"]
        for index in index_metrics:
            # Replace + with - to be compatible with legacy data in ProjMan db. /JD 2017-10-04
            sample_result_copy.tag_seq = index["IndexSequence"].replace("+", "-")
            mismatch_counts = sum_of_mismatch_counts(index["MismatchCounts"])
            if number_of_reads != 0:
                sample_result_copy.pct_tag_err = (float(mismatch_counts) / float(number_of_reads))*100
            else:
                sample_result_copy.pct_tag_err = None
            yield sample_result_copy
    else:
        sample_result_copy.tag_seq = "unknown"
        sample_result_copy.pct_tag_err = None
        yield sample_result_copy


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
            sample_id = sample_demux_result["SampleId"]

            sample_project = samplesheet.project_for_sample(sample_id, lane_nbr)
            sample_library_name = samplesheet.library_name_for_sample(sample_id, lane_nbr)

            sample_yield = float(sample_demux_result["Yield"])
            if lane_yield == 0:
                fraction_of_lane = 0
                percent_of_lane = 0
            else:
                fraction_of_lane = sample_yield / lane_yield
                percent_of_lane = fraction_of_lane * 100

            read_metrics = sample_demux_result["ReadMetrics"]

            sample_result = SampleResult(flowcell_id=flowcell_name, project_id=sample_project, sample_name=sample_name,
                                         lane_num=lane_nbr,  pct_lane=percent_of_lane, library_name=sample_library_name)

            sample_results_with_index = _add_tag_and_tag_error(sample_result, sample_demux_result)

            for sample_result_with_index in sample_results_with_index:
                for read_metric in read_metrics:
                    sample_result_with_index_copy = copy.deepcopy(sample_result_with_index)
                    read_nbr = read_metric["ReadNumber"]
                    sample_result_with_index_copy.read_num = read_nbr
                    sample_result_with_index_copy.cycles = reads_and_cycles[read_nbr]

                    if read_metric["Yield"] == 0:
                        sample_result_with_index_copy.mean_q = None
                        sample_result_with_index_copy.pct_q30 = None
                    else:
                        sample_result_with_index_copy.mean_q = read_metric["QualityScoreSum"] / read_metric["Yield"]
                        sample_result_with_index_copy.pct_q30 = (float(read_metric["YieldQ30"]) / read_metric["Yield"])*100

                    sample_result_with_index_copy.pf_clusters = total_clusters_pf * fraction_of_lane
                    yield sample_result_with_index_copy

