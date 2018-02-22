
from collections import defaultdict

from projman_filler.models.db_models import FlowcellLaneResult

def _get_mean_q_scores(lane_dict):

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
            # This will be given as a fraction to be compatible with legacy data /JD 2017-10-04
            percent_q30 = float(q30s[lane_nbr][read_nbr]) / 100
            mean_q_for_read = mean_q[read_nbr]
            yield FlowcellLaneResult(flowcell_id=flowcell_name, lane_num=lane_nbr, read_num=read_nbr,
                                     raw_density=raw_density, pf_density=pf_density, error_rate=error_rate,
                                     raw_clusters=total_clusters_raw, pf_clusters=total_clusters_pf,
                                     cycles=cycles, pct_q30=percent_q30, mean_q=mean_q_for_read)

