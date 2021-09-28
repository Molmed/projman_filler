import math
import sys

from projman_filler.models.db_models import FlowcellLaneResult


def calculate_lane_statistics(interop_service, flowcell_name, conversion_results):
    for lane in conversion_results:
        lane_nbr, \
        total_clusters_raw, \
        total_clusters_pf, \
        mean_q, \
        lane_yield, \
        sample_demux_results = lane.get_stats()

        for read_nbr in interop_service.get_reads():
            cycles = interop_service.get_cycles(read_nbr)
            error_rate = interop_service.get_error_rate(read_nbr, lane_nbr)
            raw_density = interop_service.get_raw_density(read_nbr, lane_nbr)
            pf_density = interop_service.get_pf_density(read_nbr, lane_nbr)
            percent_q30 = interop_service.get_pct_gt_q30(read_nbr, lane_nbr)
            yield FlowcellLaneResult(flowcell_id=flowcell_name, lane_num=lane_nbr, read_num=read_nbr,
                                     raw_density=raw_density, pf_density=pf_density, error_rate=error_rate,
                                     raw_clusters=total_clusters_raw, pf_clusters=total_clusters_pf,
                                     cycles=cycles, pct_q30=percent_q30)
