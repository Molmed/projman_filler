
import os
import json

from projman_filler.samplesheet import Samplesheet
from projman_filler.interop_service import InteropService
from projman_filler.lane_level_statistics import calculate_lane_statistics
from projman_filler.sample_level_statistics import calculate_sample_statistics

"""Main module."""


def calculate_stats(runfolder):

    stats_json = os.path.join(runfolder, 'Unaligned', 'Stats', 'Stats.json')
    with open(stats_json, 'r') as f:
        stats_json = json.load(f)

    reads_and_cycles = {}
    for read_info in stats_json["ReadInfosForLanes"][0]["ReadInfos"]:
        if not read_info["IsIndexedRead"]:
            reads_and_cycles[read_info["Number"]] = read_info["NumCycles"]

    flowcell_name = stats_json["Flowcell"]
    conversion_results = stats_json["ConversionResults"]

    interop = InteropService(runfolder)
    densities = interop.get_densities()
    error_rates = interop.get_error_rates()
    q30s = interop.get_q30()
    #q_scores = get_q_scores(summary)
    lane_stats = calculate_lane_statistics(flowcell_name, conversion_results, reads_and_cycles, error_rates, densities, q30s)
    for s in lane_stats:
        print(s)

    samplesheet_file = os.path.join(runfolder, "SampleSheet.csv")
    samplesheet = Samplesheet(samplesheet_file)

    sample_stats = calculate_sample_statistics(flowcell_name, conversion_results, reads_and_cycles, samplesheet)
    for s in sample_stats:
        print(s)
