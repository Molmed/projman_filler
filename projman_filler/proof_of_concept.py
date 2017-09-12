
from collections import defaultdict
import json
import os
from pprint import pprint

import xmltodict
from interop import py_interop_run_metrics, py_interop_run, py_interop_summary, py_interop_plot


# -------------------------
# The old script does this
# -------------------------

# Notes
#  - Sisyphus imports data from summary/<>/report.xml

RUNFOLDER = "data/170908_ST-E00216_0211_BH3JTGCCXY"
run_info_xml = os.path.join(RUNFOLDER, "RunInfo.xml")
demultiplexing_stats_xml_file = os.path.join(RUNFOLDER, "Unaligned", "Stats", "DemultiplexingStats.xml")

with open(run_info_xml) as f:
    run_info = xmltodict.parse(f.read())

#with open(demultiplexing_stats_xml_file) as f:
#    demultiplexing_stats = xmltodict.parse(f.read())
#
#all_projects = list(filter(lambda x: not (x["@name"] == 'default' or x["@name"] == 'all'),
#                           demultiplexing_stats["Stats"]["Flowcell"]["Project"]))

stats_json = os.path.join(RUNFOLDER, 'Unaligned', 'Stats', 'Stats.json')
with open(stats_json, 'r') as f:
    stats_json = json.load(f)

flowcell_name = stats_json["Flowcell"]
# - Check if flowcell id is already in database
# eval{	$dbhProj->do(qq(
# INSERT INTO flowcell_runfolder(flowcell_id, runfolder_name, run_date)
# VALUES('$fcId', '$rfName', '$runDate')));
# };

conversion_results = stats_json["ConversionResults"]

reads_and_cycles = {}
for read_info in stats_json["ReadInfosForLanes"][0]["ReadInfos"]:
    if not read_info["IsIndexedRead"]:
        reads_and_cycles[read_info["Number"]] = read_info["NumCycles"]


# - Insert lane result statistics
# - For each lane
#   - For each read
#	    $dbh->do(qq(
# INSERT INTO flowcell_lane_results(flowcell_id, lane_num, read_num, raw_density, pf_density,
#                                   error_rate, raw_clusters, pf_clusters, cycles, pct_q30, mean_q)
# VALUES('$fcId', $lane->{Id}, $read->{Id}, $read->{DensityRaw}, $read->{DensityPF},
#        $read->{ErrRate}, $read->{Raw}, $read->{PF}, $read->{Cycles}, $read->{PctQ30}, $read->{AvgQ}))
#);

# - Insert sample result statistics
# - For each sample
#   - For each tag
#       - For each lane
#           - For each read
#$dbh->do(qq(INSERT
#INTO
#sample_results(flowcell_id, project_id, sample_name, tag_seq, lane_num, read_num, \
#               cycles, pct_lane, pf_clusters, pct_q30, pct_tag_err, library_name, mean_q) \
#    VALUES('$fcId', '$proj', '$sample->{Id}', '$tag->{Id}', $lane->{Id}, $read->{Id}, \
#    $read->{Cycles}, $read->{PctLane}, $read->{PF}, $read->{PctQ30}, \
#    $read->{TagErr}, '$read->{LibraryName}', $read->{AvgQ})
#)
#);




summary = load_summary(RUNFOLDER)
densities = get_densities(summary)
error_rates = get_error_rates(summary)
q30s = get_q30(summary)
#q_scores = get_q_scores(summary)
#calculate_lane_statistics(conversion_results, error_rates, densities, q30s)

calculate_sample_statistics()

