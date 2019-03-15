
import os
import json

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, scoped_session

from projman_filler.samplesheet import Samplesheet
from projman_filler.interop_service import InteropService
from projman_filler.lane_level_statistics import calculate_lane_statistics
from projman_filler.sample_level_statistics import calculate_sample_statistics
from projman_filler.exceptions import FlowcellAlreadyInDb
from projman_filler.models.db_models import FlowcellRunfolder

from projman_filler.repositories.sample_results_repo import SampleResultRepo
from projman_filler.repositories.flowcell_lane_results_repo import FlowcellLaneResultsRepo
from projman_filler.repositories.flowcell_runfolder_repo import FlowcellRunfolderRepo


class ReadNumberMapper():

    def __init__(self, stats_json, atac_seq_mode=False):
        self._stats_json = stats_json
        self._atac_seq_mode = atac_seq_mode

    def get_reads_and_cycles(self):
        reads_and_cycles = {}
        for read_info in self._stats_json["ReadInfosForLanes"][0]["ReadInfos"]:
            if not read_info["IsIndexedRead"]:
                reads_and_cycles[read_info["Number"]] = read_info["NumCycles"]
        return reads_and_cycles

    def get_non_index_reads(self, interop_summary):
        non_index_reads = []
        if self._atac_seq_mode:
            stats_json_read_info = self._stats_json["ReadInfosForLanes"][0]["ReadInfos"]
            for index, read_info in enumerate(stats_json_read_info):
                if not read_info["IsIndexedRead"]:
                    non_index_reads.append(index)
        else:
            non_index_reads = []
            for read_nbr in range(interop_summary.size()):
                if not interop_summary.at(read_nbr).read().is_index():
                    non_index_reads.append(read_nbr)
        return non_index_reads


class App(object):

    def __init__(self, db_connection_string, debug=False, flowcell_runfolder_repo=None,
                 flowcell_lane_results_repo=None, sample_results_repo=None):

        engine = create_engine(db_connection_string, echo=debug)
        session_factory = scoped_session(sessionmaker())
        session_factory.configure(bind=engine)

        if not flowcell_lane_results_repo:
            self.flowcell_lane_results_repo = FlowcellLaneResultsRepo(session_factory)

        if not flowcell_runfolder_repo:
            self.flowcell_runfolder_repo = FlowcellRunfolderRepo(session_factory)

        if not sample_results_repo:
            self.sample_results_repo = SampleResultRepo(session_factory)

    def parse_rundate_from_runfolder_name(self, runfolder_name):
        return os.path.basename(runfolder_name).split("_")[0]

    def insert_runfolder_into_db(self, runfolder, bcl2fastq_stats, force=False, atac_seq_mode=False):

        stats_json = os.path.join(runfolder, bcl2fastq_stats, 'Stats.json')

        with open(stats_json, 'r') as f:
            stats_json = json.load(f)

        flowcell_name = stats_json["Flowcell"]
        conversion_results = stats_json["ConversionResults"]

        if self.flowcell_runfolder_repo.contains_flowcell(flowcell_name):
            if force:
                print("Found the specified runfolder in the db, but got a force option, so will proceed to "
                      "delete it and insert new values.")
                self.flowcell_lane_results_repo.delete_by_flowcell_name(flowcell_name)
                self.flowcell_runfolder_repo.delete_by_flowcell_name(flowcell_name)
                self.sample_results_repo.delete_by_flowcell_name(flowcell_name)
            else:
                raise FlowcellAlreadyInDb

        if atac_seq_mode:
            print("ATAC-seq mode activated. Will re-map read numbers according to settings used by bcl2fastq.")

        read_mapper = ReadNumberMapper(stats_json, atac_seq_mode)

        interop = InteropService(runfolder, read_mapper)
        densities = interop.get_densities()
        error_rates = interop.get_error_rates()
        q30s = interop.get_q30()

        reads_and_cycles = read_mapper.get_reads_and_cycles()

        lane_stats = calculate_lane_statistics(flowcell_name, conversion_results, reads_and_cycles,
                                               error_rates, densities, q30s)

        self.flowcell_lane_results_repo.add(list(lane_stats))

        samplesheet_file = os.path.join(runfolder, "SampleSheet.csv")
        samplesheet = Samplesheet(samplesheet_file)

        sample_stats = calculate_sample_statistics(flowcell_name, conversion_results, reads_and_cycles, samplesheet)
        self.sample_results_repo.add(list(sample_stats))

        runfolder_name = os.path.basename(runfolder)
        runfolder_date = self.parse_rundate_from_runfolder_name(runfolder)
        flowcell_runfolder = FlowcellRunfolder(flowcell_id=flowcell_name,
                                               runfolder_name=runfolder_name,
                                               run_date=runfolder_date)
        self.flowcell_runfolder_repo.add(flowcell_runfolder)
