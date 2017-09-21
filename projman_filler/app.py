
import os
import json

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, scoped_session

from projman_filler.samplesheet import Samplesheet
from projman_filler.interop_service import InteropService
from projman_filler.lane_level_statistics import calculate_lane_statistics
from projman_filler.sample_level_statistics import calculate_sample_statistics
from projman_filler.exceptions import FlowcellAlreadyInDb

from projman_filler.repositories.sample_results_repo import SampleResultRepo
from projman_filler.repositories.flowcell_lane_results_repo import FlowcellLaneResultsRepo
from projman_filler.repositories.flowcell_runfolder_repo import FlowcellRunfolderRepo


class App(object):

    def __init__(self, db_connection_string, flowcell_runfolder_repo=None,
                 flowcell_lane_results_repo=None, sample_results_repo=None):

        engine = create_engine(db_connection_string, echo=True)
        session_factory = scoped_session(sessionmaker())
        session_factory.configure(bind=engine)

        if not flowcell_lane_results_repo:
            self.flowcell_lane_results_repo = FlowcellLaneResultsRepo(session_factory)

        if not flowcell_runfolder_repo:
            self.flowcell_runfolder_repo = FlowcellRunfolderRepo(session_factory)

        if not sample_results_repo:
            self.sample_results_repo = SampleResultRepo(session_factory)

    @staticmethod
    def get_reads_and_cycles(stats_json):
        reads_and_cycles = {}
        for read_info in stats_json["ReadInfosForLanes"][0]["ReadInfos"]:
            if not read_info["IsIndexedRead"]:
                reads_and_cycles[read_info["Number"]] = read_info["NumCycles"]
        return reads_and_cycles

    def insert_runfolder_into_db(self, runfolder):

        stats_json = os.path.join(runfolder, 'Unaligned', 'Stats', 'Stats.json')
        with open(stats_json, 'r') as f:
            stats_json = json.load(f)

        flowcell_name = stats_json["Flowcell"]
        conversion_results = stats_json["ConversionResults"]

        if self.flowcell_runfolder_repo.contains_flowcell(flowcell_name):
            raise FlowcellAlreadyInDb

        interop = InteropService(runfolder)
        densities = interop.get_densities()
        error_rates = interop.get_error_rates()
        q30s = interop.get_q30()
        reads_and_cycles = self.get_reads_and_cycles(stats_json)

        lane_stats = calculate_lane_statistics(flowcell_name, conversion_results, reads_and_cycles,
                                               error_rates, densities, q30s)

        self.flowcell_lane_results_repo.add(list(lane_stats))

        samplesheet_file = os.path.join(runfolder, "SampleSheet.csv")
        samplesheet = Samplesheet(samplesheet_file)

        sample_stats = calculate_sample_statistics(flowcell_name, conversion_results, reads_and_cycles, samplesheet)
        self.sample_results_repo.add(list(sample_stats))
