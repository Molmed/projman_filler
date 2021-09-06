
import os

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, scoped_session

from projman_filler.samplesheet import Samplesheet
from projman_filler.lane_level_statistics import calculate_lane_statistics
from projman_filler.sample_level_statistics import calculate_sample_statistics
from projman_filler.exceptions import FlowcellAlreadyInDb
from projman_filler.models.db_models import FlowcellRunfolder
from projman_filler.bcl2fastq_run_stats_parser import Bcl2fastqRunStatsParser
from projman_filler.interop_run_stats_parser import InteropRunStatsParser

from projman_filler.repositories.sample_results_repo import SampleResultRepo
from projman_filler.repositories.flowcell_lane_results_repo import FlowcellLaneResultsRepo
from projman_filler.repositories.flowcell_runfolder_repo import FlowcellRunfolderRepo


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

    def delete_existing_flowcell_from_db(self, flowcell_name, force):
        if not self.flowcell_runfolder_repo.contains_flowcell(flowcell_name):
            return

        if force:
            print("Found the specified runfolder in the db, but got a force option, so will proceed to "
                  "delete it and insert new values.")
            self.flowcell_lane_results_repo.delete_by_flowcell_name(flowcell_name)
            self.flowcell_runfolder_repo.delete_by_flowcell_name(flowcell_name)
            self.sample_results_repo.delete_by_flowcell_name(flowcell_name)
        else:
            raise FlowcellAlreadyInDb

    def insert_flowcell_runfolder_into_db(self, runfolder, flowcell_name):
        runfolder_name = os.path.basename(runfolder)
        runfolder_date = self.parse_rundate_from_runfolder_name(runfolder)
        flowcell_runfolder = FlowcellRunfolder(flowcell_id=flowcell_name,
                                               runfolder_name=runfolder_name,
                                               run_date=runfolder_date)
        self.flowcell_runfolder_repo.add(flowcell_runfolder)

    def insert_runfolder_into_db(self, runfolder, bcl2fastq_stats_dir, force=False, atac_seq_mode=False, olink_mode=False):
        if olink_mode:
            print("Olink mode activated. Will read lane-level statistics from InterOp files instead of bcl2fastq Stats.json.")
            return self.insert_olink_runfolder_into_db(runfolder, force)

        bcl2fastq_stats = Bcl2fastqRunStatsParser(os.path.join(runfolder, bcl2fastq_stats_dir))
        flowcell_name = bcl2fastq_stats.get_flowcell_name()
        reads_and_cycles = bcl2fastq_stats.get_reads_and_cycles()
        conversion_results = bcl2fastq_stats.get_conversion_results()

        # Check if flowcell exists and should be overriden
        self.delete_existing_flowcell_from_db(flowcell_name, force)

        # For atac-seq we run bcl2fastq with special parameters declaring
        # that the second index should be interpreted as a non-index read.
        # So we allow overriding the Interop list of non-index-reads with
        # a custom list obtained from bcl2fastq stats. /ML 2021-09
        non_index_reads = None
        if atac_seq_mode:
            print("ATAC-seq mode activated. Will re-map read numbers according to settings used by bcl2fastq.")
            non_index_reads = bcl2fastq_stats.get_non_index_reads()

        interop = InteropRunStatsParser(runfolder, non_index_reads)
        lane_stats = calculate_lane_statistics(interop, flowcell_name, conversion_results)
        self.flowcell_lane_results_repo.add(list(lane_stats))

        samplesheet_file = os.path.join(runfolder, "SampleSheet.csv")
        samplesheet = Samplesheet(samplesheet_file)

        sample_stats = calculate_sample_statistics(flowcell_name, conversion_results, reads_and_cycles, samplesheet)
        self.sample_results_repo.add(list(sample_stats))

        self.insert_flowcell_runfolder_into_db(runfolder, flowcell_name)


    def insert_olink_runfolder_into_db(self, runfolder, force=False):
        interop = InteropRunStatsParser(runfolder)
        flowcell_name = interop.get_flowcell_name()

        # Check if flowcell exists and should be overriden
        self.delete_existing_flowcell_from_db(flowcell_name, force)

        conversion_results = interop.get_conversion_results()
        lane_stats = calculate_lane_statistics(interop, flowcell_name, conversion_results)

        self.flowcell_lane_results_repo.add(list(lane_stats))
        self.insert_flowcell_runfolder_into_db(runfolder, flowcell_name)
