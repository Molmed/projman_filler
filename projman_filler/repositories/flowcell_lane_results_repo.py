
from projman_filler.models.db_models import FlowcellLaneResult

class FlowcellLaneResultsRepo(object):

    def __init__(self, session_factory):
        self._session_factory = session_factory

    def add(self, flowcell_lane_results):
        session = self._session_factory()
        if isinstance(flowcell_lane_results, list):
            session.add_all(flowcell_lane_results)
        else:
            session.add(flowcell_lane_results)
        session.commit()

    def delete_by_flowcell_name(self, flowcell_name):
        session = self._session_factory()
        session.query(FlowcellLaneResult).filter(FlowcellLaneResult.flowcell_id == flowcell_name).delete()
        session.commit()
