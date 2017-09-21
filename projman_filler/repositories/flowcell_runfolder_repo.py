
from projman_filler.models.db_models import FlowcellRunfolder


class FlowcellRunfolderRepo(object):

    def __init__(self, session_factory):
        self._session_factory = session_factory

    def add(self, flowcell_runfolder):
        session = self._session_factory()
        if isinstance(flowcell_runfolder, list):
            session.add_all(flowcell_runfolder)
        else:
            session.add(flowcell_runfolder)
        session.commit()

    def contains_flowcell(self, flowcell_name):
        session = self._session_factory()
        q = session.query(FlowcellRunfolder).filter(FlowcellRunfolder.flowcell_id == flowcell_name)
        return session.query(q.exists())
