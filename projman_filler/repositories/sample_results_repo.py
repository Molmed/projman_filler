
from projman_filler.models.db_models import SampleResult


class SampleResultRepo(object):

    def __init__(self, session_factory):
        self._session_factory = session_factory

    def add(self, sample_results):
        session = self._session_factory()
        if isinstance(sample_results, list):
            session.add_all(sample_results)
        else:
            session.add(sample_results)
        session.commit()
