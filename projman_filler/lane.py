class Lane:
    def __init__(self, number, total_clusters_raw, total_clusters_pf, mean_q_scores=None, yld=None, sample_demux_results=None):
        self._number = number
        self._total_clusters_raw = int(total_clusters_raw)
        self._total_clusters_pf = int(total_clusters_pf)
        self._mean_q_scores = mean_q_scores
        self._yield = yld
        self._sample_demux_results = sample_demux_results

    def get_stats(self):
        return self._number, \
               self._total_clusters_raw, \
               self._total_clusters_pf, \
               self._mean_q_scores, \
               self._yield, \
               self._sample_demux_results
