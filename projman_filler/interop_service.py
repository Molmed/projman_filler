
from collections import defaultdict

from interop import py_interop_run_metrics, py_interop_run, py_interop_summary


class InteropService(object):

    def __init__(self, runfolder):
        run_metrics = py_interop_run_metrics.run_metrics()
        run_metrics.run_info()

        valid_to_load = py_interop_run.uchar_vector(py_interop_run.MetricCount, 0)
        py_interop_run_metrics.list_summary_metrics_to_load(valid_to_load)
        run_metrics.read(runfolder, valid_to_load)

        summary = py_interop_summary.run_summary()
        py_interop_summary.summarize_run_metrics(run_metrics, summary)

        self.summary = summary

    def get_non_index_reads(self):
        # First pick-out the reads which are not index reads
        non_index_reads = []
        for read_nbr in range(self.summary.size()):
            if not self.summary.at(read_nbr).read().is_index():
                non_index_reads.append(read_nbr)
        return non_index_reads

    def get_error_rates(self):
        lanes = self.summary.lane_count()

        non_index_reads = self.get_non_index_reads()

        error_rates = defaultdict(dict)
        # "Renumber" the reads to have their intuitive read numbers, i.e. Read 1 -> 1, Read 2 -> 2
        for new_nbr, old_read in enumerate(non_index_reads):
            for lane in range(lanes):
                read = self.summary.at(old_read).at(lane)
                error_rate = read.error_rate().mean()
                error_rates[lane+1][new_nbr+1] = error_rate
        return error_rates

    def get_densities(self):
        lanes = self.summary.lane_count()

        non_index_reads = self.get_non_index_reads()

        densities = defaultdict(dict)
        # "Renumber" the reads to have their intuitive read numbers, i.e. Read 1 -> 1, Read 2 -> 2
        for new_nbr, original_read_nbr in enumerate(non_index_reads):
            for lane in range(lanes):
                read = self.summary.at(original_read_nbr).at(lane)
                cluster_density = read.density().mean()
                pass_filter_density = read.density_pf().mean()
                densities[lane+1][new_nbr+1] = \
                    {"raw_density": cluster_density, "pass_filter_density": pass_filter_density}
        return densities

    def get_q30(self):
        lanes = self.summary.lane_count()

        non_index_reads = self.get_non_index_reads()

        q30s = defaultdict(dict)
        # "Renumber" the reads to have their intuitive read numbers, i.e. Read 1 -> 1, Read 2 -> 2
        for new_nbr, original_read_nbr in enumerate(non_index_reads):
            for lane in range(lanes):
                q30 = self.summary.at(original_read_nbr).at(lane).percent_gt_q30()
                q30s[lane+1][new_nbr+1] = q30
        return q30s

    def get_q_scores(self):
        lanes = self.summary.lane_count()

        non_index_reads = self.get_non_index_reads()

        q_scores = defaultdict(dict)
        # "Renumber" the reads to have their intuitive read numbers, i.e. Read 1 -> 1, Read 2 -> 2
        for new_nbr, original_read_nbr in enumerate(non_index_reads):
            for lane in range(lanes):
                q = self.summary.at(original_read_nbr).at(lane).q_metric()
                q_scores[lane+1][new_nbr+1] = q
        return q_scores
