
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
            for read_nbr, read_info in enumerate(stats_json_read_info):
                if not read_info["IsIndexedRead"]:
                    non_index_reads.append(read_nbr)
        else:
            non_index_reads = []
            for read_nbr in range(interop_summary.size()):
                if not interop_summary.at(read_nbr).read().is_index():
                    non_index_reads.append(read_nbr)
        return non_index_reads
