class RunStatsParserInterface:
    def get_flowcell_name(self) -> str:
        """Returns the flowcell name."""
        pass

    def get_non_index_reads(self) -> list:
        """Returns list of indices of non-indexed reads."""
        pass

    def get_conversion_results(self) -> list:
        """Returns list of Lane objects"""
        pass
