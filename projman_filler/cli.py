
"""Console script for projman_filler."""

import click

from projman_filler.stats_service import calculate_stats

@click.command("projman_filler")
@click.argument('runfolder', type=click.Path())
def main(runfolder):
    """Console script for projman_filler."""
    calculate_stats(runfolder)

if __name__ == "__main__":
    main()
