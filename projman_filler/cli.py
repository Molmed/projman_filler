
"""Console script for projman_filler."""
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, scoped_session

import click
import os
import sys

from projman_filler.app import App
from projman_filler.exceptions import FlowcellAlreadyInDb

from projman_filler import __version__ as projman_filler_version


@click.command("projman_filler")
@click.option('--force', is_flag=True)
@click.option('--debug', is_flag=True)
@click.option('--atac-seq-mode', is_flag=True)
@click.option('-b', '--bcl2fastq-stats', default="Unaligned/Stats", type=click.Path())
@click.argument('runfolder', type=click.Path())
def main(runfolder, force, atac_seq_mode, bcl2fastq_stats, debug):
    """Console script for projman_filler."""
    print("projman_filler v{}".format(projman_filler_version))
    try:
        db_connection_string = os.environ["PROJMAN_DB"]
        app = App(db_connection_string, debug)
        app.insert_runfolder_into_db(runfolder, bcl2fastq_stats, force=force, atac_seq_mode=atac_seq_mode)
    except FlowcellAlreadyInDb:
        print("Flowcell was already present in db.")
        sys.exit(1)


if __name__ == "__main__":
    main()
