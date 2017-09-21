
"""Console script for projman_filler."""
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker, scoped_session

import click
import os
import sys

from projman_filler.app import App
from projman_filler.exceptions import FlowcellAlreadyInDb

@click.command("projman_filler")
@click.option('--force', is_flag=True)
@click.argument('runfolder', type=click.Path())
def main(runfolder, force):
    """Console script for projman_filler."""
    try:
        db_connection_string = os.environ["PROJMAN_DB"]
        app = App(db_connection_string)
        app.insert_runfolder_into_db(runfolder, force)
    except FlowcellAlreadyInDb:
        print("Flowcell was already present in db.")
        sys.exit(1)

if __name__ == "__main__":
    main()
