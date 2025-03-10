ProjMan Filler
==============

App for filling ProjMan


* Free software: GNU General Public License v3


Features
--------

* This program will fill the ProjMan database, based on the `Stats.json` and `Interop` files in an Illumina runfolder

Usage
-----

To run projman_filler run:

```
export PROJMAN_DB="mssql+pymssql://<usernam>:<password>@<host>/<database>?charset=utf8"
projman_filler <path to the runfolder to insert into the db>
```

The location of the bcl2fastq statistics folder, i.e. where Stats.json is located, can be specified using the `-b|--bcl2fastq-stats` flag.
Specify the location relative to the runfolder. Default: "Unaligned/Stats"

If the runfolder has already been loaded into the database, you can add the make projman_filler remove the old results
and add the new once by adding the `--force` flag, i.e:

```
export PROJMAN_DB="mssql+pymssql://<usernam>:<password>@<host>/<database>?charset=utf8"
projman_filler --force <path to the runfolder to insert into the db>
```
For more flags, see cli.py.

Install instructions
--------------------

Installing projman_filler can be done in the following way.

Please note that projman_filler uses python 3.10, ensure this is first installed.

 * Download the projman_filler delivery
 * Install projman_filler with pip as normal, i.e. `pip install <path to source>`


Development and testing
-----------------------
Set up a Python environment using the latest version of Python 3.10 (3.10.1 as of this writing).

Then run:
```
 pip install -e . -r requirements_dev.txt
```
To run the app locally:
```
export PROJMAN_DB="mssql+pymssql://<usernam>:<password>@<host>/<database>?charset=utf8"
&& python ./projman_filler/cli.py <path to the runfolder to insert into the db>
```

Run tests with:
```
 python setup.py pytest
```

Credits
---------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage

