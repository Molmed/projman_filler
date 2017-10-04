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

If the runfolder has already been loaded into the database, you can add the make projman_filler remove the old results
and add the new once by adding the `--force` flag, i.e:

```
export PROJMAN_DB="mssql+pymssql://<usernam>:<password>@<host>/<database>?charset=utf8"
projman_filler --force <path to the runfolder to insert into the db>
```

Install instructions
--------------------

Installing projman_filler can be done in the following way.

Please note that projman_filler uses python3, ensure this is first installed.

 * Download the projman_filler delivery
 * Install the Illumina Interop library (this has to be done "manually" with pip): `pip install -f https://github.com/Illumina/interop/releases/latest interop`
 * Install projman_filler with pip as normal, i.e. `pip install <path to source>`

Credits
---------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage

