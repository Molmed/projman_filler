
import os
import unittest
from unittest.mock import patch
from click.testing import CliRunner

from projman_filler.cli import main


class TestCLI(unittest.TestCase):

    @patch("projman_filler.app.App.insert_olink_runfolder_into_db")
    def test_cli_prefers_olink_mode_even_if_demultiplexer_specified(
        self, mock_insert_olink
    ):
        runner = CliRunner()

        os.environ["PROJMAN_DB"] = "sqlite:///:memory:"

        result = runner.invoke(
            main,
            [
                "--demultiplexer", "bcl2fastq",
                "--olink-mode",
                "/path/to/runfolder"
            ],
            catch_exceptions=False,
        )

        self.assertEqual(result.exit_code, 0)

        expected_warning = (
            "Warning: Both olink_mode and demultiplexer specified, olink_mode take precedence"
        )
        self.assertIn(expected_warning, result.output)

        mock_insert_olink.assert_called_once_with("/path/to/runfolder", False)


    def test_no_olink_no_demultiplexer(
        self
    ):
        runner = CliRunner()

        os.environ["PROJMAN_DB"] = "sqlite:///:memory:"

        with self.assertRaises(ValueError) as ctx:
            runner.invoke(
                main,
                [
                    "/path/to/runfolder"
                ],
                catch_exceptions=False
            )

        self.assertEqual(str(ctx.exception), "None is not a supported demultiplexer")

