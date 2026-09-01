from click.testing import CliRunner

from inSVert import __version__
from inSVert.cli import cli
from inSVert.utils_sim import buildheader


def test_cli_and_vcf_header_report_package_version():
    result = CliRunner().invoke(cli, ["--version"])

    assert result.exit_code == 0
    assert result.output == f"inSVert, version {__version__}\n"
    assert f"##source=inSVert-{__version__}\n" in buildheader(["chr1"], [100])
