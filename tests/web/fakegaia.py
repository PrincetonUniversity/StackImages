"""A fake Gaia TAP service for tests, answering from the recorded M13 results.

It applies the query's circle and magnitude limit to the fixture rows (so split queries
return exactly the stars a real service would), and can time out like GAVO's quick
endpoint for circles larger than `timeout_radius`.
"""

from __future__ import annotations

import math
import re
from pathlib import Path

FIXTURES = Path(__file__).parent / "fixtures"
TIMEOUT_BODY = (
    b"<?xml version='1.0'?><VOTABLE><RESOURCE type='results'><INFO name='QUERY_STATUS' value='ERROR'>"
    b"Field query: Query timed out (took too long).</INFO></RESOURCE></VOTABLE>"
)
_CIRCLE = re.compile(r"DISTANCE\(\s*([-\d.]+),\s*([-\d.]+),\s*gaia\.ra,\s*gaia\.dec\)\s*<\s*([-\d.]+)")
_MAG = re.compile(r"phot_g_mean_mag\s*<\s*([-\d.]+)")


class FakeGaia:
    def __init__(self, fixture: str = "m13_gavo.csv", timeout_radius: float | None = None) -> None:
        lines = (FIXTURES / fixture).read_text().splitlines()
        self.header = lines[0]
        columns = self.header.split(",")
        i_ra, i_dec, i_mag = columns.index("ra"), columns.index("dec"), columns.index("phot_g_mean_mag")
        self.rows = [
            (line, float(f[i_ra]), float(f[i_dec]), float(f[i_mag] or "nan"))
            for line in lines[1:]
            if (f := line.split(","))
        ]
        self.timeout_radius = timeout_radius
        self.queries: list[str] = []

    def answer(self, query: str) -> tuple[int, bytes]:
        """(HTTP status, body) for an ADQL query."""
        self.queries.append(query)
        circle, mag = _CIRCLE.search(query), _MAG.search(query)
        if not circle or not mag:
            return 400, b"unsupported query"
        ra0, dec0, radius = (float(v) for v in circle.groups())
        if self.timeout_radius is not None and radius > self.timeout_radius:
            return 500, TIMEOUT_BODY
        limit = float(mag.group(1))
        kept = [
            line for line, ra, dec, g in self.rows if g < limit and _distance(ra0, dec0, ra, dec) < radius
        ]
        return 200, ("\n".join([self.header, *kept]) + "\n").encode()


def _distance(ra0: float, dec0: float, ra: float, dec: float) -> float:
    a, b = math.radians(dec0), math.radians(dec)
    h = math.sin((b - a) / 2) ** 2 + math.cos(a) * math.cos(b) * math.sin(math.radians(ra - ra0) / 2) ** 2
    return math.degrees(2 * math.asin(math.sqrt(min(h, 1.0))))
