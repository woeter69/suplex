import math
import pathlib
import sys

import pysuplex


source_root = pathlib.Path(sys.argv[1])
solver = pysuplex.Solver()
solver.read_lp(str(source_root / "data" / "examples" / "tiny.lp"))
solver.set_algorithm(pysuplex.Algorithm.PRIMAL_SIMPLEX)
status = solver.solve()

assert status == pysuplex.Status.OPTIMAL
assert math.isclose(solver.objective, 3.0, rel_tol=0.0, abs_tol=1e-7)
assert solver.solution.shape == (2,)
assert solver.dual_values.shape == (2,)
