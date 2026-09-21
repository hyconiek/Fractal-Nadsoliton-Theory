"""Run unittest cases and no-argument test functions without requiring pytest.

No pytest fixtures/parametrization are emulated. Unsupported signatures fail.
"""
import importlib.util
import inspect
from pathlib import Path
import sys
import unittest

path=Path(sys.argv[1])
spec=importlib.util.spec_from_file_location('intake_test_module',path)
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
suite=unittest.defaultTestLoader.loadTestsFromModule(module)
for name,obj in vars(module).items():
    if name.startswith('test_') and inspect.isfunction(obj):
        if inspect.signature(obj).parameters:
            raise RuntimeError('Unsupported fixture-dependent test: '+name)
        suite.addTest(unittest.FunctionTestCase(obj,description=name))
result=unittest.TextTestRunner(verbosity=2).run(suite)
raise SystemExit(not result.wasSuccessful())
