"""Run numerical and Qt regression tests in fresh Python processes."""
from concurrent.futures import ThreadPoolExecutor
import os
from pathlib import Path
import subprocess
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
TESTS = Path(__file__).resolve().parent
GUI_SHARDS = 3
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('QT_QPA_PLATFORM', 'offscreen')


def iter_cases(suite):
    for item in suite:
        if isinstance(item, unittest.TestSuite):
            yield from iter_cases(item)
        else:
            yield item


def run_suite(shard):
    sys.path.insert(0, str(ROOT / 'bin'))
    # Initialize the XLSX dependency before Qt starts its export worker.
    import openpyxl

    loader = unittest.defaultTestLoader
    if shard == 'numeric':
        modules = [path.stem for path in sorted(TESTS.glob('test_*.py')) if path.stem != 'test_gui']
        suite = loader.loadTestsFromNames(modules)
    else:
        tests = list(iter_cases(loader.loadTestsFromName('test_gui')))
        suite = unittest.TestSuite(tests[int(shard)::GUI_SHARDS])
    return unittest.TextTestRunner(verbosity=2).run(suite).wasSuccessful()


def run_child(shard):
    completed = subprocess.run([sys.executable, str(Path(__file__).resolve()), shard],
                               cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               env={**os.environ, 'PYTHONIOENCODING': 'utf-8'},
                               text=True, encoding='utf-8', timeout=240)
    return completed.returncode, completed.stdout


def run_all():
    # Long runs in a single offscreen Qt process can retain widget state between tests.
    shards = ('numeric', *(str(index) for index in range(GUI_SHARDS)))
    with ThreadPoolExecutor(max_workers=len(shards)) as pool:
        results = [pool.submit(run_child, shard) for shard in shards]
        failed = False
        for shard, future in zip(shards, results):
            try:
                status, output = future.result()
            except subprocess.TimeoutExpired as error:
                status = 1
                output = (error.output or b'').decode('utf-8', errors='replace') + '\nTimed out.\n'
            print(f'[{shard}]\n{output}', end='', flush=True)
            failed |= status != 0
    return not failed


if __name__ == '__main__':
    sys.exit(not (run_suite(sys.argv[1]) if len(sys.argv) > 1 else run_all()))
