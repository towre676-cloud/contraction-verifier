import subprocess, sys

PYEXE = sys.executable  # use the pytest interpreter (your venv)

def run_ok(args):
    r = subprocess.run([PYEXE] + args, capture_output=True, text=True)
    assert r.returncode == 0, (
        "cmd failed: {}\nstdout:\n{}\nstderr:\n{}".format(
            " ".join([PYEXE] + args), r.stdout, r.stderr
        )
    )

def test_core():
    run_ok(["-m","contraction_verifier.cli",
            "examples/core_banach/example_core.json","--type","core_banach"])

def test_ledger():
    run_ok(["-m","contraction_verifier.cli",
            "examples/ledger/example_ledger.json","--type","ledger"])

def test_coupled():
    run_ok(["-m","contraction_verifier.cli",
            "examples/coupled_streams/example_coupled.json","--type","coupled_stream"])

def test_gammaV():
    run_ok(["-m","contraction_verifier.cli",
            "examples/gamma_v/example_gammaV.json","--type","core_banach","--use-gammaV","--svd-exact"])
