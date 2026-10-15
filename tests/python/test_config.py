import json
import sys
import unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
from validate_config import validate

class TestConfiguration(unittest.TestCase):
    def setUp(self):
        self.cfg = json.loads((ROOT / "examples/single_pipe.json").read_text())

    def test_example(self):
        self.assertIs(validate(self.cfg), self.cfg)

    def test_missing_section(self):
        self.cfg.pop("boundary")
        with self.assertRaises(ValueError):
            validate(self.cfg)

    def test_bad_theta(self):
        self.cfg["numerics"]["theta"] = 1.1
        with self.assertRaises(ValueError):
            validate(self.cfg)

    def test_nonfinite_pressure(self):
        self.cfg["boundary"]["inlet_pressure_Pa"] = float("nan")
        with self.assertRaises(ValueError):
            validate(self.cfg)

    def test_boolean_cells_rejected(self):
        self.cfg["pipeline"]["cells"] = True
        with self.assertRaises(ValueError):
            validate(self.cfg)

if __name__ == "__main__":
    unittest.main()
