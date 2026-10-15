#!/usr/bin/env python3
"""E1 dependency-free subset validator; schema.json is the authoritative contract."""
import json
import math
import sys
from pathlib import Path

def validate(obj):
    if not isinstance(obj, dict):
        raise ValueError("Configuration must be an object")
    if set(obj) - {"schema_version", "pipeline", "numerics", "boundary", "metadata"}:
        raise ValueError("Unknown top-level configuration keys")
    if obj.get("schema_version") != "1.0":
        raise ValueError("Unsupported schema_version")
    rules = {
        "pipeline": {"length_m": (float, 0), "diameter_m": (float, 0),
                     "cells": (int, 2), "temperature_K": (float, 0),
                     "compressibility_Z": (float, 0)},
        "numerics": {"dt_s": (float, 0), "steps": (int, 0),
                     "theta": (float, 0)},
        "boundary": {"inlet_pressure_Pa": (float, 0),
                     "outlet_mass_flow_kg_s": (float, 0)}
    }
    for group, fields in rules.items():
        values = obj.get(group)
        if not isinstance(values, dict) or set(values) != set(fields):
            raise ValueError(f"{group} must contain exactly {sorted(fields)}")
        for key, (kind, minimum) in fields.items():
            v = values[key]
            if kind is int:
                if type(v) is not int or v < minimum:
                    raise ValueError(f"Invalid {group}.{key}")
            elif type(v) not in (float, int) or not math.isfinite(v) or (
                v <= minimum if key not in ("outlet_mass_flow_kg_s", "theta") else v < minimum
            ):
                raise ValueError(f"Invalid {group}.{key}")
    if obj["numerics"]["theta"] > 1:
        raise ValueError("theta must not exceed 1")
    if "metadata" in obj and not isinstance(obj["metadata"], dict):
        raise ValueError("metadata must be an object")
    return obj

if __name__ == "__main__":
    try:
        for filename in sys.argv[1:]:
            validate(json.loads(Path(filename).read_text(encoding="utf-8")))
            print(f"PASS {filename}")
        if len(sys.argv) == 1:
            raise ValueError("Usage: validate_config.py FILE.json [...]")
    except (ValueError, OSError, json.JSONDecodeError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        sys.exit(1)
