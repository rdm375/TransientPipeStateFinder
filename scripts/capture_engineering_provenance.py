#!/usr/bin/env python3
"""Capture verifiable source-tree provenance for engineering reports."""

import argparse
import json
import subprocess
from pathlib import Path


def git(repo, *args):
    try:
        result = subprocess.run(
            ["git", "-C", str(repo), *args],
            capture_output=True,
            text=True,
            check=True,
        )
        return result.stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return None


def inspect(repo):
    repo = Path(repo).resolve()

    revision = git(repo, "rev-parse", "HEAD")
    status = git(
        repo,
        "status",
        "--porcelain",
        "--untracked-files=all",
    )

    return {
        "path": str(repo),
        "revision": revision or "unknown",
        "working_tree_state": (
            "unknown" if status is None
            else "clean" if not status
            else "modified"
        ),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--statefinder",
        default=".",
    )
    parser.add_argument(
        "--solver",
        required=True,
    )
    parser.add_argument(
        "--output",
        required=True,
    )

    args = parser.parse_args()

    result = {
        "schema": "tpsf.build-provenance.v1",
        "statefinder": inspect(args.statefinder),
        "solver": inspect(args.solver),
        "validated_solver_baseline":
            "670ac493e4461e678ca66f2254dc3d402df3f312",
        "note": (
            "Source-tree inspection at capture time. "
            "This alone does not prove the exact binaries "
            "or object files used in an earlier build."
        ),
    }

    destination = Path(args.output)
    destination.parent.mkdir(parents=True, exist_ok=True)

    destination.write_text(
        json.dumps(result, indent=2) + "\n"
    )

    print(destination)


if __name__ == "__main__":
    main()
