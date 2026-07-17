from __future__ import annotations

import argparse

from .pipeline import run_preprocessing


def main() -> None:
    """Parse the command line and start the requested workflow."""
    # The subcommand structure leaves room for additional workflow stages later.
    parser = argparse.ArgumentParser(description="Perturb-seq preprocessing")
    sub = parser.add_subparsers(dest="command", required=True)
    preprocess = sub.add_parser("preprocess")
    preprocess.add_argument("config")
    args = parser.parse_args()
    if args.command == "preprocess":
        run_preprocessing(args.config)
