from __future__ import annotations

import argparse

from .pipeline import run_preprocessing


def main() -> None:
    parser = argparse.ArgumentParser(description="Perturb-seq preprocessing")
    sub = parser.add_subparsers(dest="command", required=True)
    preprocess = sub.add_parser("preprocess")
    preprocess.add_argument("config")
    args = parser.parse_args()
    if args.command == "preprocess":
        run_preprocessing(args.config)
