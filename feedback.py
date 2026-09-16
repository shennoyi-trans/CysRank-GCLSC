"""Validate and archive experimental or simulation feedback without mixing labels."""
import argparse
from src.paths import root_path
from src.pipeline import ingest

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', required=True, type=root_path)
    parser.add_argument('--output', required=True, type=root_path)
    ingest(parser.parse_args())
