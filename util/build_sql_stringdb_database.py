#!/usr/bin/python3

import argparse

from sql import chunked_file, build_sql_stringdb_database

def main(args):
    output_file = "{}/{}".format(args.output_dir, args.output_file)
    
    build_sql_stringdb_database(args.alias_file, args.evidence_file,
            args.actions_file, output_file, verbose=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build a database for SDBL")
    parser.add_argument("alias_file", help="example: 10090.protein.aliases.v10.5.txt.gz")
    parser.add_argument("evidence_file", help="example: 10090.protein.links.detailed.v10.5.txt.gz")
    parser.add_argument("actions_file", help="example: 10090.protein.actions.v10.5.txt.gz")
    parser.add_argument("output_file", help="database file to output")
    parser.add_argument("--output_dir", default=".", help="directory to write output into")

    args = parser.parse_args()
    main(args)
