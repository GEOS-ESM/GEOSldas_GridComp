import os
import glob
import argparse
import subprocess
import multiprocessing
from datetime import datetime, timedelta

def parse_args():
    """
    No arguments: process yesterday only (the normal NRT cron behavior).
    --date: backfill a single specific day.
    --start [--end]: backfill an inclusive range of days.
    """
    parser = argparse.ArgumentParser(
        description="Process SMOS L1C data. With no arguments, processes "
                     "yesterday's data (NRT mode). Use --date or --start/--end "
                     "to backfill missing days."
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--date', metavar='YYYYMMDD',
                        help="Backfill a single specific day.")
    group.add_argument('--start', metavar='YYYYMMDD',
                        help="First day of a backfill range (inclusive).")
    parser.add_argument('--end', metavar='YYYYMMDD',
                         help="Last day of a backfill range (inclusive). "
                              "Only valid together with --start; defaults to --start.")
    args = parser.parse_args()

    if args.end and not args.start:
        parser.error("--end requires --start")
    if args.date and args.end:
        parser.error("--date cannot be combined with --end")

    return args


def get_time_range(args):
    """
    Returns (start, end) datetimes, where `end` is exclusive (the day after
    the last day to process).
    """
    yesterday = datetime.today().replace(hour=0, minute=0, second=0, microsecond=0) - timedelta(days=1)

    if args.date:
        start = datetime.strptime(args.date, '%Y%m%d')
        end_day = start
    elif args.start:
        start = datetime.strptime(args.start, '%Y%m%d')
        end_day = datetime.strptime(args.end, '%Y%m%d') if args.end else start
    else:
        start = yesterday
        end_day = yesterday

    return start, end_day + timedelta(days=1)


