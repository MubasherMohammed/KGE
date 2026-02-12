#!/usr/bin/env python
"""MAGeCK-KGE CLI entry point.

This is a drop-in replacement for the original mageck CLI that includes
the interactive HTML report subcommand.

Copyright (c) 2015 Wei Li, Han Xu, Xiaole Liu lab
Modified 2025 with interactive Plotly HTML reports and pathway enrichment.
"""

from __future__ import print_function
import sys
import logging

from mageck.crisprFunction import *
from mageck.mageckCount import *
from mageck.pathwayFunc import *
from mageck.argsParser import *
from mageck.testVisual import *
from mageck.version import __version__
from mageck.htmlReport import mageck_report_main


def main():
    args = crisprseq_parseargs()
    logging.info('Welcome to MAGeCK-KGE v' + __version__ + '. Command: ' + args.subcmd)

    if args.subcmd == 'run' or args.subcmd == 'count':
        mageckcount_main(args)

    if args.subcmd == 'run' or args.subcmd == 'test':
        magecktest_main(args)

    if args.subcmd == 'pathway':
        mageck_pathwaytest(args)

    if args.subcmd == 'plot':
        plot_main(args)

    if args.subcmd == 'report':
        mageck_report_main(args)


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me! ;-) Bye!\n")
        sys.exit(0)
