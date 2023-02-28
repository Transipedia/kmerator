#### exit.py
import sys
import shutil
from color import *


def gracefully(args):
    ### quit properly
    if not args.keep:
        try:
            shutil.rmtree(args.tmpdir, ignore_errors=True)
        except:
            exit(f"Error deleting {args.tmpdir}.{ENDCOL}")
    sys.exit(0)
