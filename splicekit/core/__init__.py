"""
# Description
Core modules and methods of splicekit.
"""

import os
import splicekit.core as core
import splicekit

current_version = 0.3

def version():
    print("splicekit v{version} ({folder})".format(folder=os.path.dirname(splicekit.__file__), version=core.current_version))
    print("Repository: https://code.roche.com/PMDA/cross-project/spliceosome-focus-team/splicekit")
    print()
