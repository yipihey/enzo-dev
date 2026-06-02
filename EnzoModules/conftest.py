import os
import sys

# Make the `enzomodules` package importable when running pytest from this
# directory without an editable install.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
