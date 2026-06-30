import sys
import os

sys.path.append(os.path.abspath("./out/build/x64-Release"))

import spatial_engine

print(spatial_engine.hello_world())
