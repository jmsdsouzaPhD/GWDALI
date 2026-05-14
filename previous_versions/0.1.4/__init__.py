# __init__ (1)
import sys
from pathlib import Path

path = Path(__file__).parent
sys.path.insert(1,path/'GWDALI')

from GWDALI.GWDALI import *
