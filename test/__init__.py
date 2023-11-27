# add src folder of application to syspath
import sys
from os import path
sys.path.append(path.join(path.dirname(path.dirname(path.abspath(__file__))), "variant_classification"))