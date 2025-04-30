import sys
from pathlib import Path

# Ensure src/ is on the Python path for absolute imports
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
