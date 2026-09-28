from pathlib import Path
import runpy
runpy.run_path(str(Path(__file__).parents[1]/"run_smallN_campaign.py"),run_name="__main__")
