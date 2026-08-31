import os
import shutil
import logging
import multiprocessing
import yaml

from src import preprocess_nc, SCLF1C_reg2fit
from src.helper.util import *

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

# ---------------------------------------------------------
# Configuration 
# ---------------------------------------------------------
CONFIG_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),'config.yaml')
with open(CONFIG_PATH) as f:
    _cfg = yaml.safe_load(f)['paths']

EE_TO_NC_SCRIPT = _cfg['ee_to_nc_script']
SMOS_BASE_PATH  = _cfg['smos_base_path']
TMP_NC_PATH     = _cfg['tmp_nc_path']
OUT_REG_PATH    = _cfg['out_reg_path']

def run_in_isolated_process(func, *args):
    """
    Runs a function in an isolated child process and waits for it to finish.
    This prevents memory leaks or state changes in external libraries from 
    affecting the main script.
    """
    p = multiprocessing.Process(target=func, args=args)
    p.start()
    p.join()

    if p.exitcode != 0:
        logging.error(f"Process running {func.__name__} failed with exit code {p.exitcode}")
        raise RuntimeError(f"{func.__name__} failed during execution.")

    logging.info(f"Successfully completed {func.__name__} via isolated process.")

def process_ee_to_nc(date_time: datetime) -> list:
    """
    Finds .zip ee files for the given date, converts them to .nc if needed,
    and returns a list of resulting netcdf files.
    """
    # Build search pattern for EE files
    year, month, day = date_time.strftime('Y%Y'), date_time.strftime('M%m'), date_time.strftime('%d')
    date_str = date_time.strftime('%Y%m%d')

    search_pattern = os.path.join(
        SMOS_BASE_PATH, year, month,
        f'SM_*_MIR_SCLF1C_{date_str}*_724_*zip'
    )

    eeflist = glob.glob(search_pattern)

    logging.info(f"[{date_str}] Found {len(eeflist)} zip files to process.")

    if not eeflist:
        raise RuntimeError(f"No SMOS EE zip files found for {date_str}")

    # Convert EE to NC
    for fee in eeflist:
        base_name = os.path.basename(fee)[:-3] + 'nc'
        target_nc_file = os.path.join(TMP_NC_PATH, base_name)

        if not os.path.isfile(target_nc_file):
            logging.info(f"Running ee_to_nc conversion on: {fee}")
            result = subprocess.run([
                EE_TO_NC_SCRIPT,
                '--target-directory', TMP_NC_PATH,
                fee
            ], capture_output=True, text=True)

            if result.returncode != 0:
                logging.error(f"Conversion script failed for {fee}:\n{result.stderr}")
                raise RuntimeError("ee_to_nc conversion corrupted")

            logging.info(f"Done ee_to_nc conversion: {fee}")

    # Return list of resulting NC files
    nc_search_pattern = os.path.join(TMP_NC_PATH, f'SM_*_MIR_SCLF1C_{date_str}*724*.nc')
    ncflist = glob.glob(nc_search_pattern)

    if len(eeflist) != len(ncflist):
        logging.warning(f"File count mismatch: {len(eeflist)} zip files vs {len(ncflist)} nc files.")

    return ncflist

def main():
    args = parse_args()
    start_time, end_time = get_time_range(args)
    logging.info(
        f"Processing {start_time:%Y-%m-%d} through "
        f"{(end_time - timedelta(days=1)):%Y-%m-%d}"
    )

    # Ensure necessary directories exist
    os.makedirs(os.path.join(TMP_NC_PATH, 'to_delete'), exist_ok=True)
    os.makedirs(OUT_REG_PATH, exist_ok=True)

    current_date = start_time

    while current_date < end_time:
        date_str = current_date.strftime('%Y%m%d')
        logging.info(f"=== Starting processing for {date_str} ===")

        # Step 1: Convert files to NetCDF
        try:
            ncflist = process_ee_to_nc(current_date)
        except RuntimeError as e:
            logging.error(f"Halting processing for {date_str} due to error: {e}")
            current_date += timedelta(days=1)
            continue

        # Step 2: Preprocess NetCDF files into REG 
        for fnc in ncflist:
            logging.info(f"Preprocessing NC file: {fnc}")
            try:
                run_in_isolated_process(preprocess_nc, fnc, OUT_REG_PATH)
            except RuntimeError as e:
                logging.error(f"Skipping {fnc} due to error: {e}")
                continue

            # Move processed file to 'to_delete'
            dest = os.path.join(TMP_NC_PATH, 'to_delete', os.path.basename(fnc))
            shutil.move(fnc, dest)
            logging.info(f"Moved {fnc} to cleanup directory.")

        # Step 3: Run REG to FIT processing for Ascending and Descending
        next_date = current_date + timedelta(days=1)

        logging.info("Running SCLF1C_reg2fit for Ascending (_A)")
        try:
            run_in_isolated_process(SCLF1C_reg2fit, OUT_REG_PATH, current_date, next_date, '_A')
        except RuntimeError as e:
            logging.error(f"SCLF1C_reg2fit failed for Ascending (_A) on {date_str}: {e}")

        logging.info("Running SCLF1C_reg2fit for Descending (_D)")
        try:
            run_in_isolated_process(SCLF1C_reg2fit, OUT_REG_PATH, current_date, next_date, '_D')
        except RuntimeError as e:
            logging.error(f"SCLF1C_reg2fit failed for Descending (_D) on {date_str}: {e}")
        
        logging.info(f"=== Completed processing to Tb40 for {date_str} ===")
        
        # Advance to the next day
        current_date += timedelta(days=1)


if __name__ == '__main__':
    # Set the multiprocessing start method exactly once here
    try:
        multiprocessing.set_start_method('spawn')
    except RuntimeError:
        pass # Context was already set
        
    main()
