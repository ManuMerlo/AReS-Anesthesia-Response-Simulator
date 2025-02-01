# Third party library imports
import pandas as pd

# Standard library imports
from pathlib import Path


def read_data(file_name, base_dir, row_index=None, delimiter=';'):
    """
    Reads data from a file located within a directory structure.

    Parameters:
    - file_name (str): The name of the file.
    - base_dir (str or Path): The base directory (can include subdirectories).
    - row_index (int, optional): Row index to read from.

    Returns:
    - Data from the file (depends on file type and implementation).
    """

    # Convert base_dir to a Path object if it's a string
    base_dir = Path(base_dir)

    if not base_dir.is_absolute():
        base_dir = base_dir.resolve()  # Converts to absolute path

    # Construct the full file path
    file_path = base_dir / file_name

    # Check if the file exists
    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")

    df = pd.read_csv(file_path, delimiter=delimiter)

    if row_index is not None:
        row_data = df.iloc[row_index].tolist()
        return row_data
    else:
        data_list = df.values.tolist()
        return data_list
