The project is mpp library.
Usage:
from mpp import MolecularPropertiesProcessor, MolecularPropertiesProcessorIO

1. MolecularPropertiesProcessor class.
This class loads the whole csv file in memory and transforms it.
Usage:
mpproc = MolecularPropertiesProcessor(
        input_file_path: str,
        output_file_name: str,
        smiles_column: str,
        molecule_id_column: str,
        hyperthreading: Bool  # we need real physical cores for the processes
    )
    final_result = mpproc.process_data()

process_data() -> pandas.DataFrame is a method to start processing.
input_file_path - path to the input csv file.
output_file_name - file to where the transformed data will be saved in a csv format.
smiles_column - a string with which the SMILES data are marked in input csv.
molecule_id_column - a column containing unique ID's.
hyperthreading - True is system CPU has hyperthreading, False if it has not. Important for speed adjustment,
it is better to operate on physical cores, not hyperthreaded abstract ones.
Int the output there should be no duplicates.
NB: Output file is opened in append mode, so if it exists, it will grow in size.

2. MolecularPropertiesProcessorIO class.
This class should be used if the size of an input file is unknown or it is too large to load in memory.
Class does not return pandas DataFrame, but saves transformed data in output csv file. Reads input file by chunks and 
implements simultaneous file I/O and transforming the data.
Drops duplicates only within chunks, so some duplicates can be in the output.
Usage:
    mpproc = MolecularPropertiesProcessorIO(
        input_file_path="sample.csv",
        output_file_name="output.csv",
        smiles_column="clean_smiles",
        molecule_id_column="inchikey",
        chunk_size=650,  # 650 lines perform faster on my laptop, it should be adjusted to the system
        hyperthreading=True  # we need real physical cores for the processes
    )
    mpproc.process_data()

process_data() method starts the transformation.
chunk_size is a chunk reading size in lines.
Other params are the same as for MolecularPropertiesProcessor class.
NB: Output file is opened in append mode, so if it exists, it will grow in size.

Both classes do logging, they produce 3 log files:
data_log.log  contains messages regarding data issues (like bad SMILES)
process_log.log contains messages about multiprocessing.
mainlog.log is a root logger output, root also reflects on console.
Please note that the results are not sorted, 
so these two classes may produce for the same input the same results, but in a different order.