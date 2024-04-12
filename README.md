The project is mpp library.<br>
Usage:<br>
<code>from mpp import MolecularPropertiesProcessor, MolecularPropertiesProcessorIO</code><br>
<br>
<h3>1. MolecularPropertiesProcessor class.</h3>
This class loads the whole csv file in memory and transforms it.<br>
Usage:<br><br>
<code>mpproc = MolecularPropertiesProcessor(
        input_file_path: str,
        output_file_name: str,
        smiles_column: str,
        molecule_id_column: str,
        hyperthreading: Bool  # we need real physical cores for the processes
)
final_result = mpproc.process_data()</code><br><br>
process_data() -> pandas.DataFrame is a method to start processing.<br>
input_file_path - path to the input csv file.<br>
output_file_name - file to where the transformed data will be saved in a csv format.<br>
smiles_column - a string with which the SMILES data are marked in input csv.<br>
molecule_id_column - a column containing unique ID's.<br>
hyperthreading - <b>True</b> if system CPU does have hyperthreading, <b>False</b> if it doesn't. Important for speed adjustment, it is better to operate on physical cores, not hyperthreaded abstract ones.<br>
Int the output there should be no duplicates.<br>
NB: Output file will be overwritten if exists.<br><br>

<h3>2. MolecularPropertiesProcessorIO class.</h3>
This class should be used if the size of an input file is unknown or it is too large to load in memory.<br>
Class does not return pandas DataFrame, but saves transformed data in output csv file. Reads input file by chunks and 
implements simultaneous file I/O and transforming the data.<br>
Drops duplicates only within chunks, so some duplicates can be in the output.<br><br>
Usage:<br><br>
<code>mpproc = MolecularPropertiesProcessorIO(<br>
        input_file_path="sample.csv",
        output_file_name="output.csv",
        smiles_column="clean_smiles",
        molecule_id_column="inchikey",
        chunk_size=650,  # 650 lines perform faster on my laptop, it should be adjusted to the system
        hyperthreading=True  # we need real physical cores for the processes
)
mpproc.process_data()</code><br><br>
process_data() method starts the transformation.<br>
chunk_size is a chunk reading size in lines.<br>
Other params are the same as for MolecularPropertiesProcessor class.<br>
NB: Output file will be overwritten if exists.<br><br>

Both classes do logging, they produce 3 log files:<br>
<ul><li>data_log.log  contains messages regarding data issues (like bad SMILES)</li>
<li>process_log.log contains messages about multiprocessing.</li>
<li>mainlog.log is a root logger output, root also reflects on console.</li></ul><br>
Please note that the results are not sorted, so these two classes may produce the same results for the same input, but row order will differ.<br>
