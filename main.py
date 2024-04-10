from mpp import MolecularPropertiesProcessor, MolecularPropertiesProcessorIO

if __name__ == '__main__':
    # parallel model from the lectures
    mpp_instance = MolecularPropertiesProcessor(
        input_file_path="sample.csv",
        output_file_name="output.csv",
        smiles_column="clean_smiles",
        molecule_id_column="inchikey",
        hyperthreading=True  # we need real physical cores for the processes
    )
    final_result = mpp_instance.process_data()

    # reader-workers parallel model for insanely huge files
    mpp_instance = MolecularPropertiesProcessorIO(
        input_file_path="sample.csv",
        output_file_name="output.csv",
        smiles_column="clean_smiles",
        molecule_id_column="inchikey",
        chunk_size=650,  # 650 lines perform faster on my laptop, it should be adjusted to the system
        hyperthreading=True  # we need real physical cores for the processes
    )
    mpp_instance.process_data()
    # Result is in the file. Since input file supposed to be insanely huge, we cannot concat dataframe in memory.
