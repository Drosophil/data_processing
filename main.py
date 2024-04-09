from mpp import MolecularPropertiesProcessor

if __name__ == '__main__':
    mpp_instance = MolecularPropertiesProcessor(
        input_file_path="sample.csv",
        output_file_name="output.csv",
        smiles_column="clean_smiles",
        molecule_id_column="inchikey",
    )
    final_result = mpp_instance.process_data()
