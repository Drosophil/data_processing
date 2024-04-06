import re
import multiprocessing
from time import time

import pandas as pd
import numpy as np
from rdkit.Chem import AllChem, Descriptors
import math


class MoleculeProcessingException(Exception):
    '''Exception class'''
    pass


class MolecularPropertiesProcessor:
    '''Molecular data processing utility'''
    def __init__(
            self,
            input_file_path: str,
            output_file_name: str,
    ):
        start_time = time()
        self.mols_df = pd.read_csv(input_file_path)
        print(self.mols_df)
        self.output_file_name = output_file_name

        self.smiles_col = self._column_finder("^clean_smiles$")
        self.mol_name_col = self._column_finder("^inchikey$")
        print(f"File read took {time() - start_time} seconds")

    def _column_finder(self, match_str):
        matcher = re.compile(match_str, re.IGNORECASE)
        column_to_find = next(filter(matcher.match, self.mols_df.columns))
        if not column_to_find:
            raise MoleculeProcessingException(f"No {match_str} column found in a dataframe")
        return column_to_find

    def _prepare_data(self):
        print("Preparing data")
        self.mols_df = self.mols_df[
            [self.smiles_col, self.mol_name_col]
            + list(self.mols_df.columns.difference([self.smiles_col, self.mol_name_col]))
            ]
        print("dropping duplicates")
        self.mols_df.drop_duplicates(subset=self.mol_name_col, inplace=True)

    def _compute_molecule_properties_chunk(
            self,
            chunk_df: pd.DataFrame,
    ) -> pd.DataFrame:
        """ Compute molecule properties for chunk dataframe """

        print("Entering _compute_molecule_properties_chunk")
        print("started migration")

        start_time = time()

        # TODO: alter the lambda func and log corrupted rows

        chunk_df["mol"] = chunk_df[self.smiles_col].apply(lambda s: AllChem.MolFromSmiles(s))
        print("Migrated to molecules")
        chunk_df.dropna(inplace=True)  # drop rows with None values in Mol column

        mol_props_funcs = {
            "Molecular weight": lambda mol: Descriptors.MolWt(mol),
            "TPSA": lambda mol: Descriptors.TPSA(mol),
            "logP": lambda mol: Descriptors.MolLogP(mol),
            "H Acceptors": lambda mol: Descriptors.NumHAcceptors(mol),
            "H Donors": lambda mol: Descriptors.NumHDonors(mol),
            "Ring Count": lambda mol: Descriptors.RingCount(mol),
            "Lipinski pass": lambda mol: all([
                Descriptors.MolWt(mol) < 500,
                Descriptors.MolLogP(mol) < 5,
                Descriptors.NumHDonors(mol) < 5,
                Descriptors.NumHAcceptors(mol) < 10
            ])
        }
        mol_props_to_compute = list(mol_props_funcs.keys())

        print("applying mol lambda functions")

        chunk_df[mol_props_to_compute] = chunk_df.apply(
            lambda row: [mol_props_funcs[prop](row["mol"]) for prop in mol_props_to_compute],
            axis=1,
            result_type="expand"
        )

        print("Mol functions applied")

        chunk_df.drop(columns=["mol"], inplace=True)
        chunk_df.set_index(self.mol_name_col, inplace=True)

        print(f"Process finished in {time() - start_time} seconds")

        return chunk_df

    def _compute_molecule_properties(self) -> pd.DataFrame:
        """
        Compute molecule properties and fingerprints using RDKit
        in chunks
        """
        start_time = time()
        print("Entering _compute_molecule_properties()")
        # const_size_of_chunks = 5
        max_amount_of_p = multiprocessing.cpu_count() // 2  # to avoid hyperthreading
        print(f"Max CPUs = {max_amount_of_p}")

        amount_of_chunk_df = max_amount_of_p  # math.ceil(len(self.mols_df) / const_size_of_chunks)

        if amount_of_chunk_df > max_amount_of_p:
            amount_of_chunk_df = max_amount_of_p
        elif amount_of_chunk_df == 0:
            amount_of_chunk_df = 1

        print("numpy stuff")

        list_of_chunks = np.array_split(self.mols_df, amount_of_chunk_df)

        print("setting the pool")

        with multiprocessing.Pool(processes=amount_of_chunk_df) as pool:
            p_df = pool.map(self._compute_molecule_properties_chunk, list_of_chunks)

        print(f"Pool has finished, result type: {type(p_df)}")

        # list_of_p = [p for p in p_df]
        result = pd.concat(p_df)
        print(f"Main compute function worked {time() - start_time} seconds")
        return result

    def process_data(self) -> pd.DataFrame:
        '''
        main processing func
        '''
        self._prepare_data()

        mol_properties_df = self._compute_molecule_properties()
        mol_properties_df.to_csv(self.output_file_name)

        return mol_properties_df

if __name__ == '__main__':
    mpp = MolecularPropertiesProcessor(
        input_file_path="sample.csv",
        output_file_name="output.csv",
    )

    final_result = mpp.process_data()
#    final_result.to_csv(self.output_file_name)

# TODO: we want split before reading the file

# threads - share resourses
#  why - to stop the process if the error occured
#   - to share the memory

# slowest part? 1. reading the file, 2. computing the columns/manipulating with data
# Solved: slowest is 2.

# 1. reading the file:
#  - no columns that we need


# 2. computing the columns/manipulating with data
#  - would take a lot of time to process the file

# process - don't share the resourses
#  why - parallel computing
#  - won't crush the program
#

#  think about the amount of chunks
#  handling errors
#  think about the solution for the IO file
#    - set_index how we can use it to get rid of the duplicates
#    - think about how to implement the threading
#  do we need to use processes if the file is small
#  the bad data?

# TTD 1. we write the tests 2. writing the code
