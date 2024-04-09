import re
import multiprocessing
from time import time, ctime

import pandas as pd
import numpy as np
from rdkit.Chem import AllChem, Descriptors
from mylog import logger  # loading logging configuration
import math


class MoleculeProcessingException(Exception):
    '''Exception class'''
    pass

class MolecularPropertiesProcessor:
    '''Molecular data processing utility
    input parameters: input_file, output_file, smiles_column_name, molecule_id_column_name
    '''
    def __init__(
            self,
            input_file_path: str,
            output_file_name: str,
            smiles_column: str,
            molecule_id_column: str,
    ):
        self.start_time = time()
        current_time = ctime()
        logger.info(f"<<<<<<< NEW RUN >>>>>>>  {current_time}")
        logger.info(f"DATA: <<<<<<< NEW RUN >>>>>>>  {current_time}")
        logger.info(f"PROCESS: <<<<<<< NEW RUN >>>>>>> {current_time}")
        self.mols_df = pd.read_csv(input_file_path)

        self.output_file_name = output_file_name

        self.smiles_col = self._column_finder("^" + smiles_column + "$")
        self.mol_name_col = self._column_finder("^" + molecule_id_column + "$")
        logger.info(f"File read took {time() - self.start_time} seconds")

    def _column_finder(self, match_str):
        matcher = re.compile(match_str, re.IGNORECASE)
        column_to_find = next(filter(matcher.match, self.mols_df.columns))
        if not column_to_find:
            logger.error(f"DATA: No {match_str} column found in a dataframe")
            raise MoleculeProcessingException(f"DATA: No {match_str} column found in a dataframe")
        return column_to_find

    def _prepare_data(self):
        logger.info("Preparing data")
        self.mols_df = self.mols_df[
            [self.smiles_col, self.mol_name_col]
            + list(self.mols_df.columns.difference([self.smiles_col, self.mol_name_col]))
            ]
        logger.info("dropping duplicates")
        self.mols_df.drop_duplicates(subset=self.mol_name_col, inplace=True)

    def _compute_molecule_properties_chunk(
            self,
            chunk_df: pd.DataFrame,
            chunk_id: int,
    ) -> pd.DataFrame:
        """ Compute molecule properties for chunk dataframe """

        logger.info(f"PROCESS: process for chunk {chunk_id} started")
        #  logger.info("PROCESS: started migration")

        start_time = time()

        def get_mol_from_smiles_and_verify(input_s):
            '''generating mol object from smiles'''
            output_mol = AllChem.MolFromSmiles(input_s)
            if output_mol:
                return output_mol
            else:
                logger.warning(f"DATA: chunk {chunk_id} -> bad SMILES {input_s}")
                return output_mol

        chunk_df["mol"] = chunk_df[self.smiles_col].apply(lambda s: get_mol_from_smiles_and_verify(s))

        #  logger.info("PROCESS: Migrated to molecules")

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

        #  logger.info("PROCESS: applying mol lambda functions")

        chunk_df[mol_props_to_compute] = chunk_df.apply(
            lambda row: [mol_props_funcs[prop](row["mol"]) for prop in mol_props_to_compute],
            axis=1,
            result_type="expand"
        )

        #  logger.info("PROCESS: Mol functions applied")

        chunk_df.drop(columns=["mol"], inplace=True)
        chunk_df.set_index(self.mol_name_col, inplace=True)
        stop_time = time() - start_time

        logger.info(f"PROCESS: Process finished chunk {chunk_id} in {stop_time} seconds")

        return chunk_df

    def _compute_molecule_properties(self) -> pd.DataFrame:
        """
        Compute molecule properties and fingerprints using RDKit
        in chunks
        """
        start_time = time()
        logger.info("Entering _compute_molecule_properties()")
        # const_size_of_chunks = 5
        max_amount_of_p = multiprocessing.cpu_count() // 2  # to avoid hyperthreading
        logger.info(f"Max CPUs = {max_amount_of_p}")

        amount_of_chunk_df = max_amount_of_p  # one core left for the main process, it seems faster like that
        # math.ceil(len(self.mols_df) / const_size_of_chunks)

        if amount_of_chunk_df > max_amount_of_p:
            amount_of_chunk_df = max_amount_of_p
        elif amount_of_chunk_df == 0:
            amount_of_chunk_df = 1

        logger.info("numpy array split")

        list_of_chunks = np.array_split(self.mols_df, amount_of_chunk_df)

        logger.info("setting the pool")

        with multiprocessing.Pool(processes=amount_of_chunk_df) as pool:
            p_df = pool.starmap(self._compute_molecule_properties_chunk,
                            [(list_of_chunks[number-1], number) for number in range(1, amount_of_chunk_df + 1)])

        logger.info(f"Pool has finished, result type: {type(p_df)}")

        result = pd.concat(p_df)
        logger.info(f"Main compute function worked {time() - start_time} seconds")
        return result

    def process_data(self) -> pd.DataFrame:
        '''
        main processing func
        '''
        self._prepare_data()

        mol_properties_df = self._compute_molecule_properties()
        mol_properties_df.to_csv(self.output_file_name)
        stop_time = time()
        logger.info(f"PROCESS: wall execution time: {stop_time - self.start_time} sec.")

        return mol_properties_df

if __name__ == '__main__':
    logger.error("mpp module is __main__")

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
