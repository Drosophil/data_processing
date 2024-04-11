from multiprocessing import Process, Queue, cpu_count
from time import time, ctime
import pandas as pd
import re
from rdkit.Chem import AllChem, Descriptors
from mylog import logger  # loading logging configuration
import mpp.mpp
from mpp.mpp import MolecularPropertiesProcessor, MoleculeProcessingException

class MolecularPropertiesProcessorIO(MolecularPropertiesProcessor):
    '''
    Really huge file molecular data processing utility
    input parameters: input_file, output_file, smiles_column_name, molecule_id_column_name
    reads the file by chunks.
    For the very large files that do not fit in memory.
    NB: Output file overwrites any existing file with the same name.
    '''
    def __init__(
            self,
            input_file_path: str,
            output_file_name: str,
            smiles_column: str,
            molecule_id_column: str,
            chunk_size: int,
            hyperthreading=True,
    ):
        self.start_time = time()
        current_time = ctime()
        logger.info(f"<<<<<<< NEW IO CHUNKS RUN >>>>>>>  {current_time}")
        logger.info(f"DATA: <<<<<<< NEW IO CHUNKS RUN >>>>>>>  {current_time}")
        logger.info(f"PROCESS: <<<<<<< NEW IO CHUNKS RUN >>>>>>> {current_time}")
        self.hyperthreading = hyperthreading
        self.smiles_col_re = "^" + smiles_column + "$"
        self.mol_name_col_re = "^" + molecule_id_column + "$"
        self.input_file_name = input_file_path
        self.output_file_name = output_file_name
        self.max_cpu = cpu_count()
        self.chunk_size = chunk_size
        if hyperthreading:
            self.max_cpu //= 2  # we need real cores, not hyperthreading ones
        if self.max_cpu < 2:
            self.max_cpu = 2  #  for the algorythm to work: 1 core for reader, 1 for worker.
            #  it will work on 1 core, but very slow.

        self.input_chunks = Queue()   #  communication between processes
        self.result_chunks = Queue()


    @staticmethod
    def _column_finder(chunk_df: pd.DataFrame, match_str):
        '''finds a column in a dataframe'''
        matcher = re.compile(match_str, re.IGNORECASE)
        column_to_find = next(filter(matcher.match, chunk_df))
        if not column_to_find:
            logger.error(f"DATA: No {match_str} column found in a dataframe")
            raise MoleculeProcessingException(f"DATA: No {match_str} column found in a dataframe")
        return column_to_find

    def _prepare_data(self):
        '''overrided method that is useless here'''
        pass

    def _file_reader_and_writer(self):
        '''
        Reader process function.
        Reads input file, puts chunks in input queue for workers.
        If there are enough chunks in input queue, then
        gets results from output queue, writes them to an output file.
        Upon finishing the input file, issues STOP signs to workers and
        checks if they received them.
        '''
        logger.info("PROCESS: Reader|Writer started.")
        chunk_id = 0
        stop_counter = 0
        with pd.read_csv(self.input_file_name, chunksize=self.chunk_size) as reader:
            for chunk in reader:
                chunk_id += 1
                smiles_col = self._column_finder(chunk, self.smiles_col_re)
                mol_name_col = self._column_finder(chunk, self.mol_name_col_re)

                logger.info(f"Preparing data in chunk {chunk_id}")
                chunk = chunk[
                    [smiles_col, mol_name_col]
                    + list(chunk.columns.difference([smiles_col, mol_name_col]))
                    ]
                logger.info(f"dropping duplicates in chunk {chunk_id}")
                chunk.drop_duplicates(subset=mol_name_col, inplace=True)
                while stop_counter < (self.max_cpu - 1):
                    if self.input_chunks.qsize() < self.max_cpu:
                        self.input_chunks.put((1, chunk_id, chunk, smiles_col, mol_name_col))
                        break
                    else:
                        out_chunk = self.result_chunks.get()
                        if out_chunk[0]:
                            logger.info(f"PROCESS: writing chunk {out_chunk[1]}")
                            if out_chunk[1] == 1:  # if this is the first chunk, make new file with header
                                out_chunk[2].to_csv(self.output_file_name)  #  open new file, overwrite existing
                            else:
                                # no column names on append
                                out_chunk[2].to_csv(self.output_file_name, mode="a", header=False)

        logger.info("Finished reading the input file. Stopping worker processes...")
        logger.info(f"PROCESS: Issuing {self.max_cpu - 1} STOP signs...")

        for i in range(0, self.max_cpu - 1):
            self.input_chunks.put((0,0,0,0))  # stop signs for all workers
        # for this loop a timeout could be a nice option
        while (stop_counter < (self.max_cpu - 1)):  #  to ensure that all workers are stopped
            out_chunk = self.result_chunks.get()
            if out_chunk[0]:
                logger.info(f"PROCESS: writing chunk {out_chunk[1]}")
                if out_chunk[1] == 1:  # if this is the first chunk, make new file with header
                    out_chunk[2].to_csv(self.output_file_name)  # open new file, overwrite existing
                else:
                    # no column names on append
                    out_chunk[2].to_csv(self.output_file_name, mode="a", header=False)
            else:
                if stop_counter < (self.max_cpu - 1):
                    stop_counter += 1
                    logger.info(f"PROCESS: STOP sign {stop_counter} received.")
        logger.info(f"All STOP signs received: {stop_counter}={self.max_cpu - 1}.")

    def _compute_molecule_properties_chunk(self):
        '''
        Worker process function.
        Reads self.input_chunks queue, puts the result in self.result_chunks
        Upon receiving a STOP sign (0,0,0) puts it to output queue and stops.
        '''
        while True:
            start_time = time()
            chunk = self.input_chunks.get()
            if chunk[0]:
                logger.info(f"PROCESS: Started chunk {chunk[1]}.")

                def get_mol_from_smiles_and_verify(input_s):
                    '''generating mol object from smiles'''
                    output_mol = AllChem.MolFromSmiles(input_s)
                    if output_mol:
                        return output_mol
                    else:
                        logger.warning(f"DATA: chunk {chunk[1]} -> bad SMILES {input_s}")
                        return None

                chunk[2]["mol"] = chunk[2][chunk[3]].apply(lambda s: get_mol_from_smiles_and_verify(s))

                #  logger.info("PROCESS: Migrated to molecules")

                #  drop None should be done by only "mol" column

                chunk[2].dropna(subset=["mol"], inplace=True)  # drop rows with None values in Mol column


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

                chunk[2][mol_props_to_compute] = chunk[2].apply(
                    lambda row: [mol_props_funcs[prop](row["mol"]) for prop in mol_props_to_compute],
                    axis=1,
                    result_type="expand"
                )

                #  logger.info("PROCESS: Mol functions applied")

                chunk[2].drop(columns=["mol"], inplace=True)
                chunk[2].set_index(chunk[4], inplace=True)

                self.result_chunks.put((1,chunk[1],chunk[2]))  #  queuing the result

                stop_time = time() - start_time

                logger.info(f"PROCESS: Process finished chunk {chunk[1]} in {stop_time} seconds")
            else:
                break
        self.result_chunks.put((0, 0, 0))  # STOP sign received
        logger.info(f"PROCESS: worker received STOP sign. STOPPED.")

    def _compute_molecule_properties(self):
        """
        Compute molecule properties and fingerprints using RDKit
        for extremely large files, in chunks.
        Set up reader, worker processes.
        """
        logger.info("PROCESS: Setting up processes.")
        reader = Process(target=self._file_reader_and_writer)
        adders = [Process(target=self._compute_molecule_properties_chunk)
                  for i in range(0, self.max_cpu - 1)]
        logger.info("PROCESS: Start all.")
        start = time()
        reader.start()
        for proc in adders:
            proc.start()
        reader.join()
        for proc in adders:
            proc.join()
        logger.info(f"PROCESS: All closed. CPU time {time() - start} sec.")

    def process_data(self):
        '''
        main processing func to call externally
        '''

        self._compute_molecule_properties()
        stop_time = time()
        logger.info(f"PROCESS: wall execution time: {stop_time - self.start_time} sec.")

if __name__ == '__main__':
    logger.error("mppio module is __main__")

