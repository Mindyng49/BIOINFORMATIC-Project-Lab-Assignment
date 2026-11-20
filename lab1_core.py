"""
lab1_core.py
============
Lab 1 logic lifted out of ``lab1-Mindy.py`` so that the Lab 2 application can
reuse it.  Nothing here is new behaviour: the two functions are exactly the
ones written for Lab 1 (retrieve a UniProt record, basic protein analysis).
"""

from Bio import Entrez
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

EMAIL = "mindyngyu@graduate.utm.my"


def retrieve_data(protein_id: str):
    """Fetch one protein record from NCBI/UniProt by accession."""
    Entrez.email = EMAIL
    handle = Entrez.efetch(db="protein", id=protein_id, rettype="fasta", retmode="text")
    try:
        record = SeqIO.read(handle, "fasta")
    finally:
        handle.close()
    return record


def get_basic_analysis(sequence):
    """Length, amino-acid composition, molecular weight and isoelectric point."""
    seq_analysis = ProteinAnalysis(str(sequence))
    return (
        len(sequence),
        seq_analysis.count_amino_acids(),
        seq_analysis.molecular_weight(),
        seq_analysis.isoelectric_point(),
    )
