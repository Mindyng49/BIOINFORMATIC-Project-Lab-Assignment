from Bio import Entrez
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
import streamlit as st

# Function 1: Retrieve Protein Data
def retrieve_data(id):
    Entrez.email = "mindyngyu@graduate.utm.my"

    handle = Entrez.efetch(db="protein",id=id,rettype="fasta",retmode="text")

    record = SeqIO.read(handle, "fasta")
    handle.close()

    return record

# Function 2: Get Basic Analysis
def get_basic_analysis(sequence):
    seq_analysis = ProteinAnalysis(sequence)

    result1 = len(sequence)
    result2 = seq_analysis.count_amino_acids()
    result3 = seq_analysis.molecular_weight()
    result4 = seq_analysis.isoelectric_point()

    return result1, result2, result3, result4

# Streamlit User Interface
st.set_page_config(
    layout="wide"
)

st.title('Lab 1 - Mindy Ng')

with st.form("protein_search"):
    protein_id = st.text_input('Enter Uniprot ID', placeholder="e.g. P04637")
    retrieve = st.form_submit_button('Retrieve')

if retrieve:

    if protein_id.strip() != "":

        try:
            with st.spinner("Retrieving protein data..."):

                record = retrieve_data(protein_id.strip())
                sequence = record.seq

                length, composition, weight, pi = get_basic_analysis(sequence)

            st.success("Protein data retrieved successfully!")

            # Two-column layout
            left, right = st.columns([1.2, 1])

            # LEFT: Retrieved Protein
            with left:
                st.subheader("📄 Retrieved Protein")

                st.write("**Name:**", record.name)
                st.write("**Description:**", record.description)

                st.write("**Sequence:**")

                sequence = str(record.seq)

                formatted_sequence = "\n".join(
                    sequence[i:i+60] for i in range(0, len(sequence), 60)
                )

                st.code(formatted_sequence)

            # RIGHT: Basic Protein Analysis
            with right:
                st.subheader("Basic Protein Analysis")

                st.write("Sequence Length:", length)

                st.write("Amino Acid Composition:")

                composition_text = ", ".join(
                [f"{amino_acid}: {frequency}"
                for amino_acid, frequency in composition.items()]
            )   

                st.write(composition_text)

                st.write("Molecular Weight:", weight)

                st.write("Isoelectric Point:", pi)


        except Exception as e:
            st.error(f"Unable to retrieve protein data: {e}")

    else:

        st.warning("Please enter a UniProt ID.")