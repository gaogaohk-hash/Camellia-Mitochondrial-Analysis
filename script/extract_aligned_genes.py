
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation
import pandas as pd


blast_results = pd.read_csv("xiaye.txt", sep="\t", header=None)  # 替换为实际的文件名和路径


blast_results.columns = ["query_contig", "query_start", "query_end", "target_contig", "target_start", "target_end"]


blast_data = []
for _, row in blast_results.iterrows():
    blast_data.append({
        "contig": row["target_contig"],
        "start": int(row["target_start"]),
        "end": int(row["target_end"])
    })


genbank_file = "mtXIAYE.gb"  
genome = SeqIO.read(genbank_file, "genbank")


matches = []
for result in blast_data:
    for feature in genome.features:
        if feature.type == "gene":
            gene_location = feature.location
            gene_start = gene_location.start
            gene_end = gene_location.end
            gene_length = gene_end - gene_start

          
            overlap_start = max(result["start"], gene_start)
            overlap_end = min(result["end"], gene_end)
            overlap_length = max(0, overlap_end - overlap_start)


            if overlap_length >= 0.8 * gene_length:
                gene_name = feature.qualifiers.get("gene", ["unknown"])[0]  
                matches.append({
                    "contig": result["contig"],
                    "blast_start": result["start"],
                    "blast_end": result["end"],
                    "gene_name": gene_name,
                    "gene_start": gene_start,
                    "gene_end": gene_end,
                    "overlap_percentage": overlap_length / gene_length * 100  
                })


output_file = "mtXIAYEmatched.xlsx"
matched_genes_df = pd.DataFrame(matches)
matched_genes_df.to_excel(output_file, index=False)

print(f"Genes with ≥80% coverage have been saved to ‘{output_file}'")

