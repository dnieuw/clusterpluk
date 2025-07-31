import argparse
from Bio import SeqIO
import numpy as np
import re
import logging

# Setup logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def sort_fastq_by_quality(R1_input_fastq, R2_input_fastq, cluster_file, R1_output_fastq, R2_output_fastq):
    logging.info(f"Indexing R1 FASTQ: {R1_input_fastq}")
    R1_dict = SeqIO.index(R1_input_fastq, "fastq")
    logging.info(f"Indexing R2 FASTQ: {R2_input_fastq}")
    R2_dict = SeqIO.index(R2_input_fastq, "fastq")

    regex = r".*>(.*)\.\.\."
    total_clusters = 0

    logging.info(f"Starting processing: {cluster_file}")
    with open(cluster_file, "r") as clstr, open(R1_output_fastq, "w") as out1, open(R2_output_fastq, "w") as out2:
        cluster = []
        for n, line in enumerate(clstr):
            line = line.strip()
            if n == 0:
                continue
            if line.startswith(">"):
                print(f"Processing cluster: {total_clusters}", end='\r')
                total_clusters += 1
                if len(cluster) == 1:
                    best_R1, best_R2 = cluster[0]
                else:
                    best_R1, best_R2 = pluck_best_read_from_cluster(cluster)
                SeqIO.write(best_R1, out1, "fastq")
                SeqIO.write(best_R2, out2, "fastq")
                cluster = []
            else:
                match = re.match(regex, line)
                if not match:
                    logging.warning(f"Line {n} didn't match expected format: {line}")
                    continue
                read_id = match.group(1)
                try:
                    cluster.append((R1_dict[read_id], R2_dict[read_id]))
                except KeyError as e:
                    logging.error(f"Read ID not found in FASTQ files: {read_id}")

    logging.info(f"Processing complete: {total_clusters} clusters processed.")

def pluck_best_read_from_cluster(cluster):
    #logging.debug("Extracting best read from cluster of size: %d", len(cluster))
    nt_mat = np.array([read[0].seq + read[1].seq for read in cluster])
    q_mat = np.array([read[0].letter_annotations["phred_quality"] + read[1].letter_annotations["phred_quality"] for read in cluster])
    average_error = np.mean(np.power(10, -q_mat / 10), axis=1)
    error_index = np.argsort(average_error)

    uniq, index, counts = np.unique(nt_mat, return_counts=True, return_inverse=True, axis=0)
    consensus_index = error_index[np.where(index == np.argmax(counts))]
    best_read = cluster[np.max(consensus_index)]

    #logging.debug("Best read selected with average error: %.4f", average_error[np.max(consensus_index)])
    return best_read[0], best_read[1]

def main():
    parser = argparse.ArgumentParser(description="Pluck the best read from each read cluster based on consensus and phred quality")
    parser.add_argument('-i', '--input', required=True, help='Input R1 FASTQ file')
    parser.add_argument('-i2', '--input2', required=True, help='Input R2 FASTQ file')
    parser.add_argument('-c', '--cluster', required=True, help='Input .clstr file')
    parser.add_argument('-o', '--output', required=True, help='Output R1 sorted FASTQ file')
    parser.add_argument('-o2', '--output2', required=True, help='Output R2 sorted FASTQ file')

    args = parser.parse_args()
    logging.info("Starting processing")
    sort_fastq_by_quality(args.input, args.input2, args.cluster, args.output, args.output2)
    logging.info("Finished processing")

if __name__ == "__main__":
    main()

# Alternative approach:
# def sort_fastq_by_quality2(R1_input_fastq, R2_input_fastq, cluster_file, R1_output_fastq, R2_output_fastq):
#     regex = r".*>(.*)\.\.\."
#     with open(cluster_file,"r") as clstr:
#         read_dict = {}
#         cluster_list = []
#         n = 0
#         for line in clstr:
#             line = line.strip()
#             if line[0] == ">":
#                 if n > 0:
#                     cluster_list.append((n,[]))
#                 cluster_number = int(line.replace(">Cluster ",""))
#                 n = 0
#             else:
#                 read_dict[re.match(regex,line).group(1)] = cluster_number
#                 n += 1
#     print("passed cluster reading phase")
#     with open(R1_output_fastq, "w") as out1, open(R2_output_fastq, "w") as out2:
#         for R1,R2 in zip(SeqIO.parse(R1_input_fastq, "fastq"), SeqIO.parse(R2_input_fastq, "fastq")):
#             cluster_number = read_dict[R1.id]
#             cluster_list[cluster_number][1].append((R1, R2))
#             if len(cluster_list[cluster_number][1]) == cluster_list[cluster_number][0]:
#                 read1, read2 = pluck_best_read_from_cluster2(cluster_list[cluster_number][1])
#                 _ = SeqIO.write(read1, out1, "fastq")
#                 _ = SeqIO.write(read2, out2, "fastq")
#                 cluster_list[cluster_number] = []
#
# def pluck_best_read_from_cluster2(cluster):
#     #Create numpy arrays from the concatenated R1 and R2 read
#     nt_mat = np.array([read[0].seq + read[1].seq for read in cluster])
#     q_mat = np.array([read[0].letter_annotations["phred_quality"] + read[1].letter_annotations["phred_quality"] for read in cluster])
#     #Calculate the average error and get the order
#     average_error = np.mean(np.power(10,-q_mat/10), axis=1)
#     error_index = np.argsort(average_error)
#     #Count the unique read sequences in this cluster
#     uniq, index, counts = np.unique(nt_mat, return_counts=True, return_inverse=True, axis = 0)
#     #Get those indeces in the error index that are part of the major read
#     consensus_index = error_index[np.where(index==np.argmax(counts))]
#     #The read with the lowest error, thus highest sorted error index, that also has the most abundant sequence, is the best read
#     best_read = cluster[np.max(consensus_index)]
#     return(best_read[0],best_read[1])