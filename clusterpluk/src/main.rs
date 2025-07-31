use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use regex::Regex;
use log::{info, warn, error};
use bio::io::fastq;
use std::env;

fn main() {
    env_logger::init();

    let args: Vec<String> = env::args().collect();
    if args.len() != 6 {
        eprintln!("Usage: {} <R1.fastq> <R2.fastq> <clusters.clstr> <output_R1.fastq> <output_R2.fastq>", args[0]);
        std::process::exit(1);
    }

    let r1_file = &args[1];
    let r2_file = &args[2];
    let cluster_file = &args[3];
    let r1_out = &args[4];
    let r2_out = &args[5];

    sort_fastq_by_quality(r1_file, r2_file, cluster_file, r1_out, r2_out);
}

fn sort_fastq_by_quality(
    r1_input: &str,
    r2_input: &str,
    cluster_file: &str,
    r1_output: &str,
    r2_output: &str,
) {
    let r1_reader = fastq::Reader::from_file(r1_input).expect("Failed to open R1 FASTQ");
    let r2_reader = fastq::Reader::from_file(r2_input).expect("Failed to open R2 FASTQ");

    info!("\nCreating index for: {}", r1_input);
    let mut r1_index = HashMap::new();
    for result in r1_reader.records() {
        let record = result.expect("Error reading R1 record");
        r1_index.insert(record.id().to_string(), record.clone());
    }

    info!("\nCreating index for: {}", r2_input);
    let mut r2_index = HashMap::new();
    for result in r2_reader.records() {
        let record = result.expect("Error reading R2 record");
        r2_index.insert(record.id().to_string(), record.clone());
    }

    let re = Regex::new(r".*>(.*)\.\.\.").unwrap();
    let reader = BufReader::new(File::open(cluster_file).expect("Could not open cluster file"));
    let mut out1 = fastq::Writer::to_file(r1_output).expect("Failed to create output R1 file");
    let mut out2 = fastq::Writer::to_file(r2_output).expect("Failed to create output R2 file");
    // let mut out1 = File::create(r1_output).expect("Failed to create output R1 file");
    // let mut out2 = File::create(r2_output).expect("Failed to create output R2 file");

    let mut cluster: Vec<(fastq::Record, fastq::Record)> = Vec::new();
    let mut cluster_count = 0;

    for (i, line) in reader.lines().enumerate() {
        let line = line.expect("Failed to read line");
        if i == 0 {
            continue; // skip header
        }

        if line.starts_with('>') {
            cluster_count += 1;
            print!("\rProcessing cluster: {}", cluster_count);
            std::io::stdout().flush().unwrap();

            let (r1_best, r2_best) = if cluster.len() == 1 {
                cluster[0].clone()
            } else {
                pluck_best_read_from_cluster(&cluster)
            };
            out1.write_record(&r1_best).expect("Failed to write to R1 file");
            out2.write_record(&r2_best).expect("Failed to write to R2 file");

            cluster.clear();
        } else {
            if let Some(cap) = re.captures(&line) {
                let id = &cap[1];
                if let (Some(r1), Some(r2)) = (r1_index.get(id), r2_index.get(id)) {
                    cluster.push((r1.clone(), r2.clone()));
                } else {
                    error!("Read ID {} not found in FASTQ files", id);
                }
            } else {
                warn!("Malformed line: {}", line);
            }
        }
    }

    println!("\nProcessing complete. Clusters processed: {}", cluster_count);
}

fn pluck_best_read_from_cluster(cluster: &[(fastq::Record, fastq::Record)]) -> (fastq::Record, fastq::Record) {
    let mut scores = Vec::new();
    let mut seqs = Vec::new();

    for (r1, r2) in cluster {
        let full_seq: Vec<u8> = r1.seq().iter().chain(r2.seq()).copied().collect();
        let full_qual: Vec<u8> = r1.qual().iter().chain(r2.qual()).copied().collect();

        let avg_error: f64 = full_qual
            .iter()
            .map(|&q| 10f64.powf(-(q as f64 - 33.0) / 10.0))
            .sum::<f64>() / full_qual.len() as f64;

        scores.push(avg_error);
        seqs.push(full_seq);
    }

    let mut idxs: Vec<usize> = (0..scores.len()).collect();
    idxs.sort_by(|&i, &j| scores[i].partial_cmp(&scores[j]).unwrap());

    let mut counts = HashMap::new();
    for seq in &seqs {
        *counts.entry(seq.clone()).or_insert(0) += 1;
    }

    let consensus_seq = counts.into_iter().max_by_key(|&(_, c)| c).unwrap().0;
    let best_index = idxs.into_iter().find(|&i| seqs[i] == consensus_seq).unwrap();

    cluster[best_index].clone()
}