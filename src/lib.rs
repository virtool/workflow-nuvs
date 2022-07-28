use pyo3::prelude::*;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter};
use std::process::Command;
use std::path::Path;
use std::collections::HashSet;
use std::io::Write;

#[pyclass]
#[derive(Clone)]
struct Reads {
    #[pyo3(get, set)]
    pub paired: bool,
    #[pyo3(get, set)]
    pub left: String,
    #[pyo3(get, set)]
    pub right: String,
    #[pyo3(get, set)]
    pub unmapped_hosts: String,
}

#[pymethods]
impl Reads {
    #[new]
    fn new(
        paired: bool, 
        left: String, 
        right: String, 
        unmapped_hosts: String
    ) -> Self  {
        Self { paired, left, right, unmapped_hosts }
    }
}

#[derive(Debug)]
struct FastqRecord{
    header: String,
    seq: String,
    quality: String,
}


impl FastqRecord {
    fn to_fastq(&self) -> String {
        format!("{}\n{}\n{}\n{}\n", &self.header, &self.seq, "+", &self.quality)
    }
}

fn read_fastq(reader: BufReader<File>) -> Vec<FastqRecord>{
    let mut records = Vec::new();

    let mut has_plus = false;
    let mut header = "".to_string();
    let mut seq = "".to_string();

    for line in reader.lines() {
        let line = line.unwrap();
        let head = line.chars().take(1).last().unwrap();

        if head == '+' {
            has_plus = true;
            continue;
        } 

        
        if !has_plus {

            if head == '@' {
                header = line.trim().to_string();
                continue;
            }

            seq = line.trim().to_string();
            continue;
        } 

        if has_plus {
            records.push(FastqRecord{ 
                header: header.to_string(),
                seq: seq.to_string(), 
                quality: line.trim().to_string(),
            });
            has_plus = false;
        }
    }
    records
}

fn read_fastq_headers(filename: &str) -> HashSet<String> {
    let file = File::open(filename).expect("File not Found");
    let reader = BufReader::new(file);

    let mut headers = HashSet::new();
    let mut has_plus = false;

    for line in reader.lines() {
        let line = line.unwrap();

        let head = line.chars().take(1).last().unwrap();

        if head == '+' {
            has_plus = true;
            continue;
        } 

        
        if head == '@' {
            let header_start = line.trim().split(' ').collect::<Vec<&str>>()[0];
            headers.insert(header_start.to_string());
        }

        if has_plus {
            has_plus = false;
        }
    }

    headers
}


fn unzip_gz(filename: &str, target: &str) -> BufReader<File> {
    let target_file = File::create(target).expect("Could not create unzipped file");

    let mut command = Command::new("pigz");

    command
        .arg("-p")
        .arg("2")
        .arg("-dc")
        .arg("-k")
        .arg(filename)
        .stdout(target_file);
        
    command
        .output()
        .expect("Failed to decompress with pigz");


    println!("{:#?}", &command);


    let read_file = File::open(target).expect("File not found");
    BufReader::new(read_file)
}


#[pyfunction]
fn reunite_pairs(reads: Reads, work_path: String) -> PyResult<Py<Reads>> {
    let unmapped_roots: HashSet<String> = read_fastq_headers(&reads.unmapped_hosts);

    let b = Path::new(&reads.left).exists();
    let b2 = Path::new(&reads.right).exists();
    println!("{}, {}", b, b2);

    let reads_left = unzip_gz(&reads.left, &reads.left.replace(".gz", ""));
    let reads_right = unzip_gz(&reads.right, &reads.right.replace(".gz", ""));


    let mut out_1 = BufWriter::new(File::create(Path::new(&work_path).join("unmapped_1.fq")).unwrap());
    let mut out_2 = BufWriter::new(File::create(Path::new(&work_path).join("unmapped_2.fq")).unwrap());

    let records_left = read_fastq(reads_left);
    let records_right = read_fastq(reads_right);

    for entry in records_left {
        let header_split: &Vec<&str> = &entry.header.split(' ').collect();
        let header_start = header_split[0];
        if unmapped_roots.contains(header_start) {
            write!(out_1, "{}", entry.to_fastq()).unwrap();
        }
    }

    for entry in records_right {
        let header_split: &Vec<&str> = &entry.header.split(' ').collect();
        let header_start = header_split[0];
        if unmapped_roots.contains(header_start) {
            write!(out_2, "{}", entry.to_fastq()).unwrap();
        }
    }

    let gil = Python::acquire_gil();
    let py = gil.python();
    Ok(Py::new(py, reads).unwrap())
}

#[pymodule]
fn rust(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(reunite_pairs, m)?)?;
    m.add_class::<Reads>()?;
    Ok(())
}
