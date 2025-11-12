///module responsible for writing the output to a csv file

use csv::Writer;
use std::error::Error;

pub fn write_output (final_count: &Vec<u64>) -> Result<(), Box<dyn Error>> {
    let mut wtr = Writer::from_path("output.csv")?;

    let indices: Vec<u64> = (0..final_count.len() as u64).collect();
    wtr.write_record(indices.iter().map(|i| i.to_string()))?;

    wtr.write_record(final_count.iter().map(|v| v.to_string()))?;

    wtr.flush()?;
    Ok(())
}
