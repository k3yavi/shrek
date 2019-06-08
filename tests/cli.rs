use std::path::PathBuf;
use std::process::Command;

use assert_cmd::prelude::*;
use predicates::prelude::*;

#[test]
fn shrek_generate() -> Result<(), Box<std::error::Error>> {
    let mut filename = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    filename.push("tests/data/fasta.fa");

    let mut cmd = Command::cargo_bin("shrek")?;
    cmd.arg("generate");
    cmd.args(&["--fasta", filename.to_str().expect("Demo data not found")]);
    cmd.args(&["--gfa", "a"]);

    cmd.assert().success();

    Ok(())
}

#[test]
fn shrek_compare() -> Result<(), Box<std::error::Error>> {
    let mut filename_1 = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    filename_1.push("tests/data/fasta1.gfa");

    let mut filename_2 = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    filename_2.push("tests/data/fasta2.gfa");

    let mut cmd = Command::cargo_bin("shrek")?;
    cmd.arg("compare");
    cmd.args(&["-1", filename_1.to_str().expect("Demo data not found")]);
    cmd.args(&["-2", filename_2.to_str().expect("Demo data not found")]);

    cmd.assert().success();

    Ok(())
}
