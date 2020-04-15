
extern crate docopt;
extern crate bio;
extern crate rust_htslib;

#[macro_use] extern crate lazy_static;
//extern crate istring;

use std::env;
use std::process::exit;
use std::io::{Write, stderr};
use docopt::{Docopt, ArgvMap};

#[macro_use] mod common;
mod call; 
mod call2;
mod genome_reader;
mod reader_of_bams;
mod perf_profiler;

const USAGE: &str = "
Mutato is a software toolkit for variant calling.

Usage:
  variant <subcommand>

Available subcommands:
  call      Call variants based on simple thresholds (uses samtools mpileup)
  call2     Call variants based on simple thresholds (no dependencies)
  somatic   Identify somatic mutations by comparing against normal samples
";

fn main() {
	// TODO: Use match args.as_slice() { [_, "detect"] => ... } after slice
	// pattern matching stabilizes (see issue #23121).
	let args: Vec<String> = env::args().collect();

	if args.len() >= 2 && args[1] == "call" { call::main(); }
	//else if args.len() >= 2 && args[1] == "somatic" { somatic::main(); }
    else if args.len() >= 2 && args[1] == "call2" { call2::main(); }
	else { println!("{}", USAGE); exit(-1); }
}

pub fn parse_args(usage: &str) -> ArgvMap {
	Docopt::new(usage).unwrap().parse().unwrap_or_else(|_| {
		writeln!(stderr(), "Invalid arguments.\n{}", usage); exit(-1);
	})
}

// Helper methods for error reporting
trait ErrorHelper<T> {
	fn on_error(self, msg: &str) -> T;
}

impl<T> ErrorHelper<T> for Option<T> {
	fn on_error(self, msg: &str) -> T {
		match self {
			Some(x) => x,
			None => { eprintln!("ERROR: {}\n", msg); exit(-1) }
		}
	}
}

impl<T, E> ErrorHelper<T> for Result<T, E> {
	fn on_error(self, msg: &str) -> T {
		match self {
			Ok(x) => x,
			Err(_) => { eprintln!("ERROR: {}\n", msg); exit(-1) }
		}
	}
}
