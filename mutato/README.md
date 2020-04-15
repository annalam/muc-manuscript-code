# Mutato - A variant calling framework implemented in Rust

Repository currently contains two different implementations for variant calling. The first implementation works on top of samtools mpileup, and simply parses its output. The second (in-progress) implementation does not depend on samtools, but instead performs variant calling based on BAM files directly.

## Features
- High performance (analyzes xxx,xxx Mb of sequencing data / hour / thread)
- Self-contained statically linked binary, no external dependencies
- Simultaneous variant calling across multiple BAM files
- Implemented entirely in the Rust language (no unsafe code)

Installation
------------

Install Rust (version 1.31 or later). Then run the following command:
```
cargo install --force --git https://github.com/annalam/mutato
```
