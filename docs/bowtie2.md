# Bowtie2

## Installation

1. basic automatic dependency management and static linkage of `zstd` and `zlib`
   - `make static-libs && make STATIC_BUILD=1`
2. SRA (Sequence Read Archive) support
   - `make sra-deps && make USE_SRA=1.`
3. libsais support: a state-of-the-art suffix array construction algorithm that speeds up the building of the index
   - `[g]make libsais USE_SAIS_OPENMP=1 *`
   - this is important so it can be multithreaded

To build bowtie2-build with libsais first make sure
that the libsais submodule is available.

```sh
# first time cloning
git clone --recursive https://github.com/BenLangmead/bowtie2.git

# existing checkout of bowtie2
git submodule init && git submodule update
```

### Building with CMake

To build Bowtie2 with SRA and libsais support:

```sh
cmake . -D USE_SRA=1 -D USE_SAIS=1 && cmake --build .
```
