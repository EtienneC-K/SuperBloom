//! Public Rust library API for building and querying a SuperBloom index.
//!
//! This API is intentionally explicit: all parameters are provided directly in
//! `SuperBloomConfig` (no presets), then sequences/FASTA records are added with
//! `add_sequence` / `add_fasta`.
//!
//! ```no_run
//! use bloomybloom::{MinimizerMode, SuperBloom, SuperBloomConfig};
//!
//! let config = SuperBloomConfig {
//!     k: 31,
//!     m: 13,
//!     s: 21,
//!     n_hashes: 3,
//!     size_exponent: 20,
//!     block_size_exponent: 10,
//!     minimizer_mode: MinimizerMode::Simd,
//! };
//!
//! let mut bloom = SuperBloom::new(config)?;
//! let _added = bloom.add_sequence(b"ACGTACGTACGTACGTACGTACGTACGTACGT")?;
//! let frozen = bloom.into_frozen();
//! let hits = frozen.query_sequence(b"ACGTACGTACGTACGTACGTACGTACGTACGT");
//! assert!(!hits.is_empty());
//! # Ok::<(), bloomybloom::SuperBloomError>(())
//! ```

mod bloom;
pub mod decyclers;
pub mod minimizers;
pub mod super_bitvec;
pub mod utils;

use bloom::{BloomFilter, FrozenBloomFilter};
use decyclers::Decycler;
pub use minimizers::MinimizerMode;
use minimizers::selected_mins_x_pos;
use needletail::{FastxReader, parse_fastx_file};
use packed_seq::{PackedSeq, PackedSeqVec, Seq, SeqVec};
use rayon::prelude::*;
use rayon::{ThreadPool, ThreadPoolBuilder};
use std::fmt::{Display, Formatter};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::Path;
use std::sync::{Once, RwLock};
use utils::{hash_u128_to_u64, sum_vec_bool, xorshift_u64};

const SHARD_COUNT: usize = 1024;
const SUPERBLOOM_MAGIC: &[u8; 8] = b"SBLOOM01";
const PAR_BATCH_RECORDS: usize = 256;
const DEFAULT_THREADS: usize = 8;
static EXPERIMENTAL_MINIMIZER_WARNING_ONCE: Once = Once::new();

/// Full manual configuration for building a SuperBloom index.
///
/// All geometry values are explicit powers of two:
/// - `size_exponent`: bloom filter size in bits is `2^size_exponent`
/// - `block_size_exponent`: block size in bits is `2^block_size_exponent`
#[derive(Clone, Copy, Debug)]
pub struct SuperBloomConfig {
    pub k: u16,
    pub m: u16,
    pub s: u16,
    pub n_hashes: usize,
    pub size_exponent: u8,
    pub block_size_exponent: u8,
    pub minimizer_mode: MinimizerMode,
}

impl Default for SuperBloomConfig {
    fn default() -> Self {
        Self {
            k: 31,
            m: 21,
            s: 27,
            n_hashes: 8,
            size_exponent: 35,
            block_size_exponent: 9,
            minimizer_mode: MinimizerMode::Simd,
        }
    }
}

/// Summary returned by `add_fasta`.
#[derive(Clone, Copy, Debug, Default)]
pub struct AddReport {
    pub records_processed: u64,
    pub records_indexed: u64,
    pub kmers_added: u64,
}

/// Summary returned by `query_fasta`.
#[derive(Clone, Copy, Debug, Default)]
pub struct QueryReport {
    pub records_processed: u64,
    pub queried_kmers: u64,
    pub positive_kmers: u64,
}

/// Library-level errors returned by configuration, indexing, and querying APIs.
#[derive(Debug)]
pub enum SuperBloomError {
    InvalidConfig(String),
    ThreadPoolBuild(String),
    FastaOpen { path: String, message: String },
    FastaRead { path: String, message: String },
    Serialization { path: String, message: String },
    Deserialization { path: String, message: String },
    WrongMode(String),
    InternalState(String),
    PoisonedLock,
}

impl Display for SuperBloomError {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            SuperBloomError::InvalidConfig(msg) => write!(f, "invalid config: {msg}"),
            SuperBloomError::ThreadPoolBuild(msg) => {
                write!(f, "failed to build rayon thread pool: {msg}")
            }
            SuperBloomError::FastaOpen { path, message } => {
                write!(f, "failed to open FASTA/FASTQ '{path}': {message}")
            }
            SuperBloomError::FastaRead { path, message } => {
                write!(
                    f,
                    "failed to read FASTA/FASTQ record from '{path}': {message}"
                )
            }
            SuperBloomError::Serialization { path, message } => {
                write!(f, "failed to serialize superbloom to '{path}': {message}")
            }
            SuperBloomError::Deserialization { path, message } => {
                write!(
                    f,
                    "failed to deserialize superbloom from '{path}': {message}"
                )
            }
            SuperBloomError::WrongMode(message) => write!(f, "{message}"),
            SuperBloomError::InternalState(msg) => write!(f, "internal state error: {msg}"),
            SuperBloomError::PoisonedLock => write!(f, "internal mutex lock is poisoned"),
        }
    }
}

impl std::error::Error for SuperBloomError {}

fn warn_if_experimental_minimizer(mode: MinimizerMode) {
    if matches!(mode, MinimizerMode::Simd) {
        return;
    }
    EXPERIMENTAL_MINIMIZER_WARNING_ONCE.call_once(|| {
        eprintln!(
            "Warning: non-SIMD minimizer modes are experimental (selected: {:?}). \
             Prefer MinimizerMode::Simd for production workloads.",
            mode
        );
    });
}

enum BloomState {
    Mutable(BloomFilter),
    Frozen(FrozenBloomFilter),
}

/// SuperBloom index that can switch in-place between mutable (insert) and frozen (query) modes.
///
/// - `freeze()`: move shards from `Mutex<BlockShard>` to immutable `BlockShard` (no bit-array copy).
/// - `thaw()`: rebuild mutex wrappers for further insertions (no bit-array copy).
pub struct SuperBloom {
    state: RwLock<Option<BloomState>>,
    decycler: Decycler,
    config: SuperBloomConfig,
    inserted_kmers: u64,
    thread_pool: Option<ThreadPool>,
    thread_count: Option<usize>,
}

/// Read-only index for querying.
pub struct FrozenSuperBloom {
    bloom: FrozenBloomFilter,
    decycler: Decycler,
    config: SuperBloomConfig,
    inserted_kmers: u64,
    thread_pool: Option<ThreadPool>,
    thread_count: Option<usize>,
}

impl SuperBloom {
    /// Create an empty SuperBloom index from a full manual config.
    pub fn new(config: SuperBloomConfig) -> Result<Self, SuperBloomError> {
        let (size_bits, block_size_bits, nb_blocks) = resolve_geometry(config)?;
        warn_if_experimental_minimizer(config.minimizer_mode);

        let mut decycler = if matches!(config.minimizer_mode, MinimizerMode::Decycling) {
            Decycler::new(config.m)
        } else {
            Decycler::new(1)
        };
        if matches!(config.minimizer_mode, MinimizerMode::Decycling) {
            decycler.compute_blocks();
        }

        let bloom = BloomFilter::new(
            size_bits,
            config.n_hashes,
            config.k as usize,
            block_size_bits,
            nb_blocks,
        );
        let pool = ThreadPoolBuilder::new()
            .num_threads(DEFAULT_THREADS)
            .build()
            .map_err(|err| SuperBloomError::ThreadPoolBuild(err.to_string()))?;

        Ok(Self {
            state: RwLock::new(Some(BloomState::Mutable(bloom))),
            decycler,
            config,
            inserted_kmers: 0,
            thread_pool: Some(pool),
            thread_count: Some(DEFAULT_THREADS),
        })
    }

    /// Configure a dedicated Rayon thread pool for parallel indexing work.
    ///
    /// This controls the thread count used by `add_fasta`.
    /// New indexes start with a default of 8 threads.
    pub fn set_threads(&mut self, threads: usize) -> Result<(), SuperBloomError> {
        if threads == 0 {
            return Err(SuperBloomError::InvalidConfig(
                "threads must be >= 1".to_string(),
            ));
        }
        let pool = ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .map_err(|err| SuperBloomError::ThreadPoolBuild(err.to_string()))?;
        self.thread_pool = Some(pool);
        self.thread_count = Some(threads);
        Ok(())
    }

    /// Disable the dedicated thread pool and fall back to Rayon global settings.
    pub fn clear_threads(&mut self) {
        self.thread_pool = None;
        self.thread_count = None;
    }

    /// Return the user-configured thread count for indexing, if any.
    pub fn threads(&self) -> Option<usize> {
        self.thread_count
    }

    /// Switch to frozen query mode.
    ///
    /// This transition moves shards in-place and does not duplicate bit arrays.
    pub fn freeze(&self) {
        self.ensure_frozen();
    }

    /// Switch back to mutable insert mode.
    ///
    /// This transition rebuilds mutex wrappers in-place and does not duplicate bit arrays.
    pub fn thaw(&self) {
        self.ensure_mutable();
    }

    /// Returns `true` when the index is currently in frozen query mode.
    pub fn is_frozen(&self) -> bool {
        let guard = self.state.read().unwrap_or_else(|err| err.into_inner());
        matches!(guard.as_ref(), Some(BloomState::Frozen(_)))
    }

    /// Add one DNA sequence (ASCII bytes) to the index.
    ///
    /// Returns the number of k-mers added from this sequence.
    pub fn add_sequence(&mut self, sequence: &[u8]) -> Result<u64, SuperBloomError> {
        if sequence.len() < self.config.k as usize {
            return Ok(0);
        }

        self.ensure_mutable();
        let packed = PackedSeqVec::from_ascii(sequence);
        let added = self.with_mutable_bloom(|bloom| {
            insert_packed_sequence(bloom, &self.decycler, self.config, packed)
        })?;
        self.inserted_kmers = self.inserted_kmers.saturating_add(added);
        Ok(added)
    }

    /// Add all records from a FASTA/FASTQ file to the index.
    pub fn add_fasta<P: AsRef<Path>>(&mut self, path: P) -> Result<AddReport, SuperBloomError> {
        let mut report = AddReport::default();
        let path_ref = path.as_ref();
        let path_string = path_ref.display().to_string();
        let mut reader = open_fastx_reader(path_ref, &path_string)?;
        let mut batch: Vec<Vec<u8>> = Vec::with_capacity(PAR_BATCH_RECORDS);
        self.ensure_mutable();

        while let Some(record_result) = reader.next() {
            let record = record_result.map_err(|err| SuperBloomError::FastaRead {
                path: path_string.clone(),
                message: err.to_string(),
            })?;
            report.records_processed = report.records_processed.saturating_add(1);
            batch.push(record.seq().as_ref().to_vec());

            if batch.len() >= PAR_BATCH_RECORDS {
                let (records_indexed, kmers_added) = self.index_batch(&batch)?;
                report.records_indexed = report.records_indexed.saturating_add(records_indexed);
                report.kmers_added = report.kmers_added.saturating_add(kmers_added);
                batch.clear();
            }
        }

        if !batch.is_empty() {
            let (records_indexed, kmers_added) = self.index_batch(&batch)?;
            report.records_indexed = report.records_indexed.saturating_add(records_indexed);
            report.kmers_added = report.kmers_added.saturating_add(kmers_added);
        }

        self.inserted_kmers = self.inserted_kmers.saturating_add(report.kmers_added);
        Ok(report)
    }

    fn index_batch(&self, batch: &[Vec<u8>]) -> Result<(u64, u64), SuperBloomError> {
        self.with_mutable_bloom(|bloom| match &self.thread_pool {
            Some(pool) => pool
                .install(|| Self::index_batch_parallel(bloom, &self.decycler, self.config, batch)),
            None => Self::index_batch_parallel(bloom, &self.decycler, self.config, batch),
        })
    }

    fn index_batch_parallel(
        bloom: &BloomFilter,
        decycler: &Decycler,
        config: SuperBloomConfig,
        batch: &[Vec<u8>],
    ) -> Result<(u64, u64), SuperBloomError> {
        batch
            .par_iter()
            .try_fold(
                || (0u64, 0u64),
                |(records_indexed, kmers_added), sequence| {
                    if sequence.len() < config.k as usize {
                        return Ok((records_indexed, kmers_added));
                    }

                    let packed = PackedSeqVec::from_ascii(sequence);
                    let added = insert_packed_sequence(bloom, decycler, config, packed)?;
                    if added > 0 {
                        Ok((
                            records_indexed.saturating_add(1),
                            kmers_added.saturating_add(added),
                        ))
                    } else {
                        Ok((records_indexed, kmers_added))
                    }
                },
            )
            .try_reduce(
                || (0u64, 0u64),
                |(left_records, left_kmers), (right_records, right_kmers)| {
                    Ok((
                        left_records.saturating_add(right_records),
                        left_kmers.saturating_add(right_kmers),
                    ))
                },
            )
    }

    /// Query one DNA sequence while in frozen mode.
    ///
    /// Returns one boolean per k-mer window.
    pub fn query_sequence(&self, sequence: &[u8]) -> Result<Vec<bool>, SuperBloomError> {
        self.ensure_frozen();
        self.with_frozen_bloom(|frozen| {
            Ok(run_query_sequence(
                frozen,
                &self.decycler,
                self.config,
                sequence,
            ))
        })
    }

    /// Query every record from a FASTA/FASTQ file while in frozen mode.
    pub fn query_fasta<P: AsRef<Path>>(&self, path: P) -> Result<QueryReport, SuperBloomError> {
        self.ensure_frozen();
        self.with_frozen_bloom(|frozen| {
            run_query_fasta(
                frozen,
                &self.decycler,
                self.config,
                path,
                self.thread_pool.as_ref(),
            )
        })
    }

    /// Serialize current index by freezing if needed.
    pub fn save<P: AsRef<Path>>(&self, path: P) -> Result<(), SuperBloomError> {
        self.ensure_frozen();
        self.with_frozen_bloom(|frozen| {
            save_frozen_components(path, frozen, self.config, self.inserted_kmers)
        })
    }

    /// Load a serialized index in frozen query mode.
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self, SuperBloomError> {
        let frozen = FrozenSuperBloom::load(path)?;
        Ok(Self {
            state: RwLock::new(Some(BloomState::Frozen(frozen.bloom))),
            decycler: frozen.decycler,
            config: frozen.config,
            inserted_kmers: frozen.inserted_kmers,
            thread_pool: frozen.thread_pool,
            thread_count: frozen.thread_count,
        })
    }

    /// Number of k-mers inserted so far.
    pub fn inserted_kmers(&self) -> u64 {
        self.inserted_kmers
    }

    /// Access the manual build configuration.
    pub fn config(&self) -> &SuperBloomConfig {
        &self.config
    }

    fn with_mutable_bloom<T>(
        &self,
        work: impl FnOnce(&BloomFilter) -> Result<T, SuperBloomError>,
    ) -> Result<T, SuperBloomError> {
        let guard = self.state.read().unwrap_or_else(|err| err.into_inner());
        match guard.as_ref() {
            Some(BloomState::Mutable(bloom)) => work(bloom),
            Some(BloomState::Frozen(_)) => Err(SuperBloomError::WrongMode(
                "index is in frozen mode unexpectedly".to_string(),
            )),
            None => Err(SuperBloomError::InternalState(
                "bloom state is unexpectedly empty".to_string(),
            )),
        }
    }

    fn with_frozen_bloom<T>(
        &self,
        work: impl FnOnce(&FrozenBloomFilter) -> Result<T, SuperBloomError>,
    ) -> Result<T, SuperBloomError> {
        let guard = self.state.read().unwrap_or_else(|err| err.into_inner());
        match guard.as_ref() {
            Some(BloomState::Frozen(bloom)) => work(bloom),
            Some(BloomState::Mutable(_)) => Err(SuperBloomError::WrongMode(
                "index is in mutable mode unexpectedly".to_string(),
            )),
            None => Err(SuperBloomError::InternalState(
                "bloom state is unexpectedly empty".to_string(),
            )),
        }
    }

    /// Freeze this mutable index into a read-only query index.
    pub fn into_frozen(self) -> FrozenSuperBloom {
        let state_opt = self
            .state
            .into_inner()
            .unwrap_or_else(|err| err.into_inner());
        let state = state_opt.expect("internal state should be present");
        let bloom = match state {
            BloomState::Frozen(frozen) => frozen,
            BloomState::Mutable(mutable) => mutable.into_frozen(),
        };
        FrozenSuperBloom {
            bloom,
            decycler: self.decycler,
            config: self.config,
            inserted_kmers: self.inserted_kmers,
            thread_pool: self.thread_pool,
            thread_count: self.thread_count,
        }
    }

    fn ensure_frozen(&self) {
        let needs_freeze = {
            let guard = self.state.read().unwrap_or_else(|err| err.into_inner());
            matches!(guard.as_ref(), Some(BloomState::Mutable(_)))
        };
        if !needs_freeze {
            return;
        }
        let mut guard = self.state.write().unwrap_or_else(|err| err.into_inner());
        let old_state = guard.take().expect("internal state should be present");
        *guard = Some(match old_state {
            BloomState::Mutable(mutable) => BloomState::Frozen(mutable.into_frozen()),
            BloomState::Frozen(frozen) => BloomState::Frozen(frozen),
        });
    }

    fn ensure_mutable(&self) {
        let needs_thaw = {
            let guard = self.state.read().unwrap_or_else(|err| err.into_inner());
            matches!(guard.as_ref(), Some(BloomState::Frozen(_)))
        };
        if !needs_thaw {
            return;
        }
        let mut guard = self.state.write().unwrap_or_else(|err| err.into_inner());
        let old_state = guard.take().expect("internal state should be present");
        *guard = Some(match old_state {
            BloomState::Frozen(frozen) => BloomState::Mutable(frozen.into_mutable()),
            BloomState::Mutable(mutable) => BloomState::Mutable(mutable),
        });
    }
}

impl FrozenSuperBloom {
    /// Rebuild mutable wrappers to allow further insertions.
    ///
    /// This transition moves shards in-place and does not duplicate bit arrays.
    pub fn into_mutable(self) -> SuperBloom {
        SuperBloom {
            state: RwLock::new(Some(BloomState::Mutable(self.bloom.into_mutable()))),
            decycler: self.decycler,
            config: self.config,
            inserted_kmers: self.inserted_kmers,
            thread_pool: self.thread_pool,
            thread_count: self.thread_count,
        }
    }

    /// Serialize this frozen index to a binary file.
    pub fn save<P: AsRef<Path>>(&self, path: P) -> Result<(), SuperBloomError> {
        save_frozen_components(path, &self.bloom, self.config, self.inserted_kmers)
    }

    /// Load a frozen index from a file previously written with `save`.
    pub fn load<P: AsRef<Path>>(path: P) -> Result<Self, SuperBloomError> {
        let path_ref = path.as_ref();
        let path_string = path_ref.display().to_string();
        let file = File::open(path_ref).map_err(|err| SuperBloomError::Deserialization {
            path: path_string.clone(),
            message: err.to_string(),
        })?;
        let mut reader = BufReader::new(file);

        let mut magic = [0u8; 8];
        reader
            .read_exact(&mut magic)
            .map_err(|err| SuperBloomError::Deserialization {
                path: path_string.clone(),
                message: err.to_string(),
            })?;
        if &magic != SUPERBLOOM_MAGIC {
            return Err(SuperBloomError::Deserialization {
                path: path_string,
                message: "invalid file signature".to_string(),
            });
        }

        let config = read_config(&mut reader).map_err(|err| SuperBloomError::Deserialization {
            path: path_ref.display().to_string(),
            message: err.to_string(),
        })?;
        warn_if_experimental_minimizer(config.minimizer_mode);
        let mut inserted_buf = [0u8; 8];
        reader
            .read_exact(&mut inserted_buf)
            .map_err(|err| SuperBloomError::Deserialization {
                path: path_ref.display().to_string(),
                message: err.to_string(),
            })?;
        let inserted_kmers = u64::from_le_bytes(inserted_buf);

        let bloom = FrozenBloomFilter::read_from(&mut reader).map_err(|err| {
            SuperBloomError::Deserialization {
                path: path_ref.display().to_string(),
                message: err.to_string(),
            }
        })?;

        let mut decycler = if matches!(config.minimizer_mode, MinimizerMode::Decycling) {
            Decycler::new(config.m)
        } else {
            Decycler::new(1)
        };
        if matches!(config.minimizer_mode, MinimizerMode::Decycling) {
            decycler.compute_blocks();
        }
        let pool = ThreadPoolBuilder::new()
            .num_threads(DEFAULT_THREADS)
            .build()
            .map_err(|err| SuperBloomError::ThreadPoolBuild(err.to_string()))?;

        Ok(Self {
            bloom,
            decycler,
            config,
            inserted_kmers,
            thread_pool: Some(pool),
            thread_count: Some(DEFAULT_THREADS),
        })
    }

    /// Query one DNA sequence.
    ///
    /// The returned vector has one boolean per k-mer window.
    pub fn query_sequence(&self, sequence: &[u8]) -> Vec<bool> {
        run_query_sequence(&self.bloom, &self.decycler, self.config, sequence)
    }

    /// Query every record from a FASTA/FASTQ file.
    pub fn query_fasta<P: AsRef<Path>>(&self, path: P) -> Result<QueryReport, SuperBloomError> {
        run_query_fasta(
            &self.bloom,
            &self.decycler,
            self.config,
            path,
            self.thread_pool.as_ref(),
        )
    }

    /// Number of k-mers inserted before freezing.
    pub fn inserted_kmers(&self) -> u64 {
        self.inserted_kmers
    }

    /// Access the manual build configuration.
    pub fn config(&self) -> &SuperBloomConfig {
        &self.config
    }
}

fn open_fastx_reader(
    path: &Path,
    path_string: &str,
) -> Result<Box<dyn FastxReader>, SuperBloomError> {
    parse_fastx_file(path).map_err(|err| SuperBloomError::FastaOpen {
        path: path_string.to_string(),
        message: err.to_string(),
    })
}

fn run_query_sequence(
    bloom: &FrozenBloomFilter,
    decycler: &Decycler,
    config: SuperBloomConfig,
    sequence: &[u8],
) -> Vec<bool> {
    if sequence.len() < config.k as usize {
        return Vec::new();
    }

    let packed = PackedSeqVec::from_ascii(sequence);
    if config.s <= 31 {
        bloom.check_sequence(
            packed,
            config.k,
            config.m,
            config.s,
            decycler,
            config.minimizer_mode,
        )
    } else {
        bloom.check_sequence_u128(
            packed,
            config.k,
            config.m,
            config.s,
            decycler,
            config.minimizer_mode,
        )
    }
}

fn run_query_fasta<P: AsRef<Path>>(
    bloom: &FrozenBloomFilter,
    decycler: &Decycler,
    config: SuperBloomConfig,
    path: P,
    thread_pool: Option<&ThreadPool>,
) -> Result<QueryReport, SuperBloomError> {
    let mut report = QueryReport::default();
    let path_ref = path.as_ref();
    let path_string = path_ref.display().to_string();
    let mut reader = open_fastx_reader(path_ref, &path_string)?;
    let mut batch: Vec<Vec<u8>> = Vec::with_capacity(PAR_BATCH_RECORDS);

    while let Some(record_result) = reader.next() {
        let record = record_result.map_err(|err| SuperBloomError::FastaRead {
            path: path_string.clone(),
            message: err.to_string(),
        })?;
        batch.push(record.seq().as_ref().to_vec());
        if batch.len() >= PAR_BATCH_RECORDS {
            let (records_processed, queried_kmers, positive_kmers) =
                query_batch_stats(bloom, decycler, config, &batch, thread_pool);
            report.records_processed = report.records_processed.saturating_add(records_processed);
            report.queried_kmers = report.queried_kmers.saturating_add(queried_kmers);
            report.positive_kmers = report.positive_kmers.saturating_add(positive_kmers);
            batch.clear();
        }
    }

    if !batch.is_empty() {
        let (records_processed, queried_kmers, positive_kmers) =
            query_batch_stats(bloom, decycler, config, &batch, thread_pool);
        report.records_processed = report.records_processed.saturating_add(records_processed);
        report.queried_kmers = report.queried_kmers.saturating_add(queried_kmers);
        report.positive_kmers = report.positive_kmers.saturating_add(positive_kmers);
    }

    Ok(report)
}

fn query_batch_stats(
    bloom: &FrozenBloomFilter,
    decycler: &Decycler,
    config: SuperBloomConfig,
    batch: &[Vec<u8>],
    thread_pool: Option<&ThreadPool>,
) -> (u64, u64, u64) {
    match thread_pool {
        Some(pool) => pool.install(|| query_batch_stats_parallel(bloom, decycler, config, batch)),
        None => query_batch_stats_parallel(bloom, decycler, config, batch),
    }
}

fn query_batch_stats_parallel(
    bloom: &FrozenBloomFilter,
    decycler: &Decycler,
    config: SuperBloomConfig,
    batch: &[Vec<u8>],
) -> (u64, u64, u64) {
    batch
        .par_iter()
        .fold(
            || (0u64, 0u64, 0u64),
            |(records_processed, queried_kmers, positive_kmers), sequence| {
                let hits = run_query_sequence(bloom, decycler, config, sequence);
                (
                    records_processed.saturating_add(1),
                    queried_kmers.saturating_add(hits.len() as u64),
                    positive_kmers.saturating_add(sum_vec_bool(&hits) as u64),
                )
            },
        )
        .reduce(
            || (0u64, 0u64, 0u64),
            |(r1, q1, p1), (r2, q2, p2)| {
                (
                    r1.saturating_add(r2),
                    q1.saturating_add(q2),
                    p1.saturating_add(p2),
                )
            },
        )
}

fn save_frozen_components<P: AsRef<Path>>(
    path: P,
    bloom: &FrozenBloomFilter,
    config: SuperBloomConfig,
    inserted_kmers: u64,
) -> Result<(), SuperBloomError> {
    let path_ref = path.as_ref();
    let path_string = path_ref.display().to_string();
    let file = File::create(path_ref).map_err(|err| SuperBloomError::Serialization {
        path: path_string.clone(),
        message: err.to_string(),
    })?;
    let mut writer = BufWriter::new(file);

    writer
        .write_all(SUPERBLOOM_MAGIC)
        .map_err(|err| SuperBloomError::Serialization {
            path: path_string.clone(),
            message: err.to_string(),
        })?;
    write_config(&mut writer, config).map_err(|err| SuperBloomError::Serialization {
        path: path_string.clone(),
        message: err.to_string(),
    })?;
    writer
        .write_all(&inserted_kmers.to_le_bytes())
        .map_err(|err| SuperBloomError::Serialization {
            path: path_string.clone(),
            message: err.to_string(),
        })?;
    bloom
        .write_to(&mut writer)
        .map_err(|err| SuperBloomError::Serialization {
            path: path_string.clone(),
            message: err.to_string(),
        })?;
    writer
        .flush()
        .map_err(|err| SuperBloomError::Serialization {
            path: path_string,
            message: err.to_string(),
        })?;
    Ok(())
}

fn write_config<W: Write>(writer: &mut W, config: SuperBloomConfig) -> Result<(), std::io::Error> {
    writer.write_all(&config.k.to_le_bytes())?;
    writer.write_all(&config.m.to_le_bytes())?;
    writer.write_all(&config.s.to_le_bytes())?;
    write_u64(writer, config.n_hashes as u64)?;
    writer.write_all(&[config.size_exponent])?;
    writer.write_all(&[config.block_size_exponent])?;
    match config.minimizer_mode {
        MinimizerMode::Simd => {
            writer.write_all(&[0u8])?;
            writer.write_all(&0u16.to_le_bytes())?;
        }
        MinimizerMode::Decycling => {
            writer.write_all(&[1u8])?;
            writer.write_all(&0u16.to_le_bytes())?;
        }
        MinimizerMode::OpenClosed { t } => {
            writer.write_all(&[2u8])?;
            writer.write_all(&t.to_le_bytes())?;
        }
    }
    Ok(())
}

fn read_config<R: Read>(reader: &mut R) -> Result<SuperBloomConfig, std::io::Error> {
    let mut k_buf = [0u8; 2];
    let mut m_buf = [0u8; 2];
    let mut s_buf = [0u8; 2];
    reader.read_exact(&mut k_buf)?;
    reader.read_exact(&mut m_buf)?;
    reader.read_exact(&mut s_buf)?;
    let k = u16::from_le_bytes(k_buf);
    let m = u16::from_le_bytes(m_buf);
    let s = u16::from_le_bytes(s_buf);

    let n_hashes = read_u64(reader)?;
    let n_hashes = usize::try_from(n_hashes).map_err(|_| {
        std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            "n_hashes does not fit on this platform",
        )
    })?;

    let mut size_exp = [0u8; 1];
    let mut block_exp = [0u8; 1];
    let mut mode_buf = [0u8; 1];
    let mut t_buf = [0u8; 2];
    reader.read_exact(&mut size_exp)?;
    reader.read_exact(&mut block_exp)?;
    reader.read_exact(&mut mode_buf)?;
    reader.read_exact(&mut t_buf)?;
    let t = u16::from_le_bytes(t_buf);

    let minimizer_mode = match mode_buf[0] {
        0 => MinimizerMode::Simd,
        1 => MinimizerMode::Decycling,
        2 => MinimizerMode::OpenClosed { t },
        _ => {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "unknown minimizer mode tag",
            ));
        }
    };

    let config = SuperBloomConfig {
        k,
        m,
        s,
        n_hashes,
        size_exponent: size_exp[0],
        block_size_exponent: block_exp[0],
        minimizer_mode,
    };
    resolve_geometry(config)
        .map_err(|err| std::io::Error::new(std::io::ErrorKind::InvalidData, err.to_string()))?;
    Ok(config)
}

fn write_u64<W: Write>(writer: &mut W, value: u64) -> Result<(), std::io::Error> {
    writer.write_all(&value.to_le_bytes())
}

fn read_u64<R: Read>(reader: &mut R) -> Result<u64, std::io::Error> {
    let mut buf = [0u8; 8];
    reader.read_exact(&mut buf)?;
    Ok(u64::from_le_bytes(buf))
}

fn resolve_geometry(config: SuperBloomConfig) -> Result<(usize, usize, usize), SuperBloomError> {
    if config.n_hashes == 0 {
        return Err(SuperBloomError::InvalidConfig(
            "n_hashes must be greater than 0".to_string(),
        ));
    }
    if config.k == 0 {
        return Err(SuperBloomError::InvalidConfig(
            "k must be greater than 0".to_string(),
        ));
    }
    if config.m == 0 || config.m > config.k || config.m >= 32 {
        return Err(SuperBloomError::InvalidConfig(
            "m must satisfy: 1 <= m <= k and m < 32".to_string(),
        ));
    }
    if config.s == 0 || config.s > config.k || config.s >= 62 {
        return Err(SuperBloomError::InvalidConfig(
            "s must satisfy: 1 <= s <= k and s < 62".to_string(),
        ));
    }
    if let MinimizerMode::OpenClosed { t } = config.minimizer_mode
        && (t == 0 || t > config.m)
    {
        return Err(SuperBloomError::InvalidConfig(
            "open-closed minimizer requires 1 <= t <= m".to_string(),
        ));
    }
    if config.block_size_exponent > config.size_exponent {
        return Err(SuperBloomError::InvalidConfig(
            "block_size_exponent cannot be greater than size_exponent".to_string(),
        ));
    }

    let size_bits = 1usize
        .checked_shl(config.size_exponent as u32)
        .ok_or_else(|| {
            SuperBloomError::InvalidConfig(
                "size_exponent is too large for this platform".to_string(),
            )
        })?;
    let block_size_bits = 1usize
        .checked_shl(config.block_size_exponent as u32)
        .ok_or_else(|| {
            SuperBloomError::InvalidConfig(
                "block_size_exponent is too large for this platform".to_string(),
            )
        })?;
    let nb_blocks = size_bits / block_size_bits;

    if nb_blocks < SHARD_COUNT {
        return Err(SuperBloomError::InvalidConfig(format!(
            "size/block_size must be at least {SHARD_COUNT} blocks"
        )));
    }
    if nb_blocks % SHARD_COUNT != 0 {
        return Err(SuperBloomError::InvalidConfig(format!(
            "size/block_size must be a multiple of {SHARD_COUNT} blocks"
        )));
    }

    Ok((size_bits, block_size_bits, nb_blocks))
}

fn insert_packed_sequence(
    bloom: &BloomFilter,
    decycler: &Decycler,
    config: SuperBloomConfig,
    sequence: PackedSeqVec,
) -> Result<u64, SuperBloomError> {
    let total_kmers = sequence.len() + 1 - config.k as usize;
    let (super_kmers_positions, minimizer_values, packed_sequence) = selected_mins_x_pos(
        sequence,
        config.k,
        config.m,
        decycler,
        config.minimizer_mode,
    );

    if super_kmers_positions.len() != minimizer_values.len() {
        return Err(SuperBloomError::InternalState(
            "super-kmer boundaries and minimizer values do not align".to_string(),
        ));
    }

    for i in 0..super_kmers_positions.len() {
        let start_kmer = super_kmers_positions[i] as usize;
        let end_kmer = if i + 1 < super_kmers_positions.len() {
            super_kmers_positions[i + 1] as usize
        } else {
            packed_sequence.len() + 1 - config.k as usize
        };
        let hashed_minimizer = xorshift_u64(minimizer_values[i]) & (bloom.nb_blocks as u64 - 1);

        if config.s <= 31 {
            insert_super_kmer_u64(
                bloom,
                &packed_sequence,
                start_kmer,
                end_kmer,
                config.k,
                config.s,
                hashed_minimizer,
            )?;
        } else {
            insert_super_kmer_u128(
                bloom,
                &packed_sequence,
                start_kmer,
                end_kmer,
                config.k,
                config.s,
                hashed_minimizer,
            )?;
        }
    }

    Ok(total_kmers as u64)
}

fn insert_super_kmer_u64(
    bloom: &BloomFilter,
    sequence: &PackedSeqVec,
    start_kmer: usize,
    end_kmer: usize,
    k: u16,
    s: u16,
    hashed_minimizer: u64,
) -> Result<(), SuperBloomError> {
    let blocknum = (hashed_minimizer as usize) & (SHARD_COUNT - 1);
    let subblocknum = ((hashed_minimizer as usize) >> 10) & ((bloom.nb_blocks >> 10) - 1);
    let mut block = bloom.filter[blocknum]
        .lock()
        .map_err(|_| SuperBloomError::PoisonedLock)?;

    let end_smer = end_kmer + (k - s) as usize;
    for i in start_kmer..end_smer {
        let smer: PackedSeq = sequence.slice(i..i + s as usize);
        let mut hash = xorshift_u64(smer.as_u64());
        for _ in 0..bloom.n_hashes {
            let address = hash as usize & (bloom.block_size - 1);
            if !block.get(subblocknum, address) {
                block.set(subblocknum, address, true);
            }
            hash = xorshift_u64(hash);
        }
    }

    Ok(())
}

fn insert_super_kmer_u128(
    bloom: &BloomFilter,
    sequence: &PackedSeqVec,
    start_kmer: usize,
    end_kmer: usize,
    k: u16,
    s: u16,
    hashed_minimizer: u64,
) -> Result<(), SuperBloomError> {
    let blocknum = (hashed_minimizer as usize) & (SHARD_COUNT - 1);
    let subblocknum = ((hashed_minimizer as usize) >> 10) & ((bloom.nb_blocks >> 10) - 1);
    let mut block = bloom.filter[blocknum]
        .lock()
        .map_err(|_| SuperBloomError::PoisonedLock)?;

    let end_smer = end_kmer + (k - s) as usize;
    for i in start_kmer..end_smer {
        let smer: PackedSeq = sequence.slice(i..i + s as usize);
        let mut hash = hash_u128_to_u64(smer.as_u128());
        for _ in 0..bloom.n_hashes {
            let address = hash as usize & (bloom.block_size - 1);
            if !block.get(subblocknum, address) {
                block.set(subblocknum, address, true);
            }
            hash = xorshift_u64(hash);
        }
    }

    Ok(())
}
