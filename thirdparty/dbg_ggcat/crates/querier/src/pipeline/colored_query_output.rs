use crate::structs::query_colored_counters::{ColorsRange, QueryColoredCountersSerializer};
use crate::ColoredQueryOutputFormat;
use colors::colors_manager::{ColorMapReader, ColorsManager, ColorsMergeManager};
use config::{
    get_compression_level_info, get_memory_mode, ColorIndexType, SwapPriority, DEFAULT_PREFETCH_AMOUNT,
    KEEP_FILES, QUERIES_COUNT_MIN_BATCH,
};
use flate2::Compression;
use ggcat_logging::UnrecoverableErrorLogging;
use hashes::HashFunctionFactory;
use nightly_quirks::prelude::*;
use parallel_processor::buckets::readers::compressed_binary_reader::CompressedBinaryReader;
use parallel_processor::buckets::readers::BucketReader;
use parallel_processor::buckets::writers::compressed_binary_writer::CompressedBinaryWriter;
use parallel_processor::buckets::{LockFreeBucket, SingleBucket};
use parallel_processor::memory_fs::RemoveFileMode;
use parallel_processor::phase_times_monitor::PHASES_TIMES_MONITOR;
use parking_lot::{Condvar, Mutex};
use rayon::prelude::*;
use std::cmp::Reverse;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::ops::DerefMut;
use std::path::PathBuf;
use std::sync::atomic::{AtomicUsize, Ordering};

// =============================
// Optional JSONL output writer
// =============================

enum QueryOutputFileWriter {
    Plain(File),
    LZ4Compressed(lz4::Encoder<File>),
    GzipCompressed(flate2::write::GzEncoder<File>),
}

impl Write for QueryOutputFileWriter {
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        match self {
            QueryOutputFileWriter::Plain(w) => w.write(buf),
            QueryOutputFileWriter::LZ4Compressed(w) => w.write(buf),
            QueryOutputFileWriter::GzipCompressed(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> std::io::Result<()> {
        match self {
            QueryOutputFileWriter::Plain(w) => w.flush(),
            QueryOutputFileWriter::LZ4Compressed(w) => w.flush(),
            QueryOutputFileWriter::GzipCompressed(w) => w.flush(),
        }
    }
}

// =============================
// Species aggregation helpers
// =============================

#[derive(Clone)]
enum Decision {
    Species(String),
    Discard(&'static str, String),
}

#[inline]
fn species_of(name: &str) -> &str {
    name.splitn(2, '_').next().unwrap_or(name)
}

// Use the suffix after "_" as the strain identifier; retain the full string
// when no separator is present.
#[inline]
fn strain_id_of(name: &str) -> &str {
    match name.splitn(2, '_').nth(1) {
        Some(s) if !s.is_empty() => s,
        _ => name,
    }
}

// Return all candidates tied at the exact maximum support value.
#[inline]
fn winners_of_max_strict(matches: &[(String, f64)]) -> (f64, Vec<String>) {
    if matches.is_empty() {
        return (0.0, vec![]);
    }
    let mut best = 0.0f64;
    for &(_, v) in matches {
        if v > best {
            best = v;
        }
    }
    let wins = matches
        .iter()
        .filter(|(_, v)| *v == best) // strict equality
        .map(|(k, _)| k.clone())
        .collect::<Vec<_>>();
    (best, wins)
}

/// Return the decision and the set of strains tied at maximum support. The
/// winner set is meaningful only for retained reads.
#[inline]
fn decide_species_and_winners(matches: &[(String, f64)], t: f64) -> (Decision, Vec<String>) {
    let (rmax, winners) = winners_of_max_strict(matches);

    if rmax <= t {
        return (
            Decision::Discard("below_threshold", format!("Rmax={:.6} T={:.6}", rmax, t)),
            vec![],
        );
    }

    if winners.is_empty() {
        return (Decision::Discard("no_match", String::new()), vec![]);
    }

    // All maximum-supported strains must belong to one species.
    let mut sp_set: Vec<&str> = winners.iter().map(|w| species_of(w)).collect();
    sp_set.sort_unstable();
    sp_set.dedup();

    if sp_set.len() > 1 {
        return (
            Decision::Discard("multi_species", format!("winners_species={}", sp_set.join(","))),
            vec![],
        );
    }

    (Decision::Species(sp_set[0].to_string()), winners)
}

// ---------- Adaptive-threshold sampling ----------

const FIRST_N: usize = 50_000;
const PCT: f64 = 0.30; // p30
const THR_MIN: f64 = 1.0 / 3.0;
const THR_MAX: f64 = 2.0 / 3.0;

struct RmaxSampler {
    vals: Vec<f64>,
}

impl RmaxSampler {
    fn new() -> Self {
        Self {
            vals: Vec::with_capacity(FIRST_N),
        }
    }

    #[inline]
    fn try_push(&mut self, r: f64) {
        if self.vals.len() < FIRST_N {
            self.vals.push(r);
        }
    }

    #[inline]
    fn len(&self) -> usize {
        self.vals.len()
    }

    fn compute_threshold(&mut self) -> f64 {
        if self.vals.is_empty() {
            return THR_MIN;
        }

        self.vals.sort_by(|a, b| {
            a.partial_cmp(b)
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        let n = self.vals.len();
        let pos = PCT * (n as f64 - 1.0);
        let i = pos.floor() as usize;
        let j = pos.ceil() as usize;

        let p = if i == j {
            self.vals[i]
        } else {
            let frac = pos - (i as f64);
            self.vals[i] * (1.0 - frac) + self.vals[j] * frac
        };

        p.clamp(THR_MIN, THR_MAX)
    }
}

// =============================
// Strain grouping + EM (set-only, based on winners sets)
// =============================

#[derive(Clone)]
struct StoredWinners {
    idx: usize,
    winners: Vec<String>, // max-tied genome names (within the same species)
}

// Group strains with identical read-support patterns. Two strains are
// equivalent when they occur in exactly the same read-level winner sets.
fn build_indistinguishable_groups(reads: &[StoredWinners]) -> Vec<Vec<String>> {
    // strain -> list of read positions where it appears (signature)
    let mut occ: HashMap<String, Vec<usize>> = HashMap::new();

    for (pos, r) in reads.iter().enumerate() {
        for s in &r.winners {
            occ.entry(s.clone()).or_insert_with(Vec::new).push(pos);
        }
    }

    // signature(Vec<usize>) -> strains
    let mut sig2strains: HashMap<Vec<usize>, Vec<String>> = HashMap::new();
    for (strain, mut v) in occ {
        v.sort_unstable();
        sig2strains.entry(v).or_insert_with(Vec::new).push(strain);
    }

    let mut groups: Vec<Vec<String>> = sig2strains.into_values().collect();
    for g in &mut groups {
        g.sort();
    }
    // Sort group members for deterministic output.
    groups.sort_by(|a, b| a[0].cmp(&b[0]));
    groups
}

// Each read identifies a candidate-group set G_r. Apply group-level EM with
// read-to-group weights proportional to the current group abundances.
fn em_over_groups(group_candidates_per_read: &[Vec<usize>], k: usize) -> (Vec<f64>, Vec<f64>) {
    if k == 0 || group_candidates_per_read.is_empty() {
        return (vec![], vec![]);
    }

    // Initialize group abundance with a weak prior from singleton-group reads;
    // use a uniform initialization when no such reads are available.
    let mut a = vec![1.0 / (k as f64); k];
    let mut anchor = vec![0usize; k];
    for cand in group_candidates_per_read {
        if cand.len() == 1 {
            anchor[cand[0]] += 1;
        }
    }
    let mut sum = 0.0;
    for j in 0..k {
        a[j] = (anchor[j] as f64) + 1e-6;
        sum += a[j];
    }
    if sum > 0.0 {
        for x in &mut a {
            *x /= sum;
        }
    }

    // EM
    let mut counts = vec![0.0f64; k];
    const EM_MAX_ITERS: usize = 200;
    const EM_EPS: f64 = 1e-10;

    for _ in 0..EM_MAX_ITERS {
        counts.fill(0.0);

        for cand in group_candidates_per_read {
            if cand.is_empty() {
                continue;
            }
            let mut denom = 0.0;
            for &g in cand {
                denom += a[g];
            }
            if denom <= 0.0 {
                continue;
            }
            for &g in cand {
                counts[g] += a[g] / denom;
            }
        }

        let sumc: f64 = counts.iter().sum();
        if sumc <= 0.0 {
            break;
        }

        let mut maxdiff = 0.0;
        for j in 0..k {
            let newa = counts[j] / sumc;
            let d = (newa - a[j]).abs();
            if d > maxdiff {
                maxdiff = d;
            }
            a[j] = newa;
        }
        if maxdiff < EM_EPS {
            break;
        }
    }

    // final counts with converged a
    counts.fill(0.0);
    for cand in group_candidates_per_read {
        if cand.is_empty() {
            continue;
        }
        let mut denom = 0.0;
        for &g in cand {
            denom += a[g];
        }
        if denom <= 0.0 {
            continue;
        }
        for &g in cand {
            counts[g] += a[g] / denom;
        }
    }

    (a, counts)
}

// =============================
// Global aggregation state
// =============================

struct GlobalAgg {
    is_paired: bool,

    // species -> (reads_count, indices)
    species: HashMap<String, (usize, Vec<usize>)>,

    // idx -> (reason, detail)
    discard: HashMap<usize, (&'static str, String)>,

    // reason -> count
    reason_stats: HashMap<&'static str, usize>,

    // In paired mode, retain the decision and winner set until its mate arrives.
    pair_buffer: HashMap<usize, (Decision, Vec<String>)>,

    // Buffer (query index, matches) pairs until the adaptive threshold is set.
    pending: Vec<(usize, Vec<(String, f64)>)>,

    sampler: RmaxSampler,
    threshold: Option<f64>,

    // Store only retained reads: species -> [(query index, winning strains)].
    kept_winners: HashMap<String, Vec<StoredWinners>>,
}

impl GlobalAgg {
    fn new(is_paired: bool) -> Self {
        Self {
            is_paired,
            species: HashMap::new(),
            discard: HashMap::new(),
            reason_stats: HashMap::new(),
            pair_buffer: HashMap::new(),
            pending: Vec::new(),
            sampler: RmaxSampler::new(),
            threshold: None,
            kept_winners: HashMap::new(),
        }
    }

    #[inline]
    fn inc_reason(&mut self, r: &'static str, n: usize) {
        *self.reason_stats.entry(r).or_insert(0) += n;
    }

    #[inline]
    fn record_kept(&mut self, sp: &str, idx: usize, winners: Vec<String>) {
        self.kept_winners
            .entry(sp.to_string())
            .or_insert_with(Vec::new)
            .push(StoredWinners { idx, winners });
    }

    fn apply_single(&mut self, idx: usize, dec: Decision, winners: Vec<String>) {
        match dec {
            Decision::Species(sp) => {
                let entry = self.species.entry(sp.clone()).or_insert((0, Vec::new()));
                entry.0 += 1;
                entry.1.push(idx);
                self.record_kept(&sp, idx, winners);
            }
            Decision::Discard(r, d) => {
                self.discard.insert(idx, (r, d));
                self.inc_reason(r, 1);
            }
        }
    }

    // In paired mode, count both mates only when they support the same species.
    fn merge_pair(&mut self, idx: usize, dec: Decision, winners: Vec<String>) {
        let mate = if idx % 2 == 0 { idx + 1 } else { idx - 1 };

        if let Some((other_dec, other_winners)) = self.pair_buffer.remove(&mate) {
            match (dec, other_dec) {
                (Decision::Species(a), Decision::Species(b)) => {
                    if a == b {
                        let entry = self.species.entry(a.clone()).or_insert((0, Vec::new()));
                        entry.0 += 2;

                        let (i1, i2, w1, w2) = if idx % 2 == 0 {
                            (idx, mate, winners, other_winners)
                        } else {
                            (mate, idx, other_winners, winners)
                        };

                        entry.1.push(i1);
                        entry.1.push(i2);

                        self.record_kept(&a, i1, w1);
                        self.record_kept(&a, i2, w2);
                    } else {
                        let msg = format!("{} vs {}", a, b);
                        self.discard.insert(idx, ("pair_conflict", msg.clone()));
                        self.discard.insert(mate, ("pair_conflict", msg));
                        self.inc_reason("pair_conflict", 2);
                    }
                }
                (Decision::Species(_), Decision::Discard(r, d))
                | (Decision::Discard(r, d), Decision::Species(_)) => {
                    self.discard
                        .insert(idx, ("pair_missing", format!("mate={}", mate)));
                    self.discard.insert(mate, (r, d));
                    self.inc_reason("pair_missing", 1);
                    self.inc_reason(r, 1);
                }
                (Decision::Discard(r1, d1), Decision::Discard(r2, d2)) => {
                    self.discard.insert(idx, (r1, d1));
                    self.discard.insert(mate, (r2, d2));
                    self.inc_reason(r1, 1);
                    self.inc_reason(r2, 1);
                }
            }
        } else {
            self.pair_buffer.insert(idx, (dec, winners));
        }
    }

    fn apply_decision(&mut self, idx: usize, dec: Decision, winners: Vec<String>) {
        if self.is_paired {
            self.merge_pair(idx, dec, winners);
        } else {
            self.apply_single(idx, dec, winners);
        }
    }

    fn finalize_orphans(&mut self) {
        if !self.is_paired {
            self.pair_buffer.clear();
            return;
        }

        let keys: Vec<usize> = self.pair_buffer.keys().cloned().collect();
        for i in keys {
            self.discard
                .insert(i, ("pair_missing", "file_end_or_missing_mate".to_string()));
            self.inc_reason("pair_missing", 1);
        }
        self.pair_buffer.clear();
    }
}

// =============================
// Main entry point
// =============================

pub fn colored_query_output<MH: HashFunctionFactory, CX: ColorsManager>(
    colormap: &<CX::ColorsMergeManagerType as ColorsMergeManager>::GlobalColorsTableReader,
    mut colored_query_buckets: Vec<SingleBucket>,
    output_file: PathBuf,
    temp_dir: PathBuf,
    query_kmers_count: &[u64],
    colored_query_output_format: ColoredQueryOutputFormat,
    emit_jsonl: bool,
    is_paired: bool,
) -> anyhow::Result<()> {
    PHASES_TIMES_MONITOR
        .write()
        .start_phase("phase: colored query output".to_string());

    let buckets_count = colored_query_buckets.len();

    let max_bucket_queries_count = (((query_kmers_count.len() + 1) as u64)
        .nq_div_ceil(QUERIES_COUNT_MIN_BATCH)
        * QUERIES_COUNT_MIN_BATCH) as usize;

    static OPS_COUNT: AtomicUsize = AtomicUsize::new(0);
    static COL_COUNT: AtomicUsize = AtomicUsize::new(0);

    colored_query_buckets.reverse();
    let buckets_channel = Mutex::new(colored_query_buckets);

    // ===== Optional JSONL output =====
    let (maybe_query_output, output_sync_condvar, output_path_final): (
        Option<Mutex<(BufWriter<QueryOutputFileWriter>, usize)>>,
        Option<Condvar>,
        Option<PathBuf>,
    ) = if emit_jsonl {
        // Use output_path to distinguish the resolved path from output_file.
        let output_path: PathBuf = if output_file.extension().is_none() {
            output_file.with_extension("jsonl")
        } else {
            output_file.clone()
        };

        let query_output_file = File::create(&output_path).log_unrecoverable_error_with_data(
            "Cannot create output file",
            output_path.display(),
        )?;

        let query_output: Mutex<(BufWriter<QueryOutputFileWriter>, usize)> = Mutex::new((
            BufWriter::new(match output_path.extension().and_then(|e| e.to_str()) {
                Some("lz4") => QueryOutputFileWriter::LZ4Compressed(
                    lz4::EncoderBuilder::new()
                        .level(4)
                        .build(query_output_file)
                        .unwrap(),
                ),
                Some("gz") => QueryOutputFileWriter::GzipCompressed(
                    flate2::GzBuilder::new().write(query_output_file, Compression::default()),
                ),
                _ => QueryOutputFileWriter::Plain(query_output_file),
            }),
            0usize,
        ));

        (Some(query_output), Some(Condvar::new()), Some(output_path))
    } else {
        (None, None, None)
    };

    // ===== Global aggregation =====
    let global_agg = Mutex::new(GlobalAgg::new(is_paired));

    (0..rayon::current_num_threads()).into_par_iter().for_each(|_| {
        #[derive(Copy, Clone)]
        struct QueryColorListItem {
            color: ColorIndexType,
            count: u64,
            next_index: usize,
        }

        let mut queries_colors_list_pool = vec![];
        let mut queries_results = vec![(0u32, 0usize); max_bucket_queries_count];
        let mut temp_colors_list = vec![];
        let mut epoch = 0u32;

        let mut local_rmax: Vec<f64> = Vec::new();
        let mut local_pending: Vec<(usize, Vec<(String, f64)>)> = Vec::new();

        while let Some(input) = {
            let mut lock = buckets_channel.lock();
            lock.pop()
        } {
            epoch = epoch.wrapping_add(1);
            queries_colors_list_pool.clear();
            local_rmax.clear();
            local_pending.clear();

            let start_query_index = input.index as usize * max_bucket_queries_count / buckets_count;

            CompressedBinaryReader::new(
                &input.path,
                RemoveFileMode::Remove {
                    remove_fs: !KEEP_FILES.load(Ordering::Relaxed),
                },
                DEFAULT_PREFETCH_AMOUNT,
            )
            .decode_all_bucket_items::<QueryColoredCountersSerializer, _>(
                (Vec::new(), Vec::new()),
                &mut (),
                |counters, _| {
                    for query in counters.queries {
                        let (entry_epoch, colors_map_index) = &mut queries_results
                            [query.query_index as usize - start_query_index - 1];

                        if *entry_epoch != epoch {
                            *entry_epoch = epoch;
                            *colors_map_index = usize::MAX;
                        }

                        assert_eq!(counters.colors.len() % 2, 0);

                        for range in counters.colors.chunks(2) {
                            let ColorsRange::Range(range) = ColorsRange::from_slice(range);

                            OPS_COUNT.fetch_add(1, Ordering::Relaxed);
                            COL_COUNT.fetch_add(range.len(), Ordering::Relaxed);

                            for color in range {
                                queries_colors_list_pool.push(QueryColorListItem {
                                    color,
                                    count: query.count,
                                    next_index: *colors_map_index,
                                });
                                *colors_map_index = queries_colors_list_pool.len() - 1;
                            }
                        }
                    }
                },
            );

            let bucket_index = input.index;

            // Write the intermediate JSONL stream by bucket.
            let mut maybe_compressed_stream: Option<CompressedBinaryWriter> = if emit_jsonl {
                Some(CompressedBinaryWriter::new(
                    &temp_dir.join("query-data"),
                    &(
                        get_memory_mode(SwapPriority::ColoredQueryBuckets),
                        CompressedBinaryWriter::CHECKPOINT_SIZE_UNLIMITED,
                        get_compression_level_info(),
                    ),
                    bucket_index as usize,
                    &(),
                ))
            } else {
                None
            };

            let mut jsonline_buffer = vec![];

            // Process the queries in the current bucket.
            for (query, mut query_colors_list_index) in
                queries_results.iter().enumerate().filter_map(|(i, r)| {
                    if r.0 != epoch {
                        None
                    } else {
                        Some((i + start_query_index, r.1))
                    }
                })
            {
                // 1) flatten (color, count)
                temp_colors_list.clear();
                while query_colors_list_index != usize::MAX {
                    let el = &queries_colors_list_pool[query_colors_list_index];
                    temp_colors_list.push((el.color, el.count));
                    query_colors_list_index = el.next_index;
                }
                temp_colors_list.sort_unstable_by_key(|r| r.0);

                // 2) Aggregate identical colors into support fractions and build matches.
                let mut matches: Vec<(String, f64)> = Vec::new();
                for grp in temp_colors_list.nq_group_by(|a, b| a.0 == b.0) {
                    let color_index = grp[0].0;
                    let color_presence = grp.iter().map(|x| x.1).sum::<u64>();
                    let denom = query_kmers_count[query as usize].max(1) as f64;
                    let frac = (color_presence as f64) / denom;

                    let key = match colored_query_output_format {
                        ColoredQueryOutputFormat::JsonLinesWithNumbers => color_index.to_string(),
                        ColoredQueryOutputFormat::JsonLinesWithNames => {
                            colormap.get_color_name(color_index, true).to_string()
                        }
                    };

                    matches.push((key, frac));
                }

                // 3) Collect maximum support and pending queries for thresholding.
                let mut rmax = 0.0f64;
                for &(_, v) in &matches {
                    if v > rmax {
                        rmax = v;
                    }
                }
                local_rmax.push(rmax);
                local_pending.push((query, matches.clone()));

                // 4) Optionally write all matches to JSONL without further filtering.
                if let Some(ref mut stream) = maybe_compressed_stream {
                    if emit_jsonl {
                        jsonline_buffer.clear();
                        write!(jsonline_buffer, "{{\"query_index\":{}, \"matches\":{{", query).unwrap();
                        for (i, (k, v)) in matches.iter().enumerate() {
                            if i != 0 {
                                write!(jsonline_buffer, ",").unwrap();
                            }
                            write!(jsonline_buffer, "\"{}\": {:.2}", k, v).unwrap();
                        }
                        writeln!(jsonline_buffer, "}}}}").unwrap();
                        stream.write_data(&jsonline_buffer[..]);
                    }
                }
            }

            // 5) Merge threshold samples. Once FIRST_N observations are available,
            // set the threshold and make decisions bucket by bucket.
            {
                let mut g = global_agg.lock();

                if g.threshold.is_none() {
                    for r in &local_rmax {
                        g.sampler.try_push(*r);
                    }
                    g.pending.extend(local_pending.drain(..));

                    if g.sampler.len() >= FIRST_N {
                        let t = g.sampler.compute_threshold();
                        g.threshold = Some(t);

                        let pending = std::mem::take(&mut g.pending);
                        for (idx, matches) in pending {
                            let (dec, winners) = decide_species_and_winners(&matches, t);
                            g.apply_decision(idx, dec, winners);
                        }
                    }
                } else {
                    let t = g.threshold.unwrap();
                    for (idx, matches) in local_pending.drain(..) {
                        let (dec, winners) = decide_species_and_winners(&matches, t);
                        g.apply_decision(idx, dec, winners);
                    }
                }
            }

            // 6) Append this bucket's JSONL stream to the final file in bucket order.
            if let (Some(stream), Some(ref condvar), Some(ref query_output)) =
                (maybe_compressed_stream, output_sync_condvar.as_ref(), maybe_query_output.as_ref())
            {
                let stream_path = stream.get_path();
                stream.finalize();

                let mut decompress_stream = CompressedBinaryReader::new(
                    stream_path,
                    RemoveFileMode::Remove { remove_fs: true },
                    DEFAULT_PREFETCH_AMOUNT,
                );

                let mut guard = query_output.lock();
                while guard.1 != bucket_index as usize {
                    condvar.wait(&mut guard);
                }

                let (queries_file, query_write_index) = guard.deref_mut();
                std::io::copy(&mut decompress_stream.get_single_stream(), queries_file).unwrap();
                *query_write_index += 1;
                condvar.notify_all();
            }
        }
    });

    // ===== Finalization: if fewer than FIRST_N reads were observed, set the
    // threshold from all available observations and resolve pending queries. =====
    let mut g = global_agg.lock();

    let t = if let Some(t) = g.threshold {
        t
    } else {
        let t = g.sampler.compute_threshold();
        g.threshold = Some(t);

        let pending = std::mem::take(&mut g.pending);
        for (idx, matches) in pending {
            let (dec, winners) = decide_species_and_winners(&matches, t);
            g.apply_decision(idx, dec, winners);
        }
        t
    };

    // Resolve reads whose mates were not observed in paired mode.
    g.finalize_orphans();

    // Use the JSONL prefix as the base name when enabled; otherwise use output_file.
    let base = if let Some(p) = output_path_final.as_ref() {
        p.with_extension("")
    } else {
        output_file.with_extension("")
    };

    let species_out = base.with_extension("species_counts.tsv");
    let discard_out = base.with_file_name(format!(
        "{}_discard_log.tsv",
        base.file_name().unwrap().to_string_lossy()
    ));
    let summary_out = base.with_file_name(format!(
        "{}_summary.txt",
        base.file_name().unwrap().to_string_lossy()
    ));

    // Write strain-group output.
    let strain_group_map_out = base.with_extension("strain_group_map.tsv");
    let strain_group_abund_out = base.with_extension("strain_group_abundance.tsv");

    // Report distinguishable strains individually and unresolved strains as
    // semicolon-delimited groups.
    let distinguishable_out = base.with_extension("distinguishable_strains.tsv");
    let read_candidates_out = base.with_extension("read_group_candidates.tsv");


    // 1) species_counts.tsv
    {
        let mut v: Vec<_> = g
            .species
            .iter()
            .map(|(sp, (cnt, idxs))| (sp.clone(), *cnt, idxs.clone()))
            .collect();
        v.sort_by_key(|t| Reverse(t.1));

        let mut w = File::create(&species_out)?;
        for (sp, cnt, idxs) in v {
            writeln!(
                w,
                "{}\t{}\t{}",
                sp,
                cnt,
                idxs.iter()
                    .map(|x| x.to_string())
                    .collect::<Vec<_>>()
                    .join(",")
            )?;
        }
    }

    // 2) discard_log.tsv
    {
        let mut v: Vec<_> = g
            .discard
            .iter()
            .map(|(i, (r, d))| (*i, *r, d.clone()))
            .collect();
        v.sort_by_key(|t| t.0);

        let mut w = File::create(&discard_out)?;
        writeln!(w, "index\treason\tdetail")?;
        for (i, r, d) in v {
            writeln!(w, "{}\t{}\t{}", i, r, d)?;
        }
    }

    // 3) summary.txt
    {
        let kept: usize = g.species.values().map(|(n, _)| *n).sum();
        let mut w = File::create(&summary_out)?;
        writeln!(
            w,
            "threshold(T)={:.6} clipped_to=[{:.6},{:.6}]",
            t, THR_MIN, THR_MAX
        )?;
        writeln!(w, "discard_reasons:")?;
        let mut rs: Vec<_> = g.reason_stats.iter().map(|(k, v)| (*k, *v)).collect();
        rs.sort_by(|a, b| a.0.cmp(&b.0));
        for (k, v) in rs {
            writeln!(w, "  {}: {}", k, v)?;
        }
        writeln!(w, "kept_species_reads={}", kept)?;
    }

    // 4) Build strain groups, run group-level EM and write distinguishable_strains.tsv.
    {
        let mut wmap = File::create(&strain_group_map_out)?;
        let mut wab = File::create(&strain_group_abund_out)?;
        let mut wdist = File::create(&distinguishable_out)?;
        let mut wrcand = File::create(&read_candidates_out)?;

        writeln!(wmap, "species\tgroup_id\tstrains")?;
        writeln!(wab, "species\tgroup_id\tassigned_reads\trel_abundance_within_species\tstrains")?;
        writeln!(wdist, "species\tstrain_id\tassigned_reads\trel_abundance_within_species")?;

        // Candidate groups per read, restricted to reads retained in the final set.
        writeln!(wrcand, "species\tread_idx\tcandidate_groups\twinner_strains")?;

        let mut species_list: Vec<String> = g.kept_winners.keys().cloned().collect();
        species_list.sort();

        for sp in species_list {
            let reads = match g.kept_winners.get(&sp) {
                Some(v) => v,
                None => continue,
            };
            if reads.is_empty() { continue; }

            // 4.1 Group strains with identical retained-read support patterns.
            let groups = build_indistinguishable_groups(reads);
            let k = groups.len();
            if k == 0 { continue; }

            // strain -> group index
            let mut strain2g: HashMap<String, usize> = HashMap::new();
            for (gi, strains) in groups.iter().enumerate() {
                for s in strains {
                    strain2g.insert(s.clone(), gi);
                }
            }

            // Write the strain-to-group mapping.
            for (gi, strains) in groups.iter().enumerate() {
                let group_id = format!("G{}", gi + 1);
                writeln!(wmap, "{}\t{}\t{}", sp, group_id, strains.join(","))?;
            }

            // 4.2 Map each read's winner set to its candidate-group set.
            let mut cand_per_read: Vec<Vec<usize>> = Vec::with_capacity(reads.len());
            for r in reads {
                let mut c: Vec<usize> = r
                    .winners
                    .iter()
                    .filter_map(|s| strain2g.get(s).cloned())
                    .collect();
                c.sort_unstable();
                c.dedup();

                // Write the candidate-group list for each retained read.
                let cand_str = c
                    .iter()
                    .map(|gi| format!("G{}", gi + 1))
                    .collect::<Vec<_>>()
                    .join(",");

                let winners_str = r.winners.join(",");

                writeln!(wrcand, "{}\t{}\t{}\t{}", sp, r.idx, cand_str, winners_str)?;

                cand_per_read.push(c);
            }

            // 4.3 Estimate relative support at the strain-group level by EM.
            let (_a, counts) = em_over_groups(&cand_per_read, k);
            let sumc: f64 = counts.iter().sum();
            if sumc <= 0.0 { continue; }

            let mut out: Vec<(usize, f64, f64)> = Vec::new();
            for gi in 0..k {
                let assigned = counts[gi];
                let rel = assigned / sumc;
                out.push((gi, assigned, rel));
            }
            out.sort_by(|a, b| b.1.partial_cmp(&a.1).unwrap_or(std::cmp::Ordering::Equal));

            for (gi, assigned, rel) in out {
                let group_id = format!("G{}", gi + 1);
                let strains = &groups[gi];

                writeln!(
                    wab,
                    "{}\t{}\t{:.6}\t{:.8}\t{}",
                    sp, group_id, assigned, rel, strains.join(",")
                )?;

                // Emit singleton groups as individual strains and unresolved groups
                // as semicolon-delimited strain identifiers.
                let mut ids: Vec<String> = strains.iter().map(|x| strain_id_of(x).to_string()).collect();
                ids.sort();
                ids.dedup();
                let strain_id_joined = ids.join(";");

                writeln!(
                    wdist,
                    "{}\t{}\t{:.6}\t{:.8}",
                    sp, strain_id_joined, assigned, rel
                )?;
            }
        }
    }


    ggcat_logging::info!(
        "Operations count: {} vs real {}",
        OPS_COUNT.load(Ordering::Relaxed),
        COL_COUNT.load(Ordering::Relaxed)
    );

    Ok(())
}
