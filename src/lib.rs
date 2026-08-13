#![cfg_attr(target_arch = "wasm32", feature(stdarch_wasm_atomic_wait))]

#[cfg(target_arch = "wasm32")]
use skani::chain;
#[cfg(target_arch = "wasm32")]
use skani::screen;
#[cfg(target_arch = "wasm32")]
use skani::seeding;
#[cfg(target_arch = "wasm32")]
use skani::regression;
#[cfg(target_arch = "wasm32")]
use skani::params;
#[cfg(target_arch = "wasm32")]
use skani::types;

// WASM entry point (TESTING)
#[cfg(target_arch = "wasm32")]
use wasm_bindgen::prelude::*;
#[cfg(target_arch = "wasm32")]
use std::io::Cursor;
#[cfg(target_arch = "wasm32")]
use needletail::parse_fastx_reader;
#[cfg(target_arch = "wasm32")]
use rayon::prelude::*;
#[cfg(target_arch = "wasm32")]
use std::sync::Mutex;
#[cfg(target_arch = "wasm32")]

#[cfg(target_arch = "wasm32")]
pub use wasm_bindgen_rayon::init_thread_pool;

#[cfg(target_arch = "wasm32")]
fn fasta_string_to_sketch(fasta: &str, sketch_params: &params::SketchParams) -> Option<types::Sketch> {
    let mut new_sketch = types::Sketch::new(
        sketch_params.marker_c,
        sketch_params.c,
        sketch_params.k,
        "wasm_input".to_string(),
        sketch_params.use_aa,
    );
    let cursor = Cursor::new(fasta.as_bytes().to_vec());
    let reader = parse_fastx_reader(cursor);
    if reader.is_err() {
        return None;
    }
    let mut reader = reader.unwrap();
    let mut j = 0;
    while let Some(record) = reader.next() {
        if let Ok(record) = record {
            let seq = record.seq();
            if seq.len() >= params::MIN_LENGTH_CONTIG {
                new_sketch.contigs.push(String::from_utf8(record.id().to_vec()).unwrap());
                new_sketch.contig_lengths.push(seq.len() as types::GnPosition);
                new_sketch.total_sequence_length += seq.len();
                seeding::fmh_seeds(&seq, sketch_params, j as u32, &mut new_sketch, true);
                j += 1;
            }
        }
    }
    if j > 0 { Some(new_sketch) } else { None }
}




#[cfg(target_arch = "wasm32")]
use web_sys::console;

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
pub fn compare_genomes(fasta_a: &str, fasta_b: &str, k: u32, c: u32, marker_c: u32) -> String {
    console_error_panic_hook::set_once();
    
    let sketch_params = params::SketchParams::new(
    marker_c as usize,
    c as usize,
    k as usize,
    false,
    false,
);
    let command_params = params::CommandParams {
        screen: false,
        screen_val: 0.0,
        mode: params::Mode::Dist,
        out_file_name: "".to_string(),
        ref_files: vec![],
        query_files: vec![],
        refs_are_sketch: false,
        queries_are_sketch: false,
        robust: false,
        median: false,
        sparse: false,
        full_matrix: false,
        diagonal: false,
        max_results: 10000,
        individual_contig_q: false,
        individual_contig_r: false,
        min_aligned_frac: 0.0,
        both_min_aligned_frac: 0.0,
        keep_refs: false,
        est_ci: false,
        learned_ani: false,
        detailed_out: false,
        distance: false,
        rescue_small: false,
        separate_sketches: false,
        short_header: false,
    };

    console::log_1(&"Sketching A...".into());
    let ref_sketch = fasta_string_to_sketch(fasta_a, &sketch_params);
    console::log_1(&"Sketching B...".into());
    let query_sketch = fasta_string_to_sketch(fasta_b, &sketch_params);

    if ref_sketch.is_none() || query_sketch.is_none() {
        return "Error: invalid FASTA input".to_string();
    }

    let ref_sketch = ref_sketch.unwrap();
    let query_sketch = query_sketch.unwrap();

    console::log_1(&format!("ref markers: {}, query markers: {}", ref_sketch.marker_seeds.len(), query_sketch.marker_seeds.len()).into());
    let passed_screen = screen::check_markers_quickly(&ref_sketch, &query_sketch, 0.0, true);   
    console::log_1(&format!("passed_screen: {}", passed_screen).into());

let model_opt = regression::get_model(sketch_params.c, true);
let map_params = chain::map_params_from_sketch(&ref_sketch, false, &command_params, &model_opt);   
 let ani_res = chain::chain_seeds(&ref_sketch, &query_sketch, map_params);
    console::log_1(&"chain_seeds done...".into());

    format!("{{\"ani\": {}, \"af_ref\": {}, \"af_query\": {}}}", ani_res.ani, ani_res.align_fraction_query, ani_res.align_fraction_ref)
}
#[cfg(target_arch = "wasm32")]
fn default_command_params() -> params::CommandParams {
    params::CommandParams {
        screen: false,
        screen_val: 0.0,
        mode: params::Mode::Dist,
        out_file_name: "".to_string(),
        ref_files: vec![],
        query_files: vec![],
        refs_are_sketch: false,
        queries_are_sketch: false,
        robust: false,
        median: false,
        sparse: false,
        full_matrix: false,
        diagonal: false,
        max_results: 10000,
        individual_contig_q: false,
        individual_contig_r: false,
        min_aligned_frac: 0.0,
        both_min_aligned_frac: 0.0,
        keep_refs: false,
        est_ci: false,
        learned_ani: false,
        detailed_out: false,
        distance: false,
        rescue_small: false,
        separate_sketches: false,
        short_header: false,
    }
}
#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
pub fn compare_triangle(
    fastas: Vec<String>,
    names: Vec<String>,
    k: u32,
    c: u32,
    marker_c: u32,
) -> String {
    rayon::spawn(|| {}); // force pool init
let threads = rayon::current_num_threads();
console::log_1(&format!("rayon threads: {}", threads).into());
    console_error_panic_hook::set_once();
 
    if fastas.len() < 2 {
        return "Error: need at least 2 genomes for triangle comparison".to_string();
    }
 
    let sketch_params = params::SketchParams::new(
        marker_c as usize,
        c as usize,
        k as usize,
        false,
        false,
    );
    let command_params = default_command_params();
 
    // Sketch all genomes in parallel:
    console::log_1(&format!("Sketching {} genomes in parallel...", fastas.len()).into());
    let sketches: Vec<Option<types::Sketch>> = fastas
        .par_iter()
        .enumerate()
        .map(|(i, fasta)| {
            let mut s = fasta_string_to_sketch(fasta, &sketch_params)?;
            // use the provided filename as the sketch name
            s.file_name = names.get(i).cloned().unwrap_or_else(|| format!("genome_{}", i));
            Some(s)
        })
        .collect();
 
    // filter out any files that failed to sketch
    let sketches: Vec<types::Sketch> = sketches.into_iter().flatten().collect();
    let n = sketches.len();
 
    if n < 2 {
        return "Error: fewer than 2 genomes could be sketched. Check file format and minimum contig length (500bp).".to_string();
    }
 
    console::log_1(&format!("{} genomes sketched successfully.", n).into());
 
    // Build marker index for fast pre-screening:
    // Making sure we only do chain seeds on sketches that are alike 
    let kmer_to_sketch = screen::kmer_to_sketch_from_refs(&sketches);
 
    // Compare all pairs in parallel (upper triangle only): 
    // We use the same pattern as Skani's own triangle.rs:
    //   - outer par_iter over i (each genome as reference)
    //   - screen_refs finds only the j genomes worth comparing
    //   - if j > i, compare and store result (avoids duplicate pairs)
    //   - Mutex protects the shared results Vec 
    let results: Mutex<Vec<serde_json::Value>> = Mutex::new(vec![]);
    let model_opt = regression::get_model(sketch_params.c, true);
    
    (0..n - 1)
        .collect::<Vec<usize>>()
        .into_par_iter()
        .for_each(|i| {
            let ref_sketch = &sketches[i];
 
            // get only the genomes that pass the marker screen
            let screened = screen::screen_refs(
                0.0,
                &kmer_to_sketch,
                ref_sketch,
                &sketch_params,
                &sketches,
                true,
            );

            

            screened.into_par_iter().for_each(|j| {
                if j <= i {
                    return; // only upper triangle
                }
                let query_sketch = &sketches[j];
                let map_params = chain::map_params_from_sketch(
                    ref_sketch,
                    false,
                    &command_params,
                    &model_opt,
                );
                let ani_res = chain::chain_seeds(ref_sketch, query_sketch, map_params);
                if ani_res.ani > 0.0 {
                    let mut locked = results.lock().unwrap();
                    locked.push(serde_json::json!({
                        "i": i,
                        "j": j,
                        "name_i": ref_sketch.file_name,
                        "name_j": query_sketch.file_name,
                        "ani": ani_res.ani,
                        "af_query": ani_res.align_fraction_query,
                        "af_ref": ani_res.align_fraction_ref,
                    }));
                }
            });
        });
 
    let results = results.into_inner().unwrap();
    console::log_1(&format!("Triangle done — {} pairs compared.", results.len()).into());
 
    // formatting for matrix 
    let response = serde_json::json!({
        "n": n,
        "names": names,
        "results": results,
    });
 
    serde_json::to_string(&response).unwrap()
}