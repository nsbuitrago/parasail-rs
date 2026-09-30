//! Regression tests for FFI allocation ownership.
//!
//! Parasail allocates traceback and decoded CIGAR strings with its C
//! allocator. Rust must copy those strings, then release their original
//! allocations with the matching C cleanup path rather than Rust's allocator.
//!
//! This test binary installs a global allocator that tags every block it hands
//! out. It aborts if Rust attempts to free a foreign allocation, such as a
//! string returned by Parasail.

use parasail_rs::prelude::*;
use std::alloc::{GlobalAlloc, Layout, System};

struct Tagging;

const TAG: u64 = 0x7061_7261_7361_696c; // "parasail"
const HEADER: usize = 16;

unsafe impl GlobalAlloc for Tagging {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let layout =
            Layout::from_size_align(layout.size() + HEADER, layout.align().max(HEADER)).unwrap();
        let p = System.alloc(layout);
        if p.is_null() {
            return p;
        }
        (p as *mut u64).write(TAG);
        p.add(HEADER)
    }

    unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
        let base = ptr.sub(HEADER);
        if (base as *const u64).read() != TAG {
            // not allocated by the Rust allocator
            std::process::abort();
        }
        (base as *mut u64).write(0);
        let layout =
            Layout::from_size_align(layout.size() + HEADER, layout.align().max(HEADER)).unwrap();
        System.dealloc(base, layout);
    }
}

#[global_allocator]
static GLOBAL: Tagging = Tagging;

#[test]
fn traceback_strings_do_not_use_rust_allocator() {
    let matrix = Matrix::create(b"ACGT", 2, -1).unwrap();
    let aligner = Aligner::new()
        .matrix(matrix)
        .gap_open(5)
        .gap_extend(2)
        .global()
        .use_trace()
        .scan()
        .solution_width(16)
        .build();
    let (query, reference) = (b"ACGTACGTACGT", b"ACGTTCGTACGA");
    let result = aligner.align(Some(query), reference).unwrap();
    for _ in 0..3 {
        let tb = result.get_traceback_strings(query, reference).unwrap();
        assert_eq!(tb.query, "ACGTACGTACGT");
        assert_eq!(tb.reference, "ACGTTCGTACGA");
    }
}

#[test]
fn decoded_cigar_does_not_use_rust_allocator() {
    let matrix = Matrix::default();
    let aligner = Aligner::new().matrix(matrix).use_trace().build();
    let (query, reference) = (b"ACGT", b"ACGT");
    let result = aligner.align(Some(query), reference).unwrap();

    for _ in 0..3 {
        let cigar = result.get_cigar(query, reference).unwrap();
        assert_eq!(cigar, "4=");
    }
}
