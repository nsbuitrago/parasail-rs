//! Regression test: `get_traceback_strings` must not free memory owned by the
//! parasail C library with Rust's allocator.
//!
//! This test binary installs a global allocator that tags every block it
//! hands out; freeing a block it did not allocate (e.g. a string malloc'd by
//! parasail) aborts the process.

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
fn traceback_strings_are_freed_by_parasail() {
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
