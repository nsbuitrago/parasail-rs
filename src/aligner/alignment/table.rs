use std::fmt;

/// A view into a score table from an alignment result.
///
/// # Layout
/// - Rows represent positions in the query sequence (0..query_len)
/// - Columns represent positions in the reference sequence (0..ref_len)
/// - Data is stored as a flattened 1D array: `data[row * cols + col]`
///
/// # Example
/// ```rust,no_run
/// use parasail_rs::Aligner;
///
/// let query = b"ACGT";
/// let reference = b"ACGT";
/// let aligner = Aligner::new().use_table().build();
/// let result = aligner.align(Some(query), reference)?;
///
/// let table = result.get_score_table()?;
/// println!("Table dimensions: {} rows x {} cols", table.rows(), table.cols());
///
/// // Access specific cell
/// if let Some(score) = table.get(0, 0) {
///     println!("Score at (0, 0): {}", score);
/// }
///
/// // Get final alignment score
/// println!("Final score: {}", table.last());
/// # Ok::<(), Box<dyn std::error::Error>>(())
/// ```
#[derive(Debug)]
pub struct Table<'a> {
    inner: &'a [i32],
    rows: usize,
    cols: usize,
}

impl<'a> Table<'a> {
    /// Create a new Table view from a flat 1D array.
    ///
    /// # Panics
    /// Panics in debug mode if data.len() != rows * cols
    pub(crate) fn new(data: &'a [i32], rows: usize, cols: usize) -> Self {
        debug_assert_eq!(
            data.len(),
            rows * cols,
            "Table size mismatch: expected {} elements ({}x{}), got {}",
            rows * cols,
            rows,
            cols,
            data.len()
        );
        Self {
            inner: data,
            rows,
            cols,
        }
    }

    /// Get the value at the given row and column index.
    ///
    /// Returns `None` if indices are out of bounds.
    ///
    /// # Example
    /// ```rust,no_run
    /// # use parasail_rs::Aligner;
    /// # let query = b"ACGT";
    /// # let reference = b"ACGT";
    /// # let aligner = Aligner::new().use_table().build();
    /// # let result = aligner.align(Some(query), reference)?;
    /// let table = result.get_score_table()?;
    /// if let Some(score) = table.get(2, 3) {
    ///     println!("Score at (2, 3): {}", score);
    /// }
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn get(&self, row: usize, col: usize) -> Option<i32> {
        if row < self.rows && col < self.cols {
            Some(self.inner[row * self.cols + col])
        } else {
            None
        }
    }

    /// Get the number of rows (query length).
    pub fn rows(&self) -> usize {
        self.rows
    }

    /// Get the number of columns (reference length).
    pub fn cols(&self) -> usize {
        self.cols
    }

    /// Get the raw underlying data as a 1D slice.
    ///
    /// The data is stored in row-major order, so element at (row, col)
    /// can be accessed at index `row * cols + col`.
    pub fn as_slice(&self) -> &[i32] {
        self.inner
    }

    /// Get the value at the last cell (bottom-right of DP table).
    pub fn last(&self) -> i32 {
        self.inner[self.inner.len() - 1]
    }
}

impl<'a> fmt::Display for Table<'a> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Table ({}x{}):", self.rows, self.cols)?;
        for row in 0..self.rows {
            write!(f, "[")?;
            for col in 0..self.cols {
                if col > 0 {
                    write!(f, ", ")?;
                }
                write!(f, "{}", self.get(row, col).unwrap())?;
            }
            writeln!(f, "]")?;
        }
        Ok(())
    }
}
