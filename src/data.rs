use crate::builtins::*;

impl<T: Rank1> Matrix<T> {
    /// Creates a new empty `Matrix` type given
    /// # rows as `m` and # cols as `n` and storage type as `storage`
    pub fn new(m: usize, n: usize, storage: Option<Storage>) -> Self {
        let mut data: Vec<T> = Vec::new();
        data.reserve_exact(m * n);

        Self {
            data,
            n_rows: m,
            n_cols: n,
            storage_type: storage, // Default to column major
        }
    }
}

impl<T: Rank1> FromRowVector<T> for Matrix<T> {
    fn as_row_major(&mut self, data: &[T]) {
        self.data = Vec::from(data);
        self.storage_type = Some(Storage::R); // Set storage type to row major
    }

    fn as_column_major(&mut self, data: &[T]) {
        self.data = Vec::new(); // Allocate a new vector
        self.data.clear();
        self.data.reserve_exact(data.len());

        // m is lda and total size is [lda x n]
        let (m, n) = (self.n_rows, self.n_cols);

        for i in 0..m {
            for j in 0..n {
                let cm_idx = i + j * n; // column major index
                let rm_idx = j + i * m; // row major index
                self.data[cm_idx] = data[rm_idx];
            }
        }

        self.storage_type = Some(Storage::C); // Set storage type to column major
    }

    fn with_data(&mut self, data: &[T]) {
        match self.storage_type {
            Some(Storage::R) => self.as_row_major(data),
            _ => self.as_column_major(data),
        }
    }
}

impl<T: Rank1> BandMatrix<T> {
    fn new(m: usize, n: usize, storage: Option<Storage>, kl: usize, ku: usize, uplo: Uplo) -> Self {
        let matrix = Matrix::new(m, n, storage); // Make base matrix.
        Self {
            matrix,
            ku,
            kl,
            uplo,
        } // Return a BandMatrix with all metadata.
    }

    fn uplo_diags(kl: usize, ku: usize, uplo: Uplo) -> (usize, usize) {
        let kl = match uplo {
            Uplo::L => kl,
            Uplo::U => 0,
        };

        let ku = match uplo {
            Uplo::L => 0,
            Uplo::U => ku,
        };

        (kl, ku)
    }
}

impl<T: Rank1> FromRowVector<T> for BandMatrix<T> {
    fn as_column_major(&mut self, data: &[T]) {
        self.matrix.data.clear();

        let (kl, ku) = Self::uplo_diags(self.kl, self.ku, self.uplo);
        // m is lda and total size is [lda x n]
        let (m, n) = (self.matrix.n_rows, self.matrix.n_cols);

        for j in 0..n {
            let k = ku - j;
            for i in std::cmp::max(0, j - ku)..std::cmp::min(m, j + kl + 1) {
                self.matrix.data[(k + i) + j * m] = data[i + j * m];
            }
        }

        self.matrix.storage_type = Some(Storage::C);
    }

    fn as_row_major(&mut self, data: &[T]) {
        self.matrix.data.clear();

        let (kl, ku) = Self::uplo_diags(self.kl, self.ku, self.uplo);
        // m is lda and total size is [lda x n]
        let (m, n) = (self.matrix.n_rows, self.matrix.n_cols);

        for i in 0..m {
            let k = ku - i;
            for j in std::cmp::max(0, i - kl)..std::cmp::min(m, i + ku + 1) {
                self.matrix.data[(k + j) + i * n] = data[i + j * m];
            }
        }

        self.matrix.storage_type = Some(Storage::R);
    }

    fn with_data(&mut self, data: &[T]) {
        match self.matrix.storage_type {
            Some(Storage::R) => self.as_row_major(data),
            _ => self.as_column_major(data),
        }
    }
}

impl<T: Rank1> Vector<T> {
    fn new(n: usize, incx: i32) -> Self {
        let mut data = Vec::new();
        data.reserve_exact(1 + (n - 1) * incx.wrapping_abs() as usize);

        Self { data, n, incx }
    }
}

impl<T: Rank1> FromRowVector<T> for Vector<T> {
    fn as_column_major(&mut self, data: &[T]) {
        let _ = data;
        unimplemented!();
    }

    fn as_row_major(&mut self, data: &[T]) {
        let _ = data;
        unimplemented!();
    }

    fn with_data(&mut self, data: &[T]) {
        let data_size: usize = if self.incx < 0 {
            self.data.len() - 1
        } else {
            0
        };
        let incx = self.incx.wrapping_abs() as usize;

        for (i, val) in data.iter().enumerate() {
            let idx = data_size + i * incx;
            self.data[idx] = *val;
        }
    }
}
