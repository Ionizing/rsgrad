// Minimal .npy and .npz writer for numpy interoperability.
//
// Implements:
//   - write_npy_c128(path, data, shape)  → complex128 .npy
//   - write_npy_f64(path, data, shape)   → float64 .npy
//   - write_npy_i64(path, data, shape)   → int64 .npy
//   - NpzWriter: write multiple arrays to a .npz (zip) archive
//
// .npy format v1.0:
//   b'\x93NUMPY' + u8(major=1) + u8(minor=0) + u16_le(header_len) + header + data
//   header is a Python dict string, padded with spaces to a multiple of 64 bytes, ending with '\n'.

use std::io::{BufWriter, Write};
use std::fs::File;
use std::path::Path;

use anyhow::Result;
use num::complex::Complex;
use zip::write::{FileOptions, ZipWriter};
use zip::CompressionMethod;

type C128 = Complex<f64>;

// ---------------------------------------------------------------------------
// NPY header helpers
// ---------------------------------------------------------------------------

fn npy_header(dtype: &str, shape: &[u64]) -> Vec<u8> {
    let shape_str = match shape.len() {
        0 => "()".to_string(),
        1 => format!("({},)", shape[0]),
        _ => {
            let parts: Vec<String> = shape.iter().map(|x| x.to_string()).collect();
            format!("({})", parts.join(", "))
        }
    };
    // Magic + version (6 bytes) + header_len (2 bytes) = 8 bytes preamble
    // Total header area must be multiple of 64 bytes.
    // We'll format the dict, then pad to satisfy alignment.
    let dict_base = format!(
        "{{'descr': '{}', 'fortran_order': False, 'shape': {}, }}",
        dtype, shape_str
    );
    // preamble = 6 (magic) + 2 (version) + 2 (header_len field) = 10 bytes
    // header_area = len(dict) + 1 (\n) + padding spaces
    // total = 10 + header_area must be multiple of 64
    let preamble_len = 10usize;
    let dict_and_newline = dict_base.len() + 1;
    let total_min = preamble_len + dict_and_newline;
    let padding = (64 - total_min % 64) % 64;
    let header_str = format!("{}{}\n", dict_base, " ".repeat(padding));
    header_str.into_bytes()
}

fn write_npy_header(w: &mut impl Write, dtype: &str, shape: &[u64]) -> Result<()> {
    let header = npy_header(dtype, shape);
    w.write_all(b"\x93NUMPY")?;
    w.write_all(&[1u8, 0u8])?; // version 1.0
    let header_len = header.len() as u16;
    w.write_all(&header_len.to_le_bytes())?;
    w.write_all(&header)?;
    Ok(())
}

// ---------------------------------------------------------------------------
// Public write functions
// ---------------------------------------------------------------------------

/// Write complex128 data as .npy (dtype '<c16').
pub fn write_npy_c128(path: impl AsRef<Path>, data: &[C128], shape: &[u64]) -> Result<()> {
    let mut w = BufWriter::new(File::create(path.as_ref())?);
    write_npy_header(&mut w, "<c16", shape)?;
    for c in data {
        w.write_all(&c.re.to_le_bytes())?;
        w.write_all(&c.im.to_le_bytes())?;
    }
    Ok(())
}

/// Write f64 data as .npy (dtype '<f8').
pub fn write_npy_f64(path: impl AsRef<Path>, data: &[f64], shape: &[u64]) -> Result<()> {
    let mut w = BufWriter::new(File::create(path.as_ref())?);
    write_npy_header(&mut w, "<f8", shape)?;
    for &v in data {
        w.write_all(&v.to_le_bytes())?;
    }
    Ok(())
}

/// Write i64 data as .npy (dtype '<i8').
pub fn write_npy_i64(path: impl AsRef<Path>, data: &[i64], shape: &[u64]) -> Result<()> {
    let mut w = BufWriter::new(File::create(path.as_ref())?);
    write_npy_header(&mut w, "<i8", shape)?;
    for &v in data {
        w.write_all(&v.to_le_bytes())?;
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// NpzWriter: write multiple named arrays into a .npz (zip) file
// ---------------------------------------------------------------------------

pub struct NpzWriter {
    zip: ZipWriter<BufWriter<File>>,
}

impl NpzWriter {
    pub fn create(path: impl AsRef<Path>) -> Result<Self> {
        let file = BufWriter::new(File::create(path.as_ref())?);
        let zip = ZipWriter::new(file);
        Ok(Self { zip })
    }

    fn start_file(&mut self, name: &str) -> Result<()> {
        let options = FileOptions::default()
            .compression_method(CompressionMethod::Deflated);
        self.zip.start_file(format!("{}.npy", name), options)?;
        Ok(())
    }

    /// Write a complex128 array.
    pub fn write_c128(&mut self, name: &str, data: &[C128], shape: &[u64]) -> Result<()> {
        self.start_file(name)?;
        write_npy_header(&mut self.zip, "<c16", shape)?;
        for c in data {
            self.zip.write_all(&c.re.to_le_bytes())?;
            self.zip.write_all(&c.im.to_le_bytes())?;
        }
        Ok(())
    }

    /// Write a f64 array.
    pub fn write_f64(&mut self, name: &str, data: &[f64], shape: &[u64]) -> Result<()> {
        self.start_file(name)?;
        write_npy_header(&mut self.zip, "<f8", shape)?;
        for &v in data {
            self.zip.write_all(&v.to_le_bytes())?;
        }
        Ok(())
    }

    /// Write an i64 array.
    pub fn write_i64(&mut self, name: &str, data: &[i64], shape: &[u64]) -> Result<()> {
        self.start_file(name)?;
        write_npy_header(&mut self.zip, "<i8", shape)?;
        for &v in data {
            self.zip.write_all(&v.to_le_bytes())?;
        }
        Ok(())
    }

    /// Finalize the zip file.
    pub fn finish(mut self) -> Result<()> {
        self.zip.finish()?;
        Ok(())
    }
}
