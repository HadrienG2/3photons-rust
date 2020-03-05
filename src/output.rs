//! This module is in charge of outputting the final simulation results to the
//! standard output and various files

use crate::{
    config::Configuration,
    event::NUM_SPINS,
    matelems::{A, B_M, B_P, NUM_MAT_ELEMS, R_MX},
    numeric::{functions::*, reals, Float},
    resfin::{FinalResults, SP_M, SP_P},
};

use chrono;

use std::{
    fs::{File, OpenOptions},
    io::{Result, Write},
    time::Duration,
};

// Number of significant digits in file output
//
// Must print one less than the actual machine type precision to match the
// output of the C++ version of 3photons.
//
const SIG_DIGITS: usize = (reals::DIGITS - 1) as usize;

/// Output the simulation results to the console and to disk
#[allow(clippy::cast_lossless)]
pub fn dump_results(
    cfg: &Configuration,
    res_fin: &FinalResults,
    elapsed_time: Duration,
) -> Result<()> {
    // Print out some final results on stdout
    res_fin.eric();
    res_fin.fawzi();

    // Compute a timestamp of when the run ended
    let current_time = chrono::Utc::now();
    let timestamp = current_time.format("%d-%b-%y   %T").to_string();

    // Write execution timings to a file
    {
        // Prepare to write our timings into a file
        let mut tim_file = File::create("res.times")?;
        let tim_file = &mut tim_file;

        // Write a timestamp of when the run ended
        writeln_3p(tim_file, &timestamp[..])?;

        // Write program performance stats
        let elapsed_secs =
            (elapsed_time.as_secs() as Float) + 1e-9 * (elapsed_time.subsec_nanos() as Float);
        writeln_3p(tim_file, "---------------------------------------------")?;
        writeln_3p(tim_file, ("Temps ecoule", "???"))?;
        writeln_3p(tim_file, ("Temps ecoule utilisateur", elapsed_secs))?;
        writeln_3p(tim_file, ("Temps ecoule systeme", "???"))?;
        let secs_per_ev = elapsed_secs / (cfg.num_events as Float);
        writeln_3p(tim_file, ("Temps ecoule par evenement", secs_per_ev))?;
    }

    // Write main results file. Try to mimick the original C++ format as well as
    // possible to ease comparisons, even where it makes little sense.
    {
        // Shorthands
        let ev_cut = &cfg.event_cut;

        // Prepare to write our results into a file
        let mut dat_file = File::create("res.data")?;
        let dat_file = &mut dat_file;

        // Write the results to the file
        writeln_3p(dat_file, ("Nombre d'evenements", cfg.num_events))?;
        writeln_3p(dat_file, ("... apres coupure", res_fin.selected_events))?;
        writeln_3p(dat_file, ("energie dans le CdM      (GeV)", cfg.e_tot))?;
        writeln_3p(dat_file, ("coupure / cos(photon,faisceau)", ev_cut.a_cut))?;
        writeln_3p(dat_file, ("coupure / cos(photon,photon)", ev_cut.b_cut))?;
        let sin_cut = ev_cut.sin_cut;
        writeln_3p(dat_file, ("coupure / sin(normale,faisceau)", sin_cut))?;
        writeln_3p(dat_file, ("coupure sur l'energie    (GeV)", ev_cut.e_min))?;
        let inv_alpha = 1. / cfg.alpha;
        writeln_3p(dat_file, ("1/(constante de structure fine)", inv_alpha))?;
        writeln_3p(dat_file, ("1/(structure fine au pic)", 1. / cfg.alpha_z))?;
        writeln_3p(dat_file, ("facteur de conversion GeV-2/pb", cfg.convers))?;
        writeln_3p(dat_file, ("Masse du Z0              (GeV)", cfg.m_z0))?;
        writeln_3p(dat_file, ("Largeur du Z0            (GeV)", cfg.g_z0))?;
        writeln_3p(dat_file, ("Sinus^2 Theta Weinberg", cfg.sin2_w))?;
        writeln_3p(dat_file, ("Taux de branchement Z--->e+e-", cfg.br_ep_em))?;
        writeln_3p(dat_file, ("Beta plus", cfg.beta_plus))?;
        writeln_3p(dat_file, ("Beta moins", cfg.beta_minus))?;
        writeln_3p(dat_file, "---------------------------------------------")?;
        writeln_3p(dat_file, ("Section Efficace          (pb)", res_fin.sigma))?;
        let stddev_res = res_fin.sigma * res_fin.prec;
        writeln_3p(dat_file, ("Ecart-Type                (pb)", stddev_res))?;
        writeln_3p(dat_file, ("Precision Relative", res_fin.prec))?;
        writeln_3p(dat_file, "---------------------------------------------")?;
        writeln_3p(dat_file, ("Beta minimum", res_fin.beta_min))?;
        writeln_3p(dat_file, ("Stat. Significance  B+(pb-1/2)", res_fin.ss_p))?;
        let incert_ss_p = res_fin.ss_p * res_fin.inc_ss_p;
        writeln_3p(dat_file, ("Incert. Stat. Sign. B+(pb-1/2)", incert_ss_p))?;
        writeln_3p(dat_file, ("Stat. Significance  B-(pb-1/2)", res_fin.ss_m))?;
        let incert_ss_m = res_fin.ss_m * res_fin.inc_ss_m;
        writeln_3p(dat_file, ("Incert. Stat. Sign. B-(pb-1/2)", incert_ss_m))?;

        // Write more results (nature and purpose unclear in C++ code...)
        writeln!(dat_file)?;
        let decimals = (SIG_DIGITS - 1).min(7);
        for sp in 0..NUM_SPINS {
            for k in 0..NUM_MAT_ELEMS {
                writeln!(
                    dat_file,
                    "{:>3}{:>3}{:>width$.decs$e}{:>width$.decs$e}{:>width$.decs$e}",
                    sp + 1,
                    k + 1,
                    res_fin.spm2[(sp, k)],
                    res_fin.spm2[(sp, k)].abs() * res_fin.vars[(sp, k)],
                    res_fin.vars[(sp, k)],
                    width = decimals + 8,
                    decs = decimals,
                )?;
            }
            writeln!(dat_file)?;
        }
        for k in 0..NUM_MAT_ELEMS {
            // FIXME: Vectorize sums across spins once const generics enable
            //        more ergonomic small matrix manipulations.
            let spm2 = &res_fin.spm2;
            let vars = &res_fin.vars;
            let tmp1 = spm2[(SP_M, k)] + spm2[(SP_P, k)];
            let tmp2 = sqrt(
                sqr(spm2[(SP_M, k)] * vars[(SP_M, k)]) + sqr(spm2[(SP_P, k)] * vars[(SP_P, k)]),
            );
            writeln!(
                dat_file,
                "   {:>3}{:>width$.decs$e}{:>width$.decs$e}{:>width$.decs$e}",
                k + 1,
                tmp1 / 4.,
                tmp2 / 4.,
                tmp2 / abs(tmp1),
                width = decimals + 8,
                decs = decimals,
            )?;
        }
    }

    // Append the results of this run to a cumulative file
    //
    // NOTE: This part is completely broken in the C++ version, I did my best
    //       to fix it in this version.
    {
        assert_eq!(NUM_MAT_ELEMS, 5);

        let mut cum_dat_file = OpenOptions::new()
            .append(true)
            .create(true)
            .open("pil.mc")?;

        writeln!(cum_dat_file, "{}", timestamp)?;

        // FIXME: Vectorize sums across spins once const generics enable more
        //        ergonomic small matrix manipulations.
        let spm2 = &res_fin.spm2;
        let res1 = spm2[(SP_M, A)] + spm2[(SP_P, A)];
        let res2 = (spm2[(SP_M, B_P)] + spm2[(SP_P, B_P)]) * sqr(cfg.beta_plus);
        let res3 = (spm2[(SP_M, B_M)] + spm2[(SP_P, B_M)]) * sqr(cfg.beta_minus);
        let res4 = (spm2[(SP_P, R_MX)] + spm2[(SP_P, R_MX)]) * cfg.beta_plus;
        writeln!(
            cum_dat_file,
            "{} {} {} {} {} {} {}",
            cfg.e_tot,
            res1 / 4.,
            res2 / 4.,
            res3 / 4.,
            res4 / 4.,
            (res1 + res2 + res3 + res4) / 4.,
            res_fin.sigma
        )?;
    }

    // ...and we're done
    Ok(())
}

/// Text output facility that mimicks 3photons' file output styling
fn writeln_3p(file: &mut File, data: impl Write3p) -> Result<()> {
    write!(file, " ")?;
    data.write(file)?;
    writeln!(file)
}

/// Trait implemented by things which can be printed The 3photons Way (tm)
trait Write3p: Sized {
    /// Write down `self` to the output file using 3photons style
    fn write(self, file: &mut File) -> Result<()>;
}

impl Write3p for &str {
    // Strings work in the usual way
    fn write(self, file: &mut File) -> Result<()> {
        write!(file, "{}", self)
    }
}

impl Write3p for usize {
    // Integers work in the usual way too
    // FIXME: Simplify and generalize this once Rust has specialization
    fn write(self, file: &mut File) -> Result<()> {
        write!(file, "{}", self)
    }
}

impl Write3p for Float {
    // 3photons used %g for floats, this is a close approximation
    fn write(self, file: &mut File) -> Result<()> {
        write_engineering(file, self, SIG_DIGITS)
    }
}

impl<T: Write3p> Write3p for (&str, T) {
    // Key-value output that uses fixed-size columns for better readability
    fn write(self, file: &mut File) -> Result<()> {
        write!(*file, "{:<31}: ", self.0)?;
        self.1.write(file)
    }
}

/// Write a floating-point number using "engineering" notation
///
/// Analogous to the %g format of the C printf function, this method switches
/// between naive and scientific notation for floating-point numbers when the
/// number being printed becomes so small that printing leading zeroes could end
/// up larger than the scientific notation, or so large that we would be forced
/// to print more significant digits than requested.
///
fn write_engineering(writer: &mut impl Write, x: Float, sig_digits: usize) -> Result<()> {
    let mut precision = sig_digits - 1;
    if x == 0. {
        // Zero is special because you can't take its log
        write!(writer, "0")
    } else {
        // Otherwise, use log to evaluate order of magnitude
        let log_x = x.abs().log10();
        if log_x >= -3. && log_x <= (sig_digits as Float) {
            // Print using naive notation
            //
            // Since Rust's precision controls number of digits after the
            // decimal point, we must adjust it depending on magnitude in order
            // to operate at a constant number of significant digits.
            precision = (precision as isize - log_x.trunc() as isize) as usize;

            // Numbers smaller than 1 must get one extra digit since the leading
            // zero does not count as a significant digit.
            if log_x < 0. {
                precision += 1
            }

            // People don't normally expect trailing zeros or decimal point in
            // naive notation, but be careful with integer numbers...
            let str_with_zeros = format!("{:.1$}", x, precision);
            if str_with_zeros.contains('.') {
                write!(
                    writer,
                    "{}",
                    str_with_zeros.trim_end_matches('0').trim_end_matches('.')
                )
            } else {
                write!(writer, "{}", str_with_zeros)
            }
        } else {
            // Print using scientific notation
            write!(writer, "{:.1$e}", x, precision)
        }
    }
}
