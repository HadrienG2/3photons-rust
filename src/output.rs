//! This module is in charge of outputting the final simulation results to the
//! standard output and various files

use crate::{
    config::Configuration,
    event::NUM_SPINS,
    numeric::{functions::*, reals, Real},
    rescont::{A, B_M, B_P, NUM_RESULTS, R_MX},
    resfin::{FinalResults, SP_M, SP_P},
};

use chrono;

use std::{
    fs::{File, OpenOptions},
    io::{Result, Write},
    time::Duration,
};

/// Write a floating-point number using "engineering" notation
///
/// Analogous to the %g format of the C printf function, this method switches
/// between naive and scientific notation for floating-point numbers when the
/// number being printed becomes so small that printing leading zeroes could end
/// up larger than the scientific notation, or so large that we would be forced
/// to print more significant digits than requested.
///
pub fn write_engineering(writer: &mut impl Write, x: Real, sig_digits: usize) -> Result<()> {
    let mut precision = sig_digits - 1;
    if x == 0. {
        // Zero is special because you can't take its log
        write!(writer, "0")
    } else {
        // Otherwise, use log to evaluate order of magnitude
        let log_x = x.abs().log10();
        if log_x >= -3. && log_x <= (sig_digits as Real) {
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

    // Number of significant digits in file output
    //
    // Must print one less than the actual machine type precision to match the
    // output of the C++ version of 3photons.
    //
    const SIG_DIGITS: usize = (reals::DIGITS - 1) as usize;

    // Create a few closure shorthands for common file writing operations
    let write_label = |file: &mut File, label: &str| write!(*file, " {:<31}: ", label);
    let write_usize = |file: &mut File, label: &str, value: usize| {
        write_label(file, label)?;
        writeln!(*file, "{}", value)
    };
    let write_real = |file: &mut File, label: &str, value: Real| {
        write_label(file, label)?;
        write_engineering(file, value, SIG_DIGITS)?;
        writeln!(file)
    };

    // Compute a timestamp of when the run ended
    let current_time = chrono::Utc::now();
    let timestamp = current_time.format("%d-%b-%y   %T");

    // Write execution timings to a file
    {
        // Prepare to write our timings into a file
        let mut tim_file = File::create("res.times")?;

        // Write a timestamp of when the run ended
        writeln!(tim_file, " {}", timestamp)?;

        // Write program performance stats
        let elapsed_secs =
            (elapsed_time.as_secs() as Real) + 1e-9 * (elapsed_time.subsec_nanos() as Real);
        writeln!(tim_file, " ---------------------------------------------")?;
        writeln!(tim_file, " Temps ecoule                   : ???")?;
        write_real(&mut tim_file, "Temps ecoule utilisateur", elapsed_secs)?;
        writeln!(tim_file, " Temps ecoule systeme           : ???")?;
        write_real(
            &mut tim_file,
            "Temps ecoule par evenement",
            elapsed_secs / (cfg.num_events as Real),
        )?;
    }

    // Write main results file. Try to mimick the original C++ format as well as
    // possible to ease comparisons, even where it makes little sense.
    {
        // Some convenience shorthand
        let spm2 = &res_fin.spm2;
        let vars = &res_fin.vars;

        // Prepare to write our results into a file
        let mut dat_file = File::create("res.data")?;

        // Write the results to the file
        write_usize(&mut dat_file, "Nombre d'evenements", cfg.num_events)?;
        write_usize(&mut dat_file, "... apres coupure", res_fin.selected_events)?;
        write_real(&mut dat_file, "energie dans le CdM      (GeV)", cfg.e_tot)?;
        write_real(
            &mut dat_file,
            "coupure / cos(photon,faisceau)",
            cfg.event_cut.a_cut,
        )?;
        write_real(
            &mut dat_file,
            "coupure / cos(photon,photon)",
            cfg.event_cut.b_cut,
        )?;
        write_real(
            &mut dat_file,
            "coupure / sin(normale,faisceau)",
            cfg.event_cut.sin_cut,
        )?;
        write_real(
            &mut dat_file,
            "coupure sur l'energie    (GeV)",
            cfg.event_cut.e_min,
        )?;
        write_real(
            &mut dat_file,
            "1/(constante de structure fine)",
            1. / cfg.alpha,
        )?;
        write_real(&mut dat_file, "1/(structure fine au pic)", 1. / cfg.alpha_z)?;
        write_real(&mut dat_file, "facteur de conversion GeV-2/pb", cfg.convers)?;
        write_real(&mut dat_file, "Masse du Z0              (GeV)", cfg.m_z0)?;
        write_real(&mut dat_file, "Largeur du Z0            (GeV)", cfg.g_z0)?;
        write_real(&mut dat_file, "Sinus^2 Theta Weinberg", cfg.sin2_w)?;
        write_real(&mut dat_file, "Taux de branchement Z--->e+e-", cfg.br_ep_em)?;
        write_real(&mut dat_file, "Beta plus", cfg.beta_plus)?;
        write_real(&mut dat_file, "Beta moins", cfg.beta_minus)?;
        writeln!(dat_file, " ---------------------------------------------")?;
        write_real(
            &mut dat_file,
            "Section Efficace          (pb)",
            res_fin.sigma,
        )?;
        write_real(
            &mut dat_file,
            "Ecart-Type                (pb)",
            res_fin.sigma * res_fin.prec,
        )?;
        write_real(&mut dat_file, "Precision Relative", res_fin.prec)?;
        writeln!(dat_file, " ---------------------------------------------")?;
        write_real(&mut dat_file, "Beta minimum", res_fin.beta_min)?;
        write_real(
            &mut dat_file,
            "Stat. Significance  B+(pb-1/2)",
            res_fin.ss_p,
        )?;
        write_real(
            &mut dat_file,
            "Incert. Stat. Sign. B+(pb-1/2)",
            res_fin.ss_p * res_fin.inc_ss_p,
        )?;
        write_real(
            &mut dat_file,
            "Stat. Significance  B-(pb-1/2)",
            res_fin.ss_m,
        )?;
        write_real(
            &mut dat_file,
            "Incert. Stat. Sign. B-(pb-1/2)",
            res_fin.ss_m * res_fin.inc_ss_m,
        )?;

        // Write more results (nature and purpose unclear in C++ code...)
        writeln!(dat_file)?;
        let decimals = (SIG_DIGITS - 1).min(7);
        for sp in 0..NUM_SPINS {
            for k in 0..NUM_RESULTS {
                writeln!(
                    dat_file,
                    "{:>3}{:>3}{:>width$.decs$e}{:>width$.decs$e}{:>width$.decs$e}",
                    sp + 1,
                    k + 1,
                    spm2[(sp, k)],
                    spm2[(sp, k)].abs() * vars[(sp, k)],
                    vars[(sp, k)],
                    width = decimals + 8,
                    decs = decimals,
                )?;
            }
            writeln!(dat_file)?;
        }
        for k in 0..NUM_RESULTS {
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
        assert_eq!(NUM_RESULTS, 5);
        let mut cum_dat_file = OpenOptions::new()
            .append(true)
            .create(true)
            .open("pil.mc")?;
        writeln!(cum_dat_file, "{}", timestamp)?;
        let res1 = res_fin.spm2[(SP_M, A)] + res_fin.spm2[(SP_P, A)];
        let res2 = (res_fin.spm2[(SP_M, B_P)] + res_fin.spm2[(SP_P, B_P)]) * sqr(cfg.beta_plus);
        let res3 = (res_fin.spm2[(SP_M, B_M)] + res_fin.spm2[(SP_P, B_M)]) * sqr(cfg.beta_minus);
        let res4 = (res_fin.spm2[(SP_P, R_MX)] + res_fin.spm2[(SP_P, R_MX)]) * cfg.beta_plus;
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
