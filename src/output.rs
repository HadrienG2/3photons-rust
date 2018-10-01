//! This module is in charge of outputting the final simulation results to the
//! standard output and various files

use ::{
    config::Configuration,
    numeric::{
        functions::{abs, sqr, sqrt},
        Real,
    },
    rescont::{A, B_M, B_P, NUM_RESULTS, R_MX},
    resfin::{FinalResults, NUM_SPINS, SP_M, SP_P},
};

use chrono;

use std::{
  fs::{File, OpenOptions},
  io::{Result, Write},
  time::Duration,
};


// Output the simulation results to the console and to disk
pub fn dump_results(cfg: &Configuration,
                    res_fin: FinalResults,
                    elapsed_time: Duration) -> Result<()> {
    // Get current date and time for result archival purposes
    let current_time = chrono::Utc::now();
    let timestamp = current_time.format("%d-%b-%y   %T").to_string();

    // Write main results file. Try to mimick the original C++ format as well as
    // possible to ease comparisons, even where it makes little sense.
    {
        // Shorthands for spm2 and var are still convenient, but now we only
        // need read-only access to these quantities.
        let spm2 = &res_fin.spm2;
        let vars = &res_fin.vars;

        // Prepare to write our results into a file
        let mut res_file = File::create("res.dat")?;

        // TODO: Use a writer struct instead of closures

        // Create a few closure shorthands for common writing operations
        let write_i32 = |file: &mut File, title: &str, value: i32| {
          writeln!(*file, " {:<31}: {:<16}", title, value)
        };
        let write_usize = |file: &mut File, title: &str, value: usize| {
          writeln!(*file, " {:<31}: {:<16}", title, value)
        };
        let write_f = |file: &mut File, title: &str, value: Real| {
          // TODO: Use engineering notation here
          writeln!(*file, " {:<31}: {:<16}", title, value)
        };

        // Write the results to the file
        writeln!(res_file, " {}", timestamp.as_str())?;
        writeln!(res_file)?;
        write_usize(&mut res_file, "Nombre d'evenements", cfg.num_events)?;
        write_i32(&mut res_file, "... apres coupure", res_fin.selected_events)?;
        write_f(&mut res_file, "energie dans le CdM (GeV)", cfg.e_tot)?;
        write_f(&mut res_file,
                "coupure / cos(photon,faisceau)",
                cfg.event_cut.a_cut)?;
        write_f(&mut res_file,
                "coupure / cos(photon,photon)",
                cfg.event_cut.b_cut)?;
        write_f(&mut res_file,
                "coupure / sin(normale,faisceau)",
                cfg.event_cut.sin_cut)?;
        write_f(&mut res_file,
                "coupure sur l'energie (GeV)",
                cfg.event_cut.e_min)?;
        write_f(&mut res_file,
                "1/(constante de structure fine)",
                1./cfg.alpha)?;
        write_f(&mut res_file, "1/(structure fine au pic)", 1./cfg.alpha_z)?;
        write_f(&mut res_file, "facteur de conversion GeV-2/pb", cfg.convers)?;
        write_f(&mut res_file, "Masse du Z0 (GeV)", cfg.m_z0)?;
        write_f(&mut res_file, "Largeur du Z0 (GeV)", cfg.g_z0)?;
        write_f(&mut res_file, "Sinus^2 Theta Weinberg", cfg.sin2_w)?;
        write_f(&mut res_file, "Taux de branchement Z--->e+e-", cfg.br_ep_em)?;
        write_f(&mut res_file, "Beta plus", cfg.beta_plus)?;
        write_f(&mut res_file, "Beta moins", cfg.beta_minus)?;
        writeln!(res_file, " ---------------------------------------------")?;
        write_f(&mut res_file, "Section Efficace (pb)", res_fin.sigma)?;
        write_f(&mut res_file, "Ecart-Type (pb)", res_fin.sigma*res_fin.prec)?;
        write_f(&mut res_file, "Precision Relative", res_fin.prec)?;
        writeln!(res_file, " ---------------------------------------------")?;
        write_f(&mut res_file, "Beta minimum", res_fin.beta_min)?;
        write_f(&mut res_file, "Stat. Significance  B+(pb-1/2)", res_fin.ss_p)?;
        write_f(&mut res_file,
                "Incert. Stat. Sign. B+(pb-1/2)",
                res_fin.ss_p*res_fin.inc_ss_p)?;
        write_f(&mut res_file, "Stat. Significance  B-(pb-1/2)", res_fin.ss_m)?;
        write_f(&mut res_file,
                "Incert. Stat. Sign. B-(pb-1/2)",
                res_fin.ss_m*res_fin.inc_ss_m)?;

        // Write program performance stats
        let elapsed_secs = (elapsed_time.as_secs() as Real) +
                              1e-9 * (elapsed_time.subsec_nanos() as Real);
        writeln!(res_file, " ---------------------------------------------")?;
        writeln!(res_file, " Temps ecoule                   : ???")?;
        write_f(&mut res_file, "Temps ecoule utilisateur", elapsed_secs)?;
        writeln!(res_file, " Temps ecoule systeme           : ???")?;
        write_f(&mut res_file,
                "Temps ecoule par evenement",
                elapsed_secs / (cfg.num_events as Real))?;

        // Write more results (nature and purpose unclear in C++ code...)
        writeln!(res_file)?;
        for sp in 0..NUM_SPINS {
            for k in 0..NUM_RESULTS {
                writeln!(res_file,
                         "{:>3}{:>3}{:>15.7e}{:>15.7e}{:>15.7e}",
                         sp+1,
                         k+1,
                         spm2[(sp, k)],
                         spm2[(sp, k)]*vars[(sp, k)],
                         vars[(sp, k)])?;
            }
            writeln!(res_file)?;
        }
        for k in 0..NUM_RESULTS {
            let tmp1 = spm2[(SP_M, k)] + spm2[(SP_P, k)];
            let tmp2 = sqrt(sqr(spm2[(SP_M, k)]*vars[(SP_M, k)])
                              + sqr(spm2[(SP_P, k)]*vars[(SP_P, k)]));
            writeln!(res_file,
                     "   {:>3}{:>15.7e}{:>15.7e}{:>15.7e}",
                     k+1,
                     tmp1/4.,
                     tmp2/4.,
                     tmp2/abs(tmp1))?;
        }
    }

    // Print out some final results on stdout
    res_fin.eric();
    res_fin.fawzi();

    // Append the results of this run to a cumulative file
    //
    // NOTE: This part is completely broken in the C++ version, I did my best
    //       to fix it in this version.
    {
        assert_eq!(NUM_RESULTS, 5);
        let mut cum_res_file = OpenOptions::new().append(true)
                                                 .create(true)
                                                 .open("pil.mc")?;
        writeln!(cum_res_file, "{}", timestamp.as_str())?;
        let res1 = res_fin.spm2[(SP_M, A)] + res_fin.spm2[(SP_P, A)];
        let res2 = (res_fin.spm2[(SP_M, B_P)] + res_fin.spm2[(SP_P, B_P)]) *
                       sqr(cfg.beta_plus);
        let res3 = (res_fin.spm2[(SP_M, B_M)] + res_fin.spm2[(SP_P, B_M)]) *
                       sqr(cfg.beta_minus);
        let res4 = (res_fin.spm2[(SP_P, R_MX)] + res_fin.spm2[(SP_P, R_MX)]) *
                       cfg.beta_plus;
        writeln!(cum_res_file,
                 "{} {} {} {} {} {} {}",
                 cfg.e_tot,
                 res1/4.,
                 res2/4.,
                 res3/4.,
                 res4/4.,
                 (res1 + res2 + res3 + res4)/4.,
                 res_fin.sigma)?;
    }

    // ...and we're done
    Ok(())
}
