//! Mechanism for loading and sharing the simulation configuration

use crate::{evcut::EventCut, numeric::Real, Result};

use anyhow::{ensure, format_err, Context, Error};

use std::{fs::File, io::Read, str::FromStr};

/// This struct gives access to the simulation's configuration
pub struct Configuration {
    /// Number of events to be simulated
    pub num_events: usize,

    /// Collision energy at center of mass (GeV)
    pub e_tot: Real,

    /// Cuts on the angles and energies of generated photons
    pub event_cut: EventCut,

    /// Fine structure constant
    pub alpha: Real,

    /// Fine structure constant at the Z peak
    pub alpha_z: Real,

    /// Conversion factor from GeV^(-2) to pb
    pub convers: Real,

    /// Z⁰ boson mass (GeV)
    pub m_z0: Real,

    /// Z⁰ boson width (GeV)
    pub g_z0: Real,

    /// Square sine of Weinberg's Theta
    pub sin2_w: Real,

    /// Branching factor from Z to e+/e-
    pub br_ep_em: Real,

    /// Beta + (???)
    pub beta_plus: Real,

    /// Beta - (???)
    pub beta_minus: Real,

    /// Number of histogram bins (UNUSED)
    n_bin: i32,

    /// Whether intermediary results should be displayed (UNUSED)
    impr: bool,

    /// Whether results should be plotted in a histogram (UNUSED)
    plot: bool,
}
//
impl Configuration {
    /// Load the configuration from a file, check it, and print it out
    pub fn load(file_name: &str) -> Result<Self> {
        // Read out the simulation's configuration file or die trying.
        let config_str = {
            let mut config_file = File::open(file_name)?;
            let mut buffer = String::new();
            config_file.read_to_string(&mut buffer)?;
            buffer
        };

        // We will iterate over the configuration items. In 3photons' simple
        // config file format, these should be the first non-whitespace chunk of
        // text on each line. We will ignore blank lines.
        let mut config_iter = config_str
            .lines()
            .filter_map(|line| line.split_whitespace().next());

        // This closure fetches the next configuration item, tagging it with
        // the name of the configuration field which it is supposed to fill to
        // ease error reporting, and handling unexpected end-of-file too.
        let mut next_item = |name: &'static str| -> Result<ConfigItem> {
            config_iter
                .next()
                .map(|data| ConfigItem::new(name, data))
                .ok_or_else(|| format_err!("Missing configuration of {}", name))
        };

        // Decode the configuration items into concrete values
        let config = Configuration {
            num_events: next_item("num_events")?.parse::<usize>()?,
            e_tot: next_item("e_tot")?.parse::<Real>()?,
            event_cut: EventCut::new(
                next_item("a_cut")?.parse::<Real>()?,
                next_item("b_cut")?.parse::<Real>()?,
                next_item("e_min")?.parse::<Real>()?,
                next_item("sin_cut")?.parse::<Real>()?,
            ),
            alpha: next_item("alpha")?.parse::<Real>()?,
            alpha_z: next_item("alpha_z")?.parse::<Real>()?,
            convers: next_item("convers")?.parse::<Real>()?,
            m_z0: next_item("m_z0")?.parse::<Real>()?,
            g_z0: next_item("g_z0")?.parse::<Real>()?,
            sin2_w: next_item("sin2_w")?.parse::<Real>()?,
            br_ep_em: next_item("br_ep_em")?.parse::<Real>()?,
            beta_plus: next_item("beta_plus")?.parse::<Real>()?,
            beta_minus: next_item("beta_moins")?.parse::<Real>()?,
            n_bin: next_item("n_bin")?.parse::<i32>()?,
            impr: next_item("impr")?.parse_bool()?,
            plot: next_item("plot")?.parse_bool()?,
        };

        // Display it the way the C++ version used to (this eases comparisons)
        config.print();

        // A sensible simulation must run for at least one event
        ensure!(
            config.num_events > 0,
            "Invalid event count: {}",
            config.num_events
        );

        // NOTE: We don't support the original code's PAW-based plotting
        //       features, so we make sure that it was not enabled.
        ensure!(!config.plot, "Plotting is not supported by this version");

        // NOTE: We do not support the initial code's debugging feature which
        //       displays all intermediary results during sampling. Such a
        //       feature should be set up at build time to avoid run-time costs.
        ensure!(
            !config.impr,
            "Individual result printing is not supported. This debugging feature has a run-time \
             performance cost even when unused. It should be implemented at compile-time instead."
        );

        // If nothing bad occured, we can now return the configuration
        Ok(config)
    }

    /// Display the configuration, following formatting of the original version
    pub fn print(&self) {
        println!("ITOT           : {}", self.num_events);
        println!("ETOT           : {}", self.e_tot);
        println!("oCutpar.ACUT   : {}", self.event_cut.a_cut);
        println!("oCutpar.BCUT   : {}", self.event_cut.b_cut);
        println!("oCutpar.EMIN   : {}", self.event_cut.e_min);
        println!("oCutpar.SINCUT : {}", self.event_cut.sin_cut);
        println!("ALPHA          : {}", self.alpha);
        println!("ALPHAZ         : {}", self.alpha_z);
        println!("CONVERS        : {}", self.convers);
        println!("oParam.MZ0     : {}", self.m_z0);
        println!("oParam.GZ0     : {}", self.g_z0);
        println!("SIN2W          : {}", self.sin2_w);
        println!("BREPEM         : {}", self.br_ep_em);
        println!("BETAPLUS       : {}", self.beta_plus);
        println!("BETAMOINS      : {}", self.beta_minus);
        println!("NBIN           : {}", self.n_bin);
        println!("oParam.IMPR    : {}", self.impr);
        println!("PLOT           : {}", self.plot);
    }
}

/// A value from the configuration file, tagged with the struct field which it
/// is supposed to map for error reporting purposes.
struct ConfigItem<'a> {
    name: &'static str,
    data: &'a str,
}
//
impl<'a> ConfigItem<'a> {
    /// Build a config item from a struct field tag and raw iterator data
    fn new(name: &'static str, data: &'a str) -> Self {
        Self { name, data }
    }

    /// Parse this data using Rust's standard parsing logic
    fn parse<T: FromStr>(self) -> Result<T>
    where
        <T as FromStr>::Err: ::std::error::Error + Send + Sync + 'static,
    {
        self.data
            .parse::<T>()
            .map_err(Error::new)
            .context(format!("Could not parse configuration of {}", self.name))
            .map_err(|e| e.into())
    }

    /// Parse this data using special logic which handles Fortran's bool syntax
    fn parse_bool(self) -> Result<bool> {
        match self.data.to_lowercase().as_str() {
            // Handle FORTRAN booleans as a special case
            ".true." => Ok(true),
            ".false." => Ok(false),
            // Delegate other booleans to the standard Rust parser
            _ => self.parse::<bool>(),
        }
    }
}
