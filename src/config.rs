//! Mechanism for loading and sharing the simulation configuration

use evcut::EventCut;
use numeric::Real;
use std::fs::File;
use std::io::Read;


// This struct loads and gives access to the simulation's configuration
pub struct Configuration {
    /// Number of events to be simulated
    pub num_events: i32,

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
    // Load the configuration from a file, check it, and print it out
    pub fn new(file_name: &str) -> Result<Self> {
        // ### LOAD CONFIGURATION ###

        // Read out the simulation configuration file or die trying.
        let mut config_file = File::open(file_name)?;
        let mut config_str = String::new();
        config_file.read_to_string(&mut config_str)?;

        // We will iterate over the configuration items. In 3photons' simple
        // config file format, these should be the first non-whitespace chunk of
        // text on each line. We will ignore blank lines.
        let mut config_iter =
            config_str.lines()
                      .filter_map(|line| line.split_whitespace()
                                             .next());

        // To allow better error reporting, configuration items will be tagged
        // with the name of the variable that they are mapped to.
        struct ConfigItem<'a> {
            name: &'static str,
            value: Result<&'a str>,
        };

        // Fetch the next configuration item, in textual form
        let mut next_item = |name: &'static str| -> ConfigItem {
            ConfigItem {
                name: name,
                value: config_iter.next()
                                  .ok_or(ErrorKind::Missing(name).into())
            }
        };

        // Now that we have textual configuration items, we can parse them.
        // The process is mostly the same for integers and reals, in fact we
        // could have used a single generic implementation here if Rust allowed
        // us to write generic closures.
        let as_i32 = |item: ConfigItem| -> Result<i32> {
            let item_name = item.name;
            item.value?
                .parse::<i32>()
                .map_err(|e| ErrorKind::Unreadable(item_name,
                                                   Box::new(e)).into())
        };
        let as_real = |item: ConfigItem| -> Result<Real> {
            let item_name = item.name;
            item.value?
                .parse::<Real>()
                .map_err(|e| ErrorKind::Unreadable(item_name,
                                                   Box::new(e)).into())
        };

        // Booleans are a bit special: for compatibility with the original
        // 3photons code, we also support FORTRAN syntax for them.
        let as_bool = |item: ConfigItem| -> Result<bool> {
            let item_name = item.name;
            let item_value_lower = item.value?.to_lowercase();
            match item_value_lower.as_str() {
                // Handle FORTRAN booleans as a special case
                ".true." => Ok(true),
                ".false." => Ok(false),
                // Delegate other booleans to the standard Rust parser
                other =>
                    other.parse::<bool>()
                         .map_err(|e| ErrorKind::Unreadable(item_name,
                                                            Box::new(e)).into())
            }
        };

        // Decode the configuration items into concrete values
        let config = Configuration {
            num_events: as_i32(next_item("num_events"))?,
            e_tot: as_real(next_item("e_tot"))?,
            event_cut: EventCut::new(as_real(next_item("a_cut"))?,
                                     as_real(next_item("b_cut"))?,
                                     as_real(next_item("e_min"))?,
                                     as_real(next_item("sin_cut"))?),
            alpha: as_real(next_item("alpha"))?,
            alpha_z: as_real(next_item("alpha_z"))?,
            convers: as_real(next_item("convers"))?,
            m_z0: as_real(next_item("m_z0"))?,
            g_z0: as_real(next_item("g_z0"))?,
            sin2_w: as_real(next_item("sin2_w"))?,
            br_ep_em: as_real(next_item("br_ep_em"))?,
            beta_plus: as_real(next_item("beta_plus"))?,
            beta_minus: as_real(next_item("beta_moins"))?,
            n_bin: as_i32(next_item("n_bin"))?,
            impr: as_bool(next_item("impr"))?,
            plot: as_bool(next_item("plot"))?,
        };


        // ### DISPLAY IT ###

        // Print out the configuration that we extracted from the file, in the
        // same format as used by the original C++ version (this eases comparisons)
        println!("ITOT           : {}", config.num_events);
        println!("ETOT           : {}", config.e_tot);
        println!("oCutpar.ACUT   : {}", config.event_cut.a_cut);
        println!("oCutpar.BCUT   : {}", config.event_cut.b_cut);
        println!("oCutpar.EMIN   : {}", config.event_cut.e_min);
        println!("oCutpar.SINCUT : {}", config.event_cut.sin_cut);
        println!("ALPHA          : {}", config.alpha);
        println!("ALPHAZ         : {}", config.alpha_z);
        println!("CONVERS        : {}", config.convers);
        println!("oParam.MZ0     : {}", config.m_z0);
        println!("oParam.GZ0     : {}", config.g_z0);
        println!("SIN2W          : {}", config.sin2_w);
        println!("BREPEM         : {}", config.br_ep_em);
        println!("BETAPLUS       : {}", config.beta_plus);
        println!("BETAMOINS      : {}", config.beta_minus);
        println!("NBIN           : {}", config.n_bin);
        println!("oParam.IMPR    : {}", config.impr);
        println!("PLOT           : {}", config.plot);


        // ### CHECK UNSUPPORTED FEATURES ###

        // NOTE: This is where the FORTRAN code would setup PAW for plotting.
        //       We don't support plotting, so we only check that it's disabled.
        if config.plot { return Err(ErrorKind::Unsupported("plot").into()); }

        // NOTE: We do not support the initial code's debugging feature which
        //       displays all intermediary results during sampling. That feature
        //       should be configured at compile time to avoid run-time costs.
        if config.impr { return Err(ErrorKind::Unsupported("impr").into()); }

        // If nothing bad occured, we can now return the configuration
        Ok(config)
    }
}


// Here are the various things that can go wrong while loading the configuration
mod config_errors {
    error_chain!{
        foreign_links{
            // Failed to load the simulation configuration
            Io(::std::io::Error);
        }

        errors {
            // The configuration file is missing some field
            Missing(field: &'static str) {
                description("Missing configuration field")
                display("Missing configuration field: '{}'", field)
            }

            // A configuration field could not be parsed
            Unreadable(field: &'static str,
                       parse_error: Box<::std::error::Error + Send>) {
                description("Failed to parse configuration field")
                display("Failed to parse configuration for field '{}' ({})",
                        field,
                        parse_error)
            }

            // The value of a configuration field is valid, but unsupported
            Unsupported(field: &'static str) {
                description("Unsupported configuration")
                display("Unsupported configuration for field '{}'", field)
            }
        }
    }
}
//
pub use self::config_errors::*;
