use std::env;
use std::path::PathBuf;

use directories::ProjectDirs;

const PREFS_FILE_NAME: &str = "hyperspeedcube";
const PREFS_FILE_EXTENSION: &str = "yaml";

const LUA_DIR_NAME: &str = "lua";

lazy_static! {
    pub(crate) static ref PATHS: Option<AppPaths> = app_paths();
}

/// Paths to external files read by Hyperspeedcube.
pub struct AppPaths {
    /// Path to the Hyperspeedcube preferences file.
    pub prefs_file: PathBuf,
    /// Path to the Hyperspeedcube Lua directory.
    pub lua_dir: PathBuf,
}

/// Returns the app paths.
///
/// - For dev builds, uses the project directory.
/// - For official release builds in portable mode (the default on Windows &
///   Linux), uses the directory of the current executable.
/// - For official release builds in nonportable mode (the default on macOS),
///   uses the system directories.
///
/// If the preferred behavior (portable vs. nonportable) fails, then this
/// function falls back on the other.
fn app_paths() -> Option<AppPaths> {
    match is_nonportable() {
        true => nonportable_paths().or_else(portable_paths),
        false => portable_paths().or_else(nonportable_paths),
    }
}

fn nonportable_paths() -> Option<AppPaths> {
    match ProjectDirs::from("", "", "Hyperspeedcube") {
        Some(dirs) => {
            log::info!("Using nonportable paths");
            Some(AppPaths {
                prefs_file: dirs
                    .config_dir()
                    .join(format!("{PREFS_FILE_NAME}.{PREFS_FILE_EXTENSION}")),
                lua_dir: dirs.data_dir().join(LUA_DIR_NAME),
            })
        }
        None => {
            log::error!("Error getting nonportable directories");
            None
        }
    }
}

fn portable_paths() -> Option<AppPaths> {
    match portable_dir() {
        Some(dir) => {
            log::info!("Using portable paths");
            Some(AppPaths {
                prefs_file: dir.join(format!("{PREFS_FILE_NAME}.{PREFS_FILE_EXTENSION}")),
                lua_dir: dir.join(LUA_DIR_NAME),
            })
        }
        None => {
            log::error!("Error getting portable directory");
            None
        }
    }
}

fn portable_dir() -> Option<PathBuf> {
    let exe_path = env::current_exe().ok()?.canonicalize().ok()?;
    if crate::IS_OFFICIAL_BUILD {
        // `/hyperspeedcube.exe`
        Some(exe_path.parent()?.to_path_buf())
    } else {
        // `/target/debug/hyperspeedcube.exe`
        Some(exe_path.parent()?.parent()?.parent()?.to_path_buf())
    }
}

fn is_nonportable() -> bool {
    if crate::IS_OFFICIAL_BUILD && cfg!(target_os = "macos") {
        // If we are in a macOS app package, then we are always nonportable
        // because macOS doesn't allow storing files in the same directory
        // as the executable.
        true
    } else if let Some(mut p) = portable_dir() {
        // Otherwise, check whether the `nonportable` file exists in the same
        // directory as the executable.
        p.push("nonportable");
        p.exists()
    } else {
        false
    }
}