use std::sync::Arc;

use parking_lot::{Mutex, MutexGuard};

use super::*;
use crate::PuzzleBuilder;

lua_userdata_value_conversion_wrapper! {
    #[name = "puzzlebuilder"]
    pub struct LuaPuzzleBuilder(Arc<Mutex<Option<PuzzleBuilder>>>);
}

impl LuaUserData for LuaNamedUserData<Arc<Mutex<Option<PuzzleBuilder>>>> {
    fn add_methods<'lua, T: LuaUserDataMethods<'lua, Self>>(methods: &mut T) {
        // TODO
    }
}

impl LuaPuzzleBuilder {
    pub fn lock(&self) -> MutexGuard<'_, Option<PuzzleBuilder>> {
        self.0.lock()
    }

    pub fn get(lua: LuaContext<'_>) -> LuaResult<Self> {
        lua.globals()
            .get("PUZZLE")
            .map_err(|_| LuaError::external("no puzzle being built"))
    }
    pub fn with<R>(
        lua: LuaContext<'_>,
        f: impl FnOnce(&mut PuzzleBuilder) -> LuaResult<R>,
    ) -> LuaResult<R> {
        let mutex = Self::get(lua)?;
        let mut mutex_guard = mutex.lock();
        let puzzle_builder = mutex_guard
            .as_mut()
            .ok_or_else(|| LuaError::external("no puzzle being built"))?;
        f(puzzle_builder)
    }
    pub fn take(lua: LuaContext<'_>) -> LuaResult<PuzzleBuilder> {
        Self::get(lua)?
            .lock()
            .take()
            .ok_or_else(|| LuaError::external("no puzzle being bulit"))
    }
}
